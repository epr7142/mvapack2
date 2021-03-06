## Copyright (C) 2013 University of Nebraska Board of Regents.
## Written by Bradley Worley <bradley.worley@huskers.unl.edu>.
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; If not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @anchor{permtest}
## @deftypefn {Function File} {@var{S} =} permtest (@var{mdl})
## @deftypefnx {Function File} {@var{S} =} permtest (@var{mdl}, @var{n})
## Uses response permutation testing of a PLS or OPLS model to calculate
## the significance of the model parameters. The output structure @var{S}
## contains the following information:
##
## @code{S.n}: number of permutations performed. @*
## @code{S.r}: permutation Y-correlation coefficients. @*
## @code{S.Rsq.orig}: original @math{R^2} value. @*
## @code{S.Qsq.orig}: original @math{Q^2} value. @*
## @code{S.Rsq.perm}: permutation @math{R^2} values. @*
## @code{S.Qsq.perm}: permutation @math{Q^2} values. @*
## @code{S.Rsq.t}: @math{R^2} t-statistic. @*
## @code{S.Qsq.t}: @math{Q^2} t-statistic. @*
## @code{S.Rsq.p}: @math{R^2} t-test p-value. @*
## @code{S.Qsq.p}: @math{Q^2} t-test p-value. @*
##
## An optional second argument @var{n} may be passed to set the number of
## permutations to execute. The default number of permutations is 100.
##
## The structure @var{S} generated by this function may be passed to
## @ref{permscatter} or @ref{permdensity} for visualization of the
## permutation test results.
## @end deftypefn

function S = permtest (mdl, n, doperm)
  % check for proper arguments.
  if (!any(nargin == [1 : 2]) || nargout != 1 || !isstruct(mdl))
    % improper arguments. print the usage statement.
    print_usage();
  end

  % check that the model is a pls or opls model.
  if (!strcmp(mdl.type, 'pls') && ...
      !strcmp(mdl.type, 'opls') && ...
      !strcmp(mdl.type, 'lda'))
    % invalid model type. throw an exception.
    error('permtest: model type must be supervised');
  end

  % check if a second argument was supplied.
  if (nargin < 2 || isempty(n))
    % use the default number of permutations: 100.
    n = 100;
  else
    % ensure the argument is a real scalar.
    if (!isscalar(n) || !isreal(n))
      % throw an exception.
      error('permtest: number of permutations must be a real scalar');
    end
  end

  % check if a third argument was supplied.
  if (nargin < 3 || isempty(doperm))
    % use the default value.
    doperm = true;
  end

  % check that the model has a builder function handle.
  if (!isfield(mdl, 'builder') || !is_function_handle(mdl.builder))
    % model has no builder function handle. throw an exception.
    error('permtest: model lacks information on how to rebuild itself');
  end

  % extract the original data matrix from the model.
  if (ismultiblock(mdl))
    % extract the multiblock array.
    X = mdl.Xmb;
  else
    % extract the data matrix.
    X = mdl.X0;
  end

  % extract the original response matrix from the model.
  Y = mdl.Y0;

  % extract the component count from the model.
  if (strcmp(mdl.type, 'opls'))
    % the component count will contain information about how to exactly
    % rebuild the desired combination of predictive and orthogonal
    % components.
    Aout = [mdl.Ap; mdl.no];
  else
    % the component count will be a simple scalar.
    Aout = mdl.A;
  end

  % initialize the validation arrays.
  S.Rsq.orig = mdl.Rsq.Y.cum;
  S.Qsq.orig = mdl.Qsq.cum;
  S.Rsq.perm = zeros(n, 1);
  S.Qsq.perm = zeros(n, 1);
  S.r = zeros(n, 1);
  S.n = n;

  % check if the model supports CV-ANOVA.
  if (strcmp(mdl.type, 'pls') || strcmp(mdl.type, 'opls'))
    % run the test.
    CVA = cvanova(mdl);

    % initialize the extra validation arrays.
    S.CVA.F.orig = CVA.F;
    S.CVA.p.orig = CVA.p;
    S.CVA.F.perm = zeros(n, 1);
    S.CVA.p.perm = zeros(n, 1);
  end

  % run through the permutations.
  for idx = 1 : n
    % shuffle the response matrix.
    if (doperm)
      Yshuf = shuffle(Y);
    else
      Yshuf = Y;
    end

    % build the model.
    mdlPerm = mdl.builder(X, Yshuf, mdl.scaling, [], Aout);

    % calculate the correlation between the original and shuffled response.
    S.r(idx) = mean(diag(abs(corr(Y, Yshuf))));

    % add the cross-validation statistics to the lot.
    S.Rsq.perm(idx) = mdlPerm.Rsq.Y.cum;
    S.Qsq.perm(idx) = mdlPerm.Qsq.cum;

    % check if the model supports CV-ANOVA.
    if (strcmp(mdl.type, 'pls') || strcmp(mdl.type, 'opls'))
      % run the test.
      CVA = cvanova(mdlPerm);

      % add the results to the lot.
      S.CVA.F.perm(idx) = CVA.F;
      S.CVA.p.perm(idx) = CVA.p;
    end

    % clear the permuted model.
    clear mdlPerm;
  end

  % non-parametrically calculate the p-values.
  S.Rsq.p = sum(S.Rsq.perm >= S.Rsq.orig) / n;
  S.Qsq.p = sum(S.Qsq.perm >= S.Qsq.orig) / n;

  % calculate the slopes of the fit lines.
  S.Rsq.m = [S.r - 1] \ (S.Rsq.perm - S.Rsq.orig);
  S.Qsq.m = [S.r - 1] \ (S.Qsq.perm - S.Qsq.orig);

  % calculate the two points of the fit lines.
  S.Rsq.xy = [[0; 1], [S.Rsq.orig - S.Rsq.m; S.Rsq.orig]];
  S.Qsq.xy = [[0; 1], [S.Qsq.orig - S.Qsq.m; S.Qsq.orig]];

  % store the number of model components.
  S.A = mdl.A;
end

