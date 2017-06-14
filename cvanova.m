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
## @anchor{cvanova}
## @deftypefn {Function File} {@var{S} =} cvanova (@var{mdl})
## Uses an already performed internal cross-validation of a PLS or
## OPLS model to calculate a CV-ANOVA p-value. The output structure @var{S}
## contains the following information:
##
## @code{S.SS}: cross-validated sum of squares. @*
## @code{S.MS}: cross-validated mean square errors. @*
## @code{S.DF}: degrees of freedom. @*
## @code{S.F}: f-statistic from cross-validated mean square errors. @*
## @code{S.p}: p-value indicating how well the model fits the data. @*
##
## @quotation
## L. Eriksson, J. Trygg, S. Wold. `CV-ANOVA for significance testing of
## of PLS and OPLS models.' J. Chemometrics 2008(22): 594-600.
## @end quotation
## @end deftypefn

function S = cvanova (mdl)
  % check for proper arguments.
  if (nargin != 1 || !isstruct(mdl))
    % improper arguments. print the usage statement.
    print_usage();
  end

  % check that the model is a pls or opls model.
  if (!strcmp(mdl.type, 'pls') && !strcmp(mdl.type, 'opls'))
    % invalid model type. throw an exception.
    error('cvanova: model type must be either pls or opls');
  end

  % check that the model has a cross-validation structure.
  if (!isfield(mdl, 'cv') || !isfield(mdl, 'CV') || isempty(mdl.CV))
    % throw an error.
    error('cvanova: model has no cross-validation information');
  end

  % initialize the total, fit and residual sum of squares matrices.
  SStot = zeros(mdl.cv.niter, 1);
  SSfit = zeros(mdl.cv.niter, 1);
  SSerr = zeros(mdl.cv.niter, 1);

  % initialize the total, fit and residual mean square error matrices.
  MStot = zeros(mdl.cv.niter, 1);
  MSfit = zeros(mdl.cv.niter, 1);
  MSerr = zeros(mdl.cv.niter, 1);

  % calculate the degrees of freedom of the model.
  DFtot = mdl.N - 1;
  DFfit = 2 * mdl.A;
  DFerr = DFtot - DFfit;

  % loop through the number of cross-validation iterations.
  for i = 1 : mdl.cv.niter
    % compute the total and fit sums of squares.
    SStot(i) = sumsq(vec(mdl.Y));
    SSfit(i) = sumsq(vec(mdl.CV{i}.Yh));
  end

  % calculate the error sum of squares.
  SSerr = SStot - SSfit;

  % calculate the mean square errors.
  MStot = SStot ./ DFtot;
  MSfit = SSfit ./ DFfit;
  MSerr = SSerr ./ DFerr;

  % store the sum of squares values in the output structure.
  S.SS.tot = median(SStot);
  S.SS.fit = median(SSfit);
  S.SS.err = median(SSerr);

  % store the mean square errors in the output structure.
  S.MS.tot = median(MStot);
  S.MS.fit = median(MSfit);
  S.MS.err = median(MSerr);

  % store the degrees of freedom in the output structure.
  S.DF.tot = DFtot;
  S.DF.fit = DFfit;
  S.DF.err = DFerr;

  % calculate the F-statistic from the mean square error ratio.
  S.F = S.MS.fit / S.MS.err;

  % calculate the p-value from the F cumulative distribution function.
  S.p = 1 - fcdf(S.F, S.DF.fit, S.DF.err);
end

