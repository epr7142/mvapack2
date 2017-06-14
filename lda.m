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
## @anchor{lda}
## @tex
## Linear Discriminant Analysis (LDA) generates a linear model of an input
## data matrix $X$ and response matrix $Y$ such that the model components
## capture a maximal amount of between-class variation. The model projects
## the (non-singular) input data into a discrimination space as follows:
## $$ T = X P $$
## In other words, the model is structured as follows:
## $$ X = T P' + E $$
## Where $P$ is a matrix of the first $A$ significant eigenvectors of the
## Fisher LDA matrix $S$:
## $$ P D P^{-1} = S = S_W^{-1} S_B $$
## And $S_W$ is the within-class covariance matrix and $S_B$ is the
## between-class covariance matrix. New observations are classified based on
## which class mean they fall closest to after projection by $P$.
##
## The number of components $A$ is chosen by
## cross-validation such that the model's cumulative $R^2_X$ and
## $R^2_Y$ metrics are greater than 0.01 and the model components'
## cumulative $Q^2$ metrics are strictly increasing.
## @end tex
## @deftypefn {Function File} {@var{mdl} =} lda (@var{X}, @var{Y})
## @deftypefnx {Function File} {@var{mdl} =} lda (@var{X}, @var{Y}, @var{scalefn})
## @deftypefnx {Function File} {@var{mdl} =} lda (@var{X}, @var{Y}, @var{scalefn}, @var{ncv})
## @deftypefnx {Function File} {@var{mdl} =} lda (@var{X}, @var{Y}, @var{scalefn}, @var{ncv}, @var{aout})
## Performs Linear Discriminant Analysis (LDA) on the input data matrix
## @var{X} and the input response matrix @var{Y} using an eigendecomposition of
## the within-class and between-class covariance matrix ratio. More information
## may be found here:
##
## @quotation
## W. Hardle, L. Simar. `Applied Multivariate Statistical Analysis, 2nd ed.',
## Springer-Verlag, Berlin Heidelberg, 2007.
## @end quotation
##
## The optional argument @var{scalefn} may be set to the function handle of
## an appropriate scaling routine. The default scaling function is @ref{suv}.
##
## The optional argument @var{ncv} may be used to specify the number of
## cross-validation iterations and/or groups. If @var{ncv} is a scalar, it
## will be used to set the number of CV groups, with 10 CV iterations. If
## @var{ncv} is a 2-element vector, its elements will be used to set the
## number of CV groups and iterations, respectively. The default is to use
## 7 groups and 10 iterations.
##
## The final optional argument @var{aout} may be used if a specific number of
## model components is desired. @strong{CAUTION:} using @var{aout} will not
## guarantee that the resultant model components are statistically
## significant!
## @end deftypefn

function mdl = lda (X, Y, scalefn, ncv, aout)
  % check for the minimum number of arguments.
  if (nargin < 1 || !ismatrix(X))
    % print the usage statement.
    print_usage();
  end

  % set the type of the output structure.
  mdl.type = 'lda';
  mdl.builder = @lda;
  mdl.classifier = @ldaclassify;

  % check if a second argument was passed.
  if (nargin >= 3 && !isempty(scalefn) && is_function_handle(scalefn))
    % store the scaling function handle.
    mdl.scaling = scalefn;
  else
    % store a default scaling function handle.
    mdl.scaling = @suv;
  end

  % check if a cross-validation argument was passed.
  if (nargin >= 4 && !isempty(ncv))
    % vector or scalar?
    if (isvector(ncv) && length(ncv) == 2)
      % use the argument as the groups and iterations.
      mdl.cv.ngroup = ncv(1);
      mdl.cv.niter = ncv(2);
    elseif (isscalar(ncv))
      % use the argument as the number of groups.
      mdl.cv.ngroup = ncv;
      mdl.cv.niter = 10;
    else
      % output an error.
       error('lda: invalid construction of cross-validation argument');
    end
  else
    % no. use the default: 7 groups and 10 iterations.
    mdl.cv.ngroup = 7;
    mdl.cv.niter = 10;
  end

  % check for valid cross-validation parameters.
  if (mdl.cv.ngroup < 2 || mdl.cv.niter < 1)
    % throw an exception.
    error('lda: invalid cross-validation parameters');
  end

  % check if a target component count was passed.
  if (nargin < 5 || isempty(aout))
    % set the target count to zero, indicating a reliance on validation.
    aout = 0;
  end

  % check whether the model is of the discriminant type.
  if (!isclasses(Y))
    % throw an exception.
    error('lda: response matrix must contain class information');
  end

  % initialize the discrimination status.
  mdl.isda = true;

  % store the original data and response matrices.
  mdl.X0 = X;
  mdl.Y0 = Y;

  % store the scaled data and response matrices.
  [mdl.X, mdl.mean.X, mdl.scale.X] = mdl.scaling(X);
  [mdl.Y, mdl.mean.Y, mdl.scale.Y] = mdl.scaling(Y);

  % store the matrix dimensions.
  [mdl.N, mdl.K] = size(mdl.X);
  mdl.M = columns(mdl.Y);

  % ensure the covariance matrices will not be singular.
  if (mdl.N <= mdl.K)
    % whoops. this probably won't work well.
    error('lda: input data matrix must be of full rank');
  end

  % ensure the cv group count is less than the observation count.
  if (mdl.cv.ngroup > mdl.N)
    % whoops. throw an exception.
    error('lda: number of cross-validation groups exceeds observation count');
  end

  % ensure the observation count matches the response count.
  if (mdl.N != rows(mdl.Y))
    % whoops. throw an exception.
    error('lda: row counts of data and response matrices do not match');
  end

  % initialize the model fields.
  mdl.A = 0;
  mdl.U = [];
  mdl.T = [];
  mdl.P = [];
  mdl.E = mdl.X;
  mdl.F = mdl.Y0;

  % initialize model eigenvalue, R-squared and Q-squared values.
  mdl.lambda = [];
  mdl.Rsq.X.comp = [];
  mdl.Rsq.Y.comp = [];
  mdl.Qsq.comp = [];
  mdl.Qsq.dev = [];

  % initialize the cross-validation structure.
  mdl.CV = cvldainit(mdl);

  % calculate the maximum number of components.
  Amax = max([aout; min([mdl.N - 1; mdl.K - 1; mdl.M - 1])]);

  % calculate the eigendecomposition of the covariance matrix quotient.
  [P, D, mdl.U] = ldacomp(mdl.X, mdl.Y0);

  % calculate the cross-validated decomposition.
  [QsqCV, QstdCV, mdl.CV] = cvlda(mdl, Amax);

  % loop through all possibly significant components.
  for a = 1 : Amax
    % extract the new component of weights and coefficients.
    p = P(:,a);

    % compose a set of dimension-limited parameters.
    Pq = P(:,1:a);
    Dq = D(1:a);

    % build the estimated data matrix.
    Xhat = mdl.X * Pq;

    % build the estimated response matrix.
    Yhat = ldaclassify(Pq, mdl.U, mdl.X);

    % calculate the eigenvalue and R-squared values.
    lambda = D(a);
    RsqX = ssratio(Xhat, mdl.X * P);
    RsqY = 1 - ssratio(mdl.Y0 - Yhat, mdl.Y0);

    % see if we need to subtract the last values.
    if (a > 1)
      % yeah.
      RsqX -= sum(mdl.Rsq.X.comp);
      RsqY -= sum(mdl.Rsq.Y.comp);
    end

    % calculate the cross-validation Q-squared value.
    Qsq = min([RsqY; QsqCV(a)]);
    Qstd = QstdCV(a);

    % see if we have a target component count.
    if (aout)
      % see if we've reached it.
      if (mdl.A >= aout)
        % we have. break the loop.
        break;
      end
    else
      % check if the current component is significant.
      if (RsqX < 0.01 || RsqY < 0.00 || ...
          sum([mdl.Qsq.comp; Qsq]) < sum(mdl.Qsq.comp))
        % no. break from the loop.
        break;
      end
    end

    % store the new eigenvalue and R-squared values.
    mdl.lambda = [mdl.lambda; lambda];
    mdl.Rsq.X.comp = [mdl.Rsq.X.comp; RsqX];
    mdl.Rsq.Y.comp = [mdl.Rsq.Y.comp; RsqY];
    mdl.Qsq.comp = [mdl.Qsq.comp; Qsq];
    mdl.Qsq.dev = [mdl.Qsq.dev; Qstd];

    % store the new component.
    mdl.P = [mdl.P, p];

    % increment the component count.
    mdl.A++;
  end

  % check if any components are significant.
  if (mdl.A == 0)
    % no. throw an exception.
    error('lda: no significant components found');
  end

  % calculate the final scores for the model.
  mdl.T = mdl.X * mdl.P;

  % store the cumulative R-squared and Q-squared values.
  mdl.Rsq.X.cum = sum(mdl.Rsq.X.comp);
  mdl.Rsq.Y.cum = sum(mdl.Rsq.Y.comp);
  mdl.Qsq.cum = sum(mdl.Qsq.comp);
end

