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
## @anchor{opls}
## @tex
## Orthogonal Projections to Latent Structures (OPLS) generates a linear
## model of an input data matrix $X$ and response matrix $Y$ such that
## the model components are orthogonal and capture both a maximal amount
## of correlation and anti-correlation between $X$ and $Y$. The OPLS
## model is structured as follows:
## $$ X = T_p P_p' + T_o P_o' + E
##      = \sum_{a=1}^{A_p} t_{p,a} p_{p,a}'
##      + \sum_{b=1}^{A_o} t_{o,b} p_{o,b}' + E $$
## $$ Y = U C' + F = \sum_{a=1}^{A_p} u_a c_a' + F $$
## Where $T_p$ is an $N \times A_p$ matrix of predictive scores, $T_o$ is
## an $N \times A_o$ matrix of orthogonal scores, $P_p$ is a $K \times A_p$
## matrix of predictive $X$-loadings, $P_o$ is a $K \times A_o$ matrix of
## orthogonal $X$-loadings, $E$ is an $N \times K$ matrix of $X$ residuals,
## $U$ is an $N \times A_p$ matrix of $Y$-scores, $C$ is an $M \times A_p$
## matrix of $Y$-weights, and $F$ is an $N \times M$ matrix of $Y$
## residuals.
##
## The number of predictive components $A_p$ is chosen by
## cross-validation such that the model's cumulative $R^2_{X_p}$ and
## $R^2_Y$ metrics are greater than 0.01 and the model components'
## cumulative $Q^2$ metrics are strictly increasing. At least one
## orthogonal component is simultaneously extracted with each predictive
## component during fitting.
##
## OPLS scores in $T_p$ and $T_o$ are low-rank approxiations of predictive
## and orthogonal variation in observations in $X$, respectively. The OPLS
## loadings in $P_p$ and $P_o$ fill similar roles for the variables in
## $X$.
## @end tex
## @deftypefn {Function File} {@var{mdl} =} opls (@var{X}, @var{Y})
## @deftypefnx {Function File} {@var{mdl} =} opls (@var{X}, @var{Y}, @var{scalefn})
## @deftypefnx {Function File} {@var{mdl} =} opls (@var{X}, @var{Y}, @var{scalefn}, @var{ncv})
## @deftypefnx {Function File} {@var{mdl} =} opls (@var{X}, @var{Y}, @var{scalefn}, @var{ncv}, @var{aout})
## @deftypefnx {Function File} {@var{mdl} =} opls (@var{X}, @var{Y}, @var{scalefn}, @var{ncv}, @var{aout}, @var{w})
## Performs Orthogonal Projection to Latent Structures (OPLS) on the input
## data matrix @var{X} and response matrix @var{Y} using the iterative
## NIPALS method, described here:
##
## @quotation
## J. Trygg, S. Wold, `Orthogonal Projections to Latent Structures (O-PLS)',
## J. Chemometrics, 2002(16): 119-128.
## @end quotation
##
## The optional argument @var{scalefn} may be set to the function handle of
## an appropriate scaling routine. The default scaling function is
## @ref{spareto}.
##
## The optional argument @var{ncv} may be used to specify the number of
## cross-validation iterations and/or groups. If @var{ncv} is a scalar, it
## will be used to set the number of CV groups, with 10 CV iterations. If
## @var{ncv} is a 2-element vector, its elements will be used to set the
## number of CV groups and iterations, respectively. The default is to use
## 7 groups and 10 iterations.
##
## The final optional argument @var{aout} may be used if a specific number of
## model components is desired. If @var{aout} is specified as a scalar, then
## the model will be built with that number of predictive components. If
## @var{aout} is a vector, the first element should hold the desired number of
## predictive components. The remaining elements in @var{aout} should hold the
## desired number of orthogonal components for each predictive component.
## @strong{CAUTION:} using @var{aout} will not guarantee that the resultant
## model components are statistically significant!
## @end deftypefn

function mdl = opls (X, Y, scalefn, ncv, aout, w)
  % check for the minimum argument count.
  if (nargin < 2 || !ismatrix(X))
    % print the usage.
    print_usage();
  end

  % set the type of the output structure.
  mdl.type = 'opls';
  mdl.builder = @opls;
  mdl.classifier = @oplsclassify;

  % check if the third argument was passed.
  if (nargin >= 3 && !isempty(scalefn) && is_function_handle(scalefn))
    % store the passed scaling function handle.
    mdl.scaling = scalefn;
  else
    % store the default scaling function handle.
    mdl.scaling = @spareto;
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
       error('opls: invalid construction of cross-validation argument');
    end
  else
    % no. use the default: 7 groups and 10 iterations.
    mdl.cv.ngroup = 7;
    mdl.cv.niter = 10;
  end

  % check for valid cross-validation parameters.
  if (mdl.cv.ngroup < 2 || mdl.cv.niter < 1)
    % throw an exception.
    error('opls: invalid cross-validation parameters');
  end

  % check if a target component count was passed.
  if (nargin < 5 || isempty(aout))
    % set the target count to zero, indicating a reliance on validation.
    aout = 0;
  else
    % check the type of the component count.
    if (isscalar(aout))
      % rebuild the count to reflect no restriction on orthogonal components.
      aout = [aout; zeros(aout, 1)];
    elseif (isvector(aout))
      % ensure the vector is well formed.
      if (length(aout) != aout(1) + 1)
        % throw an exception.
        error('opls: improperly specified target component count vector');
      end
    else
      % throw an exception.
      error('opls: invalid target component count');
    end
  end

  % check if a weighting vector was passed.
  if (nargin < 6 || isempty(w))
    % set a unit weight.
    w = ones(1, columns(X));
  end

  % store whether the model is of the discriminant type.
  mdl.isda = isclasses(Y);

  % store the original data and response matrices.
  mdl.X0 = X;
  mdl.Y0 = Y;

  % store the data matrix variable weights.
  mdl.weight.X = w;
  mdl.weight.Y = ones(1, columns(Y));

  % store the scaled data and response matrices.
  [mdl.X, mdl.mean.X, mdl.scale.X] = mdl.scaling(X, w);
  [mdl.Y, mdl.mean.Y, mdl.scale.Y] = mdl.scaling(Y);

  % store the matrix dimensions.
  [mdl.N, mdl.K] = size(mdl.X);
  mdl.M = columns(mdl.Y);

  % initialize the model fields.
  mdl.no = [];
  mdl.Ap = 0;
  mdl.Ao = 0;
  mdl.A = 0;

  % ensure the cv group count is less than the observation count.
  if (mdl.cv.ngroup > mdl.N)
    % whoops. throw an exception.
    error('opls: number of cross-validation groups exceeds observation count');
  end

  % ensure the observation count matches the response count.
  if (mdl.N != rows(mdl.Y))
    % whoops. throw an exception.
    error('opls: row counts of data and response matrices do not match');
  end

  % initialize more model fields.
  mdl.Wp = [];
  mdl.Tp = [];
  mdl.Pp = [];
  mdl.Wo = [];
  mdl.To = [];
  mdl.Po = [];
  mdl.U = [];
  mdl.C = [];
  mdl.E = mdl.X;
  mdl.F = mdl.Y;

  % initialize the eigenvalue, R-squared and Q-squared values.
  mdl.iter = [];
  mdl.lambda = [];
  mdl.Rsq.Xp.comp = [];
  mdl.Rsq.Xo.comp = [];
  mdl.Rsq.Y.comp = [];
  mdl.Qsq.comp = [];
  mdl.Qsq.dev = [];

  % initialize the cross-validation structure.
  mdl.CV = cvoplsinit(mdl);

  % calculate the Y-orthogonal subspace of X.
  V = orthspace(mdl.X, mdl.Y);

  % calculate the maximum number of model components.
  Amax = max([aout(1); min([mdl.N / 2; mdl.K / 2; mdl.M])]);

  % loop until the maximum component count.
  for a = 1 : Amax
    % get the desired number of orthogonal components.
    if (a + 1 > length(aout))
      % no more specified.
      aorth = 0;
    else
      % still more specified.
      aorth = aout(a + 1);
    end

    % calculate the next component.
    [w, t, p, u, c, Wo, To, Po, iter] = oplscomp(mdl.E, mdl.F, V, aorth);

    % calculate the current component contribution.
    TPp = t * p';
    TPo = To * Po';
    TC = t * c';

    % calculate the eigenvalue, R-squared and Q-squared values.
    lambda = t' * t;
    RsqXp = ssratio(TPp, mdl.X);
    RsqXo = [];
    for ao = 1 : columns(Wo)
      RsqXo = [RsqXo; ssratio(To(:,ao) * Po(:,ao)', mdl.X)];
    end
    RsqY = ssratio(TC, mdl.Y);

    % calculate the new Q-squared value.
    [Qsq, Qstd, mdl.CV] = cvopls(mdl, V, columns(Wo));

    % see if we have a target component count.
    if (aout)
      % see if we've reached it.
      if (mdl.Ap >= aout(1))
        % we have. break the loop.
        break;
      end
    else
      % see if the current component is significant.
      if (RsqXp < 0.01 || RsqY < 0.01 ||
          sum([mdl.Qsq.comp; Qsq]) < sum(mdl.Qsq.comp))
        % not significant. break execution here.
        break;
      end
    end

    % store the new component information.
    mdl.Wp = [mdl.Wp, w];
    mdl.Tp = [mdl.Tp, t];
    mdl.Pp = [mdl.Pp, p];
    mdl.Wo = [mdl.Wo, Wo];
    mdl.To = [mdl.To, To];
    mdl.Po = [mdl.Po, Po];
    mdl.U = [mdl.U, u];
    mdl.C = [mdl.C, c];
    mdl.iter = [mdl.iter; iter];

    % store the eigenvalue, R-squared and Q-squared values.
    mdl.lambda = [mdl.lambda; lambda];
    mdl.Rsq.Xp.comp = [mdl.Rsq.Xp.comp; RsqXp];
    mdl.Rsq.Xo.comp = [mdl.Rsq.Xo.comp; RsqXo];
    mdl.Rsq.Y.comp = [mdl.Rsq.Y.comp; RsqY];
    mdl.Qsq.comp = [mdl.Qsq.comp; Qsq];
    mdl.Qsq.dev = [mdl.Qsq.dev; Qstd];

    % calculate the new residuals and component counts.
    mdl.E = mdl.E - TPp - TPo;
    mdl.F = mdl.F - TC;
    mdl.Ap++;

    % calculate the new orthogonal component count.
    mdl.no = [mdl.no; columns(Wo)];
    mdl.Ao += columns(Wo);
  end

  % store the total number of components.
  mdl.A = mdl.Ap + mdl.Ao;

  % calculate the original un-deflated weights and regression coefficients.
  mdl.Wps = mdl.Wp * inv(mdl.Pp' * mdl.Wp);
  mdl.Wos = mdl.Wo * inv(mdl.Po' * mdl.Wo);
  mdl.B = mdl.Wps * mdl.C';

  % concatenate the predictive and orthogonal components.
  mdl.W = [mdl.Wp, mdl.Wo];
  mdl.Ws = [mdl.Wps, mdl.Wos];
  mdl.T = [mdl.Tp, mdl.To];
  mdl.P = [mdl.Pp, mdl.Po];

  % calculate the leverages of the observations.
  mdl.Hx = diag(mdl.Tp * inv(mdl.Tp' * mdl.Tp) * mdl.Tp');
  mdl.Hy = diag(mdl.U * inv(mdl.U' * mdl.U) * mdl.U');

  % store the cumulative R-squared and Q-squared values.
  mdl.Rsq.Xp.cum = sum(mdl.Rsq.Xp.comp);
  mdl.Rsq.Xo.cum = sum(mdl.Rsq.Xo.comp);
  mdl.Rsq.Y.cum = sum(mdl.Rsq.Y.comp);
  mdl.Qsq.cum = sum(mdl.Qsq.comp);

  % check if any components were significant.
  if (mdl.Ap == 0)
    % no. throw an exception.
    error('opls: no significant predictive components found');
  end
end

