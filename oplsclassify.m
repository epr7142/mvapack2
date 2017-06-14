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
## @anchor{oplsclassify}
## @deftypefn {Function File} {@var{Y} =} oplsclassify (@var{mdl}, @var{X})
## @deftypefnx {Function File} {[@var{Y}, @var{T}] =} oplsclassify (@var{mdl}, @var{X})
## Predicts responses @var{Y} from one or more observations @var{X} based on
## the OPLS model provided in @var{mdl}. The observations in @var{X} are
## filtered based on orthogonal variation, transformed by the regression
## coefficients (@var{B}) and classified based on sum of squares to the
## model classes.
##
## @strong{NOTE:} this function is not meant to be used directly. If you
## want to use a OPLS model to classify new observations, use @ref{classify}.
## @end deftypefn

function [Y, T] = oplsclassify (mdl, X)
  % check the type of arguments.
  if (nargin != 2 || !isstruct(mdl) || !ismatrix(X))
    % invalid arguments. throw an exception.
    print_usage();
  end

  % get the number of observations to classify.
  n = rows(X);

  % initialize for transformation of the new observation matrix.
  E = X;
  To = [];
  Tp = [];
  ao = 1;

  % loop through the predictive model components.
  for ap = 1 : mdl.Ap
    % get the current predictive component.
    Wp = mdl.Wp(ap, :);
    Pp = mdl.Pp(ap, :);

    % get the orthogonal components for the current predictive component.
    Wo = mdl.Wo(:, ao : ao + mdl.no(ap) - 1);
    Po = mdl.Po(:, ao : ao + mdl.no(ap) - 1);

    % calculate the new orthogonal scores.
    to = (E * Wo) * diag(1 ./ diag(Wo' * Wo));
    To = [To, to];

    % calculate the new OSC-filtered observations.
    E -= to * Po';

    % calculate the new predictive scores.
    tp = (E * Wp) * diag(1 ./ diag(Wp' * Wp));
    Tp = [Tp, tp];

    % calculate the new residual.
    E -= tp * Pp';

    % increment to the next set of orthogonal components.
    ao += mdl.no(ap);
  end

  % subtract out orthogonal variation from the observations.
  Xp = X - To * mdl.Po';

  % transform the observations and means into the Y-space.
  Yhat = backscale(Xp * mdl.B, mdl.mean.Y, mdl.scale.Y);

  % initialize the output matrix.
  Y = [];
  YM = eye(mdl.M);

  % locate the maximum response for each observation.
  [dmax, Mmax] = max(Yhat, [], 2);

  % generate a response value for all observations simultaneously.
  Y = YM(Mmax,:);

  % see if a matrix of scores was requested.
  if (nargout >= 2)
    % yes. return both the predictive and orthogonal scores.
    T.p = Tp;
    T.o = To;
  end
end

