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
## @anchor{ldaclassify}
## @deftypefn {Function File} {@var{Y} =} ldaclassify (@var{mdl}, @var{X})
## @deftypefnx {Function File} {@var{Y} =} ldaclassify (@var{P}, @var{U}, @var{X})
## @deftypefnx {Function File} {[@var{Y}, @var{T}] =} ldaclassify (@var{mdl}, @var{X})
## @deftypefnx {Function File} {[@var{Y}, @var{T}] =} ldaclassify (@var{P}, @var{U}, @var{X})
## Predicts responses @var{Y} from one or more observations @var{X} based on
## the LDA either provided in @var{mdl} or as @var{P} and @var{U}. The
## observations in @var{X} are transformed into the discriminant space and
## classified based on Euclidean distances to the model classes.
##
## @strong{NOTE:} this function is not meant to be used directly. If you
## want to use an LDA model to classify new observations, use @ref{classify}.
## @end deftypefn

function [Y, T] = ldaclassify (a, b, c)
  % check the type of arguments.
  if (nargin == 2 && isstruct(a) && ismatrix(b))
    % get the model structure.
    mdl = a;

    % get the P and U matrices.
    P = mdl.P;
    U = mdl.U;

    % get the X matrix.
    X = b;
  elseif (nargin == 3 && ismatrix(a) && ismatrix(b) && ismatrix(c))
    % get the P and U matrices.
    P = a;
    U = b;

    % get the X matrix.
    X = c;
  else
    % throw an exception.
    print_usage();
  end

  % get the weight matrix size.
  [K, A] = size(P);
  M = rows(U);

  % get the number of observations to classify.
  n = rows(X);

  % transform the observations and means into the LD space.
  Xhat = X * P;
  Uhat = U * P;

  % initialize the output matrix.
  Y = [];
  YM = eye(M);

  % loop through the transformed observations.
  for i = 1 : n
    % calculate the distances of the observation to the means.
    d = sumsq(ones(M, 1) * Xhat(i,:) - Uhat, 2);

    % locate the minimum distance.
    [dmin, Mmin] = min(d);

    % generate a response value for the observation.
    Y = [Y; YM(Mmin,:)];
  end

  % see if a matrix of scores was requested.
  if (nargout >= 2)
    % yes. return the scores.
    T = Xhat;
  end
end

