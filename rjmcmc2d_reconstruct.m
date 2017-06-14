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
## @anchor{rjmcmc2d_reconstruct}
## @deftypefn {Function File} {@var{D} =} rjmcmc2d_reconstruct (@var{t}, @var{omega}, @var{rho})
## @deftypefnx {Function File} {@var{yhat} =} rjmcmc2d_reconstruct (@var{t}, @var{omega}, @var{rho}, @var{a})
## @deftypefnx {Function File} {@var{yhat} =} rjmcmc2d_reconstruct (@var{t}, @var{P})
## Reconstructs a two-dimensional time-domain free induction decay from
## parameters estimated by @ref{rjmcmc2d}. If only frequencies (@var{omega})
## and decay rates (@var{rho}) are provided, a basis of signals (@var{D}) is
## returned. If amplitudes are also provided in @var{a}, a final signal
## estimate will be returned.
## @end deftypefn

function x = rjmcmc2d_reconstruct (t, omega, rho, a)
  % check the arguments.
  if (!any(nargin == [2 : 4]) || nargout != 1)
    % invalid arguments.
    print_usage();
  end

  % check the type of the time array.
  if (!iscell(t) || length(t) != 2 || ...
      !isvector(t{1}) || !isreal(t{1}) || ...
      !isvector(t{2}) || !isreal(t{2}))
    % invalid type. throw an exception.
    error(['rjmcmc2d_reconstruct: ', ...
           'time axis must be an array of two real vectors']);
  end

  % check if only two input arguments were provided.
  if (nargin == 2 && !isempty(omega) && ismatrix(omega))
    % 'rename' the input matrix.
    P = omega;

    % grab the parameters from the signal matrix.
    omega = P(:,2:3);
    rho = P(:,4:5);
    a = P(:,1);
  end

  % check the type of the frequency data.
  if (!ismatrix(omega) || !isreal(omega) || columns(omega) != 2)
    % invalid type. throw an exception.
    error(['rjmcmc2d_reconstruct: ', ...
           'frequencies must be in a real two-column matrix']);
  end

  % get the number of signals to model.
  K = rows(omega);

  % check the type of the decay rate data.
  if (!ismatrix(rho) || !isreal(rho) || columns(rho) != 2 || rows(rho) != K)
    % invalid type. throw an exception.
    error(['rjmcmc2d_reconstruct: ', ...
           'decay rates must be in a real two-column matrix']);
  end

  % calculate the basis signals.
  D = __rjmcmc2d_basis(t, omega, rho);

  % check if we're returning a basis matrix or a signal vector.
  if (nargin == 2 || (nargin >= 4 && !isempty(a)))
    % check the type of the amplitude data.
    if (!isvector(a) || length(a) != K)
      % invalid type. throw an exception.
      error('rjmcmc2d_reconstruct: amplitudes must be in a vector');
    end

    % reshape the amplitudes.
    a = reshape(a, K, 1);

    % initialize the signal matrix.
    x = zeros(2 * rows(D{1, 1}), columns(D{1, 1}));

    % loop through the signals in the basis.
    for k = 1 : K
      % compute the cosine- and sine-modulated components.
      ycos = D{k, 1} .* a(k);
      ysin = D{k, 2} .* a(k);

      % join the components together.
      y = zeros(size(x));
      y(2 : 2 : end, :) = ycos;
      y(1 : 2 : end - 1, :) = ysin;

      % sum the current signal contribution into the output matrix.
      x += y;
    end
  else
    % initialize the basis cell array.
    x = cell(K, 1);

    % loop through each element of the array.
    for k = 1 : K
      % extract the relevant basis matrix.
      Dcosk = D{k, 1};
      Dsink = D{k, 2};

      % initialize the current basis matrix.
      D = zeros(2 * rows(Dcosk), columns(Dcosk));

      % interleave the cosine- and sine-modulated matrices together.
      D(2 : 2 : end, :) = Dcosk;
      D(1 : 2 : end - 1, :) = Dsink;

      % store the calculated basis matrix.
      x{k} = D;
    end
  end
end

