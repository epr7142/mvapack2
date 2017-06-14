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
## @anchor{rjmcmc1d_reconstruct}
## @deftypefn {Function File} {@var{D} =} rjmcmc1d_reconstruct (@var{t}, @var{omega}, @var{rho})
## @deftypefnx {Function File} {@var{yhat} =} rjmcmc1d_reconstruct (@var{t}, @var{omega}, @var{rho}, @var{a})
## @deftypefnx {Function File} {@var{yhat} =} rjmcmc1d_reconstruct (@var{t}, @var{P})
## Reconstructs a one-dimensional time-domain free induction decay from
## parameters estimated by @ref{rjmcmc1d}. If only frequencies (@var{omega})
## and decay rates (@var{rho}) are provided, a basis of signals (@var{D}) is
## returned. If amplitudes are also provided in @var{a}, a final signal
## estimate will be returned.
## @end deftypefn

function x = rjmcmc1d_reconstruct (t, omega, rho, a)
  % check the arguments.
  if (!any(nargin == [2 : 4]) || nargout != 1)
    % invalid arguments.
    print_usage();
  end

  % check the type of the time axis.
  if (!isvector(t) || !isreal(t))
    % invalid type. throw an exception.
    error('rjmcmc1d_reconstruct: time axis must be a real vector');
  end

  % check if only two input arguments were provided.
  if (nargin == 2 && !isempty(omega) && ismatrix(omega))
    % 'rename' the input matrix.
    P = omega;

    % grab the parameters from the signal matrix.
    omega = P(:,2);
    rho = P(:,3);
    a = P(:,1);
  end

  % check the type of the frequency data.
  if (!isvector(omega) || !isreal(omega))
    % invalid type. throw an exception.
    error('rjmcmc1d_reconstruct: frequencies must be in a real vector');
  end

  % get the number of signals to model.
  K = length(omega);

  % check the type of the decay rate data.
  if (!isvector(rho) || !isreal(rho) || length(rho) != K)
    % invalid type. throw an exception.
    error('rjmcmc1d_reconstruct: decay rates must be in a real vector');
  end

  % reshape the time axis vector.
  t = reshape(t, length(t), 1);

  % reshape the frequencies and decay rates.
  omega = reshape(omega, K, 1);
  rho = reshape(rho, K, 1);

  % calculate the basis signals.
  D = __rjmcmc1d_basis(t, omega, rho);

  % check if we're returning a basis matrix or a signal vector.
  if (nargin == 2 || (nargin >= 4 && !isempty(a)))
    % check the type of the amplitude data.
    if (!isvector(a) || length(a) != K)
      % invalid type. throw an exception.
      error('rjmcmc1d_reconstruct: amplitudes must be in a vector');
    end

    % reshape the amplitudes.
    a = reshape(a, K, 1);

    % return the signal vector.
    x = D * a;
  else
    % return the basis matrix.
    x = D;
  end
end

