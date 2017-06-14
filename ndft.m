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
## @anchor{ndft}
## @deftypefn {Function File} {@var{y} =} ndft (@var{x}, @var{t})
## @deftypefnx {Function File} {[@var{y}, @var{w}] =} ndft (@var{x}, @var{t})
## @deftypefnx {Function File} {@var{y} =} ndft (@var{x}, @var{t}, @var{w})
## Compute the non-uniform discrete Fourier transform of @var{x} using the
## brute-force (@math{O(N^2)}) NDFT algorithm.
##
## The NDFT is calculated along the first non-singleton dimension of the
## array. Thus if @var{x} is a matrix, @code{ndft (x, t)} computes the NDFT
## for each column of @var{x}.
##
## An optional second return value @var{w} may be requested that will hold the
## frequency domain axis, in whatever units correspond to those in @var{t}.
## Alternatively, a third argument (also @var{w}) may be passed to explicitly
## specify the frequencies at which to compute the NDFT.
## @end deftypefn

function [y, w] = ndft (x, t, w)
  % check the arguments.
  if (!any(nargin == [2 : 3]) || !any(nargout == [1 : 2]))
    % invalid arguments.
    print_usage();
  end

  % check the type of the time axis.
  if (!isvector(t) || !isreal(t))
    % invalid type. throw an exception.
    error('ndft: time axis is expected to be a real vector');
  end

  % reshape the time axis vector.
  t = reshape(t, length(t), 1);

  % check if a frequency axis was provided.
  if (nargin >= 3 && !isempty(w))
    % check the type of the argument.
    if (!isvector(w) || !isreal(w))
      % invalid argument. throw an exception.
      error('ndft: frequency axis is expected to be a real vector');
    end
  else
    % build a best guess of the frequency axis.
    w = linspace(0, 1 / min(diff(t)), length(t))';
  end

  % check the type of the input data.
  if (isvector(x))
    % check the lengths.
    if (length(x) != length(t))
      % invalid lengths. throw an exception.
      error('ndft: lengths of data and time vectors must match');
    end

    % reshape the data vector.
    x = reshape(x, length(x), 1);

    % calculate the NDFT.
    y = __ndft(x, t, w);
  elseif (ismatrix(x))
    % check the dimensions.
    if (rows(x) != length(t))
      % invalid lengths. throw an exception.
      error('ndft: data matrix row count and time axis length must match');
    end

    % initialize the output matrix.
    y = zeros(size(x));

    % loop through the data columns.
    for n = 1 : columns(x)
      % calculate the NDFT for the current column.
      y(:,n) = __ndft(x(:,n), t, w);
    end
  end
end

