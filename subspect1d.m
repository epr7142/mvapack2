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
## @anchor{subspect1d}
## @deftypefn {Function File} {[@var{ssub}, @var{absub}] =} subspect1d (@var{s}, @var{ab}, @var{parms}, @var{Fmin}, @var{Fmax})
## Extracts a band of frequencies (@math{[F_{min},F_{max}]}) from @var{s} into
## @var{ssub}.
##
## It is highly recommended that you use @ref{subspect} instead of calling this
## function directly.
## @end deftypefn

function [ssub, absub] = subspect1d (s, ab, parms, Fmin, Fmax)
  % check if the number of expected arguments was passed.
  if (nargin != 5 || nargout != 2)
    % print the usage statement.
    print_usage();
  end

  % check the time axis argument.
  if (!isvector(ab) || !isreal(ab))
    % invalid type. throw an exception.
    error('subspect1d: abscissa must be a real vector');
  end

  % check the parameter argument.
  if (!isstruct(parms))
    % invalid type. throw an exception.
    error('subspect1d: parameters must be a structure');
  end

  % check the frequency boundary arguments.
  if (!isscalar(Fmin) || !isscalar(Fmax) || !isreal(Fmin) || !isreal(Fmax))
    % invalid type. throw an exception.
    error('subspect1d: frequency bounds must be real scalar values');
  end

  % build an 'ROI' matrix from the frequencies.
  R = [Fmin, Fmax];

  % find the indices at which to crop.
  v = sort(findnearest(ab, R));

  % apply the extraction operations.
  if (isvector(s))
    % return the cropped vectors.
    ssub = s(v(1) : v(2));
    absub = ab(v(1) : v(2));
  elseif (ismatrix(s))
    % return the cropped matrix and abscissa vector.
    ssub = s(:, v(1) : v(2));
    absub = ab(v(1) : v(2));
  else
    % throw an exception.
    error('subspect1d: input data must be a vector or a matrix');
  end
end

