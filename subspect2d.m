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
## @anchor{subspect2d}
## @deftypefn {Function File} {[@var{ssub}, @var{absub}] =} subspect2d (@var{s}, @var{ab}, @var{parms}, @var{Fmin}, @var{Fmax})
## Extracts a band of frequencies (@math{[F_{min},F_{max}]}) from @var{s} into
## @var{ssub}.
##
## It is highly recommended that you use @ref{subspect} instead of calling this
## function directly.
## @end deftypefn

function [ssub, absub] = subspect2d (s, ab, parms, Fmin, Fmax)
  % check if the number of expected arguments was passed.
  if (nargin != 5 || nargout != 2)
    % print the usage statement.
    print_usage();
  end

  % check the time axis argument.
  if (!iscell(ab) || length(ab) != 2 || ...
      !isreal(ab{1}) || !isreal(ab{2}))
    % invalid type. throw an exception.
    error('subspect2d: abscissa must be an array of two real vectors');
  end

  % check the parameter argument.
  if (!iscell(parms) || length(parms) < 2 || ...
      !isstruct(parms{1}) || !isstruct(parms{2}))
    % invalid type. throw an exception.
    error('subspect2d: parameters must be an array of two structures');
  end

  % check the frequency boundary arguments.
  if (!isvector(Fmin) || !isreal(Fmin) || length(Fmin) != 2 ||...
      !isvector(Fmax) || !isreal(Fmax) || length(Fmax) != 2)
    % invalid type. throw an exception.
    error('subspect2d: frequency bounds must be real two-vector values');
  end

  % build an 'ROI' matrix from the frequencies.
  Fmin = reshape(Fmin, 2, 1);
  Fmax = reshape(Fmax, 2, 1);
  R = [Fmin, Fmax];

  % find the indices at which to crop.
  v = sort([findnearest(ab{1}, R(1,:)); findnearest(ab{2}, R(2,:))], 2);

  % return the abscissa cell array.
  absub = cell(2, 1);
  absub{1} = ab{1}(v(1,1) : v(1,2));
  absub{2} = ab{2}(v(2,1) : v(2,2));

  % apply the extraction operations.
  if (ismatrix(s))
    % return the cropped matrix.
    ssub = s(v(2,1) : v(2,2), v(1,1) : v(1,2));
  elseif (iscell(s))
    % initialize the output cell array.
    ssub = cell(length(s), 1);

    % loop through the matrices in the cell array.
    for idx = 1 : length(s)
      % store the cropped matrix.
      ssub{idx} = s{idx}(v(2,1) : v(2,2), v(1,1) : v(1,2));
    end
  else
    % throw an exception.
    error('subspect2d: input data must be a matrix or a cell array');
  end
end

