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
## @anchor{phase1d}
## @deftypefn {Function File} {@var{sp} =} phase1d (@var{s}, @var{phc0}, @var{phc1})
## Corrects the phase of a one-dimensional Fourier-transformed spectrum or
## spectral dataset with a zero order correction @var{phc0} and a first order
## correction @var{phc1}.
##
## Instead of using this function directly, it is recommended that you use
## @ref{phase}.
## @end deftypefn

function sp = phase1d (s, phc0, phc1)
  % check for the correct number of arguments.
  if (nargin != 3 || nargout != 1)
    % output an error.
    print_usage();
  end

  % check the type of the input data.
  if (isvector(s))
    % check the type of the phase correction arguments.
    if (!isscalar(phc0) || !isscalar(phc1))
      % throw an exception.
      error('phase1d: for vector "s", corrections must be scalars');
    end

    % build the abscissa from the input vector.
    ab = linspace(-1, 1, length(s))';

    % reshape the input vector into a column vector.
    s = reshape(s, length(s), 1);

    % build the phase correction vector and apply it to the input vector.
    phc = exp(-i .* (pi / 180) .* (phc0 + phc1 .* (ab - median(ab))));
    sp = phc .* s;
  elseif (ismatrix(s))
    % initialize the output matrix.
    sp = s;

    % build the abscissa from the input matrix.
    ab = linspace(-1, 1, columns(s))';

    % check if the phase values are scalars.
    if (isscalar(phc0) && isscalar(phc1))
      % up-convert them to vectors having the same size as the matrix.
      phc0 = ones(rows(s), 1) .* phc0;
      phc1 = ones(rows(s), 1) .* phc1;
    end

    % check the phase value types and lengths.
    if (!isvector(phc0) || length(phc0) != rows(s) || ...
        !isvector(phc1) || length(phc1) != rows(s))
      % throw an exception.
      error('phase1d: for matrix "s", corrections must be scalars or vectors');
    end

    % loop through the matrix rows.
    for n = 1 : rows(s)
      % build the phase correction vector.
      phc = exp(-i .* (pi / 180) .* (phc0(n) + phc1(n) .* (ab - median(ab))));

      % apply the correction value to the matrix row.
      sp(n,:) = phc' .* s(n,:);
    end
  else
    % invalid data matrix type. throw an exception.
    error('phase1d: "s" must be either a vector or a matrix');
  end
end

