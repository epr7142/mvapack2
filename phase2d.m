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
## @anchor{phase2d}
## @deftypefn {Function File} {@var{sp} =} phase2d (@var{s}, @var{phc0}, @var{phc1})
## Corrects the phase of a two-dimensional Fourier-transformed spectral
## matrix or cell array with a zero order correction @var{phc0} and a
## first order correction @var{phc1}.
##
## Instead of using this function directly, it is recommended that you use
## @ref{phase}.
## @end deftypefn

function sp = phase2d (s, phc0, phc1)
  % check for the correct number of arguments.
  if (nargin != 3 || nargout != 1)
    % output an error.
    print_usage();
  end

  % check the type of the input data.
  if (ismatrix(s))
    % check the type of the phase correction arguments.
    if (iscell(phc0) || !isvector(phc0) || length(phc0) != 2 || ...
        iscell(phc1) || !isvector(phc1) || length(phc1) != 2)
      % throw an exception.
      error('phase2d: for matrix "s", corrections must be two-vectors');
    end

    % build the abscissas from the input matrix.
    ab = cell(2, 1);
    ab{1} = linspace(-1, 1, columns(s))';
    ab{2} = linspace(-1, 1, rows(s) / 2)';

    % build the phase correction vectors for each dimension.
    phc = cell(2, 1);
    for idx = 1 : 2
      % compute the current dimension's phase correction vector.
      phc{idx} = phc0(idx) + phc1(idx) .* (ab{idx} - median(ab{idx}));
      phc{idx} = exp(-i .* (pi / 180) .* phc{idx});
    end

    % apply the f2 phase correction values to the matrix rows.
    % (this takes care of the direct dimension)
    sp = (ones(rows(s), 1) * phc{1}') .* s;

    % de-interlace, phase correct, and re-interlace the input matrix.
    % (this is for the indirect dimension)
    [A, B] = states(sp);
    A = (phc{2} * ones(1, columns(A))) .* A;
    B = (phc{2} * ones(1, columns(B))) .* B;
    sp = states(A, B);
  elseif (iscell(s))
    % initialize the output cell array.
    sp = cell(size(s));

    % check the type of the phase correction arguments.
    if (isvector(phc0) && length(phc0) == 2 && ...
        isvector(phc1) && length(phc1) == 2)
      % use identical phase correction values for each spectrum.
      phc0 = ones(length(s), 1) * reshape(phc0, 1, 2);
      phc1 = ones(length(s), 1) * reshape(phc1, 1, 2);
    end

    % ensure the phase correction arguments are suitable now.
    if (!(ismatrix(phc0) && columns(phc0) == 2 && rows(phc0) == length(s) && ...
          ismatrix(phc1) && columns(phc1) == 2 && rows(phc1) == length(s)))
      % throw an exception.
      error('phase2d: for cell "s", corrections must be vectors or matrices');
    end

    % phase correct each matrix in the cell array.
    for idx = 1 : length(s)
      % phase correct the matrix.
      sp{idx} = phase2d(s{idx}, phc0(idx,:), phc1(idx,:));
    end
  else
    % invalid data matrix type. throw an exception.
    error('phase2d: "s" must be either a matrix or a cell array');
  end
end

