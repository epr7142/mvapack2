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
## @anchor{autophase2d}
## @deftypefn {Function File} {[@var{sp}, @var{phc0}, @var{phc1}] =} autophase2d (@var{s}, @var{objective})
## Corrects the phase of a two-dimensional Fourier-transformed spectrum or
## spectral dataset found by simplex optimization. The method uses a whitening
## objective during optimization (@xref{simplex_whiten}).
##
## Instead of using this function directly, it is recommended that you use
## @ref{autophase}.
## @end deftypefn

function [sp, phc0, phc1] = autophase2d (s, objective)
  % check for proper arguments.
  if (nargin != 2 || nargout != 3)
    % improper arguments. print the usage statement.
    print_usage();
  end

  % initialize the dimension count.
  n = 2;

  % check the type of the input data.
  if (ismatrix(s))
    % initialize the phase correction matrix.
    phc = [0, 0; 0, 0];

    % run the phase correction once over each dimension.
    phc = __simplexphase(s, objective, false, false, phc, n, 1);
    phc = __simplexphase(s, objective, false, false, phc, n, 2);

    % return the phase correction values.
    phc0 = phc(:,1);
    phc1 = phc(:,2);

    % return the phased spectrum.
    sp = phase2d(s, phc0, phc1);
  elseif (iscell(s))
    % initialize the output phase correction matrices.
    phc0 = zeros(length(s), 2);
    phc1 = zeros(length(s), 2);

    % loop over all spectra in the array.
    for idx = 1 : length(s)
      % initialize the phase correction matrix.
      phc = [0, 0; 0, 0];

      % run the phase correction once over each dimension.
      phc = __simplexphase(s{idx}, objective, false, false, phc, n, 1);
      phc = __simplexphase(s{idx}, objective, false, false, phc, n, 2);

      % store the phase correction values.
      phc0(idx,:) = reshape(phc(:,1), 2, 1);
      phc1(idx,:) = reshape(phc(:,2), 2, 1);
    end

    % return the phased spectral cell array.
    sp = phase2d(s, phc0, phc1);
  else
    % throw an exception.
    error('autophase2d: input data must be a matrix or a cell array');
  end
end

