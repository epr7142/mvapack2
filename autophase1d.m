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
## @anchor{autophase1d}
## @deftypefn {Function File} {[@var{sp}, @var{phc0}, @var{phc1}] =} autophase1d (@var{s}, @var{objective})
## Corrects the phase of a one-dimensional Fourier-transformed spectrum or
## spectral dataset found by simplex optimization. The method uses an entropy
## minimization objective during optimization (@xref{simplex_entropy}).
##
## Instead of using this function directly, it is recommended that you use
## @ref{autophase}.
## @end deftypefn

function [sp, phc0, phc1] = autophase1d (s, objective)
  % check for proper arguments.
  if (nargin != 2 || nargout != 3)
    % improper arguments. print the usage statement.
    print_usage();
  end

  % initialize the dimension count and the target dimension.
  n = 1;
  d = 1;

  % check the type of the input data.
  if (isvector(s))
    % initialize the phase correction matrix.
    phc = [0, 0];

    % run the phase correction.
    phc = __simplexphase(s, objective, false, false, phc, n, d);

    % return the phase correction values.
    phc0 = phc(1);
    phc1 = phc(2);

    % return the phased spectrum.
    sp = phase1d(s, phc0, phc1);
  elseif (ismatrix(s))
    % initialize the output phase correction vectors.
    phc0 = zeros(rows(s), 1);
    phc1 = zeros(rows(s), 1);

    % loop over all spectra in the matrix.
    for idx = 1 : rows(s)
      % initialize the phase correction matrix.
      phc = [0, 0];

      % run the phase correction.
      phc = __simplexphase(s(idx,:)', objective, false, false, phc, n, d);

      % store the phase correction values.
      phc0(idx) = phc(1);
      phc1(idx) = phc(2);
    end

    % return the phased spectral data matrix.
    sp = phase1d(s, phc0, phc1);
  else
    % throw an exception.
    error('autophase1d: input data must be a vector or a scalar');
  end
end

