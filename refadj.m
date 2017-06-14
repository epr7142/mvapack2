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
## @anchor{refadj}
## @deftypefn {Function File} {@var{ppmadj} =} refadj (@var{ppm}, @var{oldcs}, @var{newcs})
## Shifts the values in the chemical shift vector to move a reference position
## to zero chemical shift. This operates only on the abscissa vector, not on
## the data matrix.
##
## If you need to reference individual spectra within a data matrix to a
## common chemical shift axis, use @ref{coshift}.
## @end deftypefn

function ppmadj = refadj (ppm, oldcs, newcs)
  % check for proper arguments.
  if (nargin != 3 || nargout != 1 || !isvector(ppm) || ...
      !isscalar(oldcs) || !isscalar(newcs))
    % improper arguments. print the usage statement.
    print_usage();
  end

  % subtract the reference value from the chemical shift vector.
  ppmadj = ppm + (newcs - oldcs);
end

