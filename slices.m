## Copyright (C) 2014 University of Nebraska Board of Regents.
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
## @anchor{slices}
## @deftypefn {Function File} {[@var{Xcos}, @var{Xsin}] =} slices (@var{X})
## @deftypefnx {Function File} {@var{Xcos} =} slices (@var{X})
## De-interlaces States/Haberkorn/Ruben cosine- and sine-modulated rows
## of a complex matrix @var{X} into the complex matrices @var{Xcos} and
## @var{Xsin}. See @ref{states} for more information.
## @end deftypefn

function [Xcos, Xsin] = slices (X)
  % check for proper arguments.
  if (nargin != 1 || !any(nargout == [1 : 2]))
    % improper arguments. print the usage statement.
    print_usage();
  end

  % extract the sine and cosine submatrices.
  Xcos = X(2 : 2 : end, :);
  Xsin = X(1 : 2 : end - 1, :);
end

