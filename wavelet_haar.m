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
## @anchor{wavelet_haar}
## @deftypefn {Function File} {@var{y} =} wavelet_haar (@var{t})
## Calculates the Haar mother wavelet.
## @end deftypefn

function y = wavelet_haar (t)
  % check if the number of expected arguments was passed.
  if (nargin != 1 || (!isscalar(t) && !isvector(t)))
    % print the usage statement.
    print_usage();
  end

  % calculate the values of the mother wavelet.
  tnan = ifelse(t >= 0 & t < 1, t, nan);
  ia = find(ifelse(tnan < 0.5, 1, 0));
  ib = find(ifelse(tnan >= 0.5, 1, 0));
  y = zeros(size(t));
  y(ia) = 1;
  y(ib) = -1;
end

