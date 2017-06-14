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
## @anchor{integralsplot}
## @deftypefn {Function File} {} integralsplot (@var{x}, @var{ab}, @var{Ix}, @var{Iab})
## Overlays integral curves generated from @ref{integrals} on a spectral
## line plot.
## @end deftypefn

function integralsplot (x, ab, Ix, Iab)
  % check the number of input arguments.
  if (nargin != 4 || !isvector(x) || !isvector(ab) || ...
      !isvector(Ix) || !isvector(Iab))
    % unknown. throw an exception.
    print_usage();
  end

  % scale the integrals to the spectral range.
  Ix = Ix .* (range(x) / range(Ix));

  % find the jump locations in the integral abscissa.
  kstops = findjumps(Iab);
  kstops = [1; kstops; length(Iab) + 1];

  % start the plot.
  figure();
  hold on;

  % plot the spectrum.
  plot(ab, x);

  % loop through the jumps.
  for k = 1 : length(kstops) - 1
    % get the bounds.
    k1 = kstops(k);
    k2 = kstops(k + 1) - 1;

    % slice the integral data.
    Ii = Ix(k1 : k2);
    Iabi = Iab(k1 : k2);

    % add the curve.
    plot(Iabi, Ii, 'color', [0, 0, 0]);
  end

  % finish the figure.
  hold off;
end

