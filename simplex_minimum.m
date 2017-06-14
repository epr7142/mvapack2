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
## @anchor{simplex_minimum}
## @deftypefn {Function File} {@var{obj} =} simplex_minimum (@var{s}, @var{phc})
## Objective function for @ref{autophase} that maximizes the lowest real
## spectral point.
## @end deftypefn

function obj = simplex_minimum (s, phc)
  % check the type of the spectral data.
  if (isvector(s))
    % 1D: find the minimum real point of the phased spectrum.
    sp = phase1d(s, phc(1), phc(2));
    obj = min(real(sp));
  elseif (ismatrix(s))
    % 2D: find the minimum real point of the phased spectrum.
    sp = phase2d(s, phc(:,1), phc(:,2));
    obj = min(min(real(states(sp))));
  end
end

