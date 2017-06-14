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
## @anchor{simplex_integral}
## @deftypefn {Function File} {@var{obj} =} simplex_integral (@var{s}, @var{phc})
## Objective function for @ref{autophase} that maximizes the integrated
## area of the real component of the spectrum.
## @end deftypefn

function obj = simplex_integral (s, phc)
  % check the type of the spectral data.
  if (isvector(s))
    % 1D: calculate the integral of the real phased spectrum.
    sp = phase1d(s, phc(1), phc(2));
    obj = trapz(real(sp));
  elseif (ismatrix(s))
    % 2D: calculate the integral of the real phased spectrum.
    sp = phase2d(s, phc(:,1), phc(:,2));
    obj = trapz(trapz(real(states(sp))));
  end
end

