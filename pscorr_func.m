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
## @anchor{pscorr_func}
## @deftypefn {Function File} {@var{y} =} pscorr_func (@var{x}, @var{p})
## Objective function used by @ref{pscorr} to calculate optimal phasing and
## scatter correction of a set of NMR spectra. This is probably useless
## on its own.
## @end deftypefn

function y = pscorr_func (x, p)
  % build the phased and normalized objective input value.
  y = p(3) .* real(phase1d(complex(x(:,1), x(:,2)), p(1), p(2)));
end

