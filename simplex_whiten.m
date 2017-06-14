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
## @anchor{simplex_whiten}
## @deftypefn {Function File} {@var{obj} =} simplex_whiten (@var{s}, @var{phc})
## Objective function for @ref{autophase} that minimizes the number of
## `colored pixels' of the real spectral component:
##
## @quotation
## G. Balacco, C. Cobas. `Automatic phase correction of 2D NMR spectra by a
## whitening method'. Mag. Res. Chem., 2009. 47: 322-327.
## @end quotation
## @end deftypefn

function obj = simplex_whiten (s, phc)
  % check the type of the spectral data.
  if (isvector(s))
    % 1D: calculate the phased spectrum.
    sp = phase1d(s, phc(1), phc(2));
    x = real(sp);
    A = abs(sp);
    t = mean(A);

    % compute the number of points outside the threshold.
    N = sum(x > t) + sum(x < -t);
    obj = -N;
  elseif (ismatrix(s))
    % 2D: calculate the phased spectrum.
    sp = states(phase2d(s, phc(:,1), phc(:,2)));
    x = real(sp);
    A = abs(sp);
    t = mean(vec(A));

    % compute the number of points outside the threshold.
    N = sum(vec(x) > t) + sum(vec(x) < -t);
    obj = -N;
  end
end

