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
## @anchor{simplex_entropy}
## @deftypefn {Function File} {@var{obj} =} simplex_entropy (@var{s}, @var{phc})
## Objective function for @ref{autophase} that minimizes the entropy of
## the first derivative of the real spectral component:
##
## @quotation
## L. Chen, Z. Weng, L. Y. Goh, M. Garland. `An efficient algorithm for
## automatic phase correction of NMR spectra based on entropy minimization'.
## J. Magn. Res., 2002. 158(2002): 164-168.
## @end quotation
## @end deftypefn

function obj = simplex_entropy (s, phc)
  % check the type of the spectral data.
  if (isvector(s))
    % 1D: calculate the phased spectrum.
    sp = phase1d(s, phc(1), phc(2));
    x = real(sp);

    % calculate the histogram and penalty terms.
    g = range(real(s));
    drv = abs(diff(x));
    hst = drv ./ sum(drv);
    pen = g * sum((x < 0) .* (x .^ 2));

    % return the objective function.
    obj = sum(hst .* log(hst)) - pen;
  elseif (ismatrix(s))
    % 2D: calculate the phased spectrum.
    sp = phase2d(s, phc(:,1), phc(:,2));
    x = real(states(sp));
    g = range(x, 2);

    % calculate the objective over each row of the real phased spectrum.
    objs = zeros(rows(x), 1);
    for idx = 1 : rows(x)
      % calculate the histogram and penalty terms.
      drv = abs(diff(x(idx,:)'));
      hst = drv ./ sum(drv);
      pen = g(idx) * sum((x(idx,:)' < 0) .* (x(idx,:)' .^ 2));

      % calculate the current row objective function.
      objs(idx) = sum(hst .* log(hst)) - pen;
    end

    % return the geometric mean of all row objective values.
    obj = mean(objs, 'g');
  end
end

