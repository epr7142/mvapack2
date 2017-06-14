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
## @anchor{histmatch_func}
## @deftypefn {Function File} {@var{err} =} histmatch_func (@var{alpha}, @var{Ht}, @var{Xs}, @var{Zmap})
## This is the objective function for @ref{histmatch}. You'll never call this
## function directly, if you're writing sane code.
## @end deftypefn

function err = histmatch_func (alpha, Ht, Xs, Zmap)
  % calculate a histogram of the log-transformed scaled spectra.
  Zs = log2(abs(alpha .* Xs) + 1);
  Hs = hist(Zs, Zmap);

  % return the difference between the two histograms as a scalar error.
  D = Ht - Hs;
  err = D * D';
end

