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
## @anchor{phase}
## @deftypefn {Function File} {@var{sp} =} phase (@var{s}, @var{parms}, @var{phc0}, @var{phc1})
## Corrects the phase of a Fourier-transformed spectrum with a zero order
## correction @var{phc0} and a first order correction @var{phc1}. The first
## order correction is performed as a function of chemical shift offset from
## the spectral center frequency. The arguments @var{phc0} and @var{phc1}
## may be scalars or vectors, where the @var{i}-th elements of the vectors
## correct the @var{i}-th row of @var{s}.
## @end deftypefn

function sp = phase (s, parms, phc0, phc1)
  % check for the correct number of arguments.
  if (nargin != 4 || nargout != 1)
    % output an error.
    print_usage();
  end

  % get the dimensions in the dataset.
  nd = nmrdims(s, parms);

  % check whether we are phasing one- or two-dimensional data.
  if (nd == 1)
    % run a one-dimensional phase correction.
    sp = phase1d(s, phc0, phc1);
  elseif (nd == 2)
    % run a two-dimensional phase correction.
    sp = phase2d(s, phc0, phc1);
  else
    % throw an exception.
    error('phase: input data must be one- or two-dimensional');
  end
end

