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
## @anchor{th_soft}
## @deftypefn {Function File} {[@var{y}, @var{xth}] =} th_soft (@var{x}, @var{tau})
## Compute the soft thresholded vector @var{y} of an input vector @var{x}.
##
## If @var{tau} is not provided, a default value of @code{0.02} will be
## assumed.
## @end deftypefn

function [y, xth] = th_soft (x, tau)
  % check for proper arguments.
  if (!any(nargin == [1 : 2]) || !any(nargout == [1 : 2]))
    % improper arguments. print the usage statement.
    print_usage();
  end

  % see if a threshold value was provided.
  if (nargin < 2 || isempty(tau) || !isscalar(tau))
    % no. set the default value.
    tau = 0.02;
  end

  % compute the threshold value.
  th = 1 - tau;

  % find the indices that are above the threshold.
  keepers = (abs(x) > th * max(abs(x)));

  % compute the soft thresholded vector.
  y = keepers .* x;
  xth = x - tau .* y;
end

