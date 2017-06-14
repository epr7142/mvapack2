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
## @anchor{gausswindow}
## @deftypefn {Function File} {@var{w} =} gausswindow (@var{t}, @var{lb})
## Calculate window coefficients for gaussian windowing of a signal. The
## parameter @var{lb} is the line-broadening factor (in Hz) to apply to
## the signal having time points in @var{t}.
##
## The default line-broadening factor is 0.3 Hz.
## @tex
##
## The equation applied to the free induction decay is as follows:
## $$ wfid(t) = fid(t) \cdot e^{-(lb \cdot t)^2} $$
## @end tex
##
## Instead of using this function directly, it is recommended that you use
## @ref{apodize}.
## @end deftypefn

function w = gausswindow (t, lb)
  % check for the minimum number of arguments.
  if (!any(nargin == [1 : 2]) || nargout != 1)
    % print the usage statement.
    print_usage();
  end

  % check for a passed line-broadening factor.
  if (nargin >= 2 && !isempty(lb))
    % yes. check the type of the factor.
    if (!isscalar(lb) || !isreal(lb))
      % invalid argument. throw an exception.
      error('gausswindow: line-broadening factor must be a real scalar');
    end
  else
    % no. use the default line-broadening factor.
    lb = 0.3;
  end

  % compute the window coefficients.
  w = exp(-(lb .* t) .^ 2);
end

