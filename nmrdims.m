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
## @anchor{nmrdims}
## @deftypefn {Function File} {@var{n} =} nmrdims (@var{x}, @var{parms})
## Determines whether the dataset given by @var{x} and @var{parms} is
## one- or two-dimensional, based on the expected data types that each
## dimensionality of data may assume. The data @var{x} may be real or
## complex, time or frequency domain, @emph{etc}.
## @end deftypefn

function n = nmrdims (x, parms)
  % check if the number of expected arguments was passed.
  if (nargin != 2 || nargout != 1)
    % print the usage statement.
    print_usage();
  end

  % check whether we are manipulating one- or two-dimensional data.
  if ((ismatrix(x) && (isstruct(parms) || length(parms) == 1)) || ...
      (isvector(x) && !iscell(x)))
    % one-dimensional.
    n = 1;
  elseif ((ismatrix(x) && iscell(parms) && length(parms) >= 2) || ...
          iscell(x))
    % two-dimensional.
    n = 2;
  else
    % umm... not sure.
    n = 0;
  end
end

