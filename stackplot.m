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
## @anchor{stackplot}
## @deftypefn {Function File} {} stackplot (@var{X})
## @deftypefnx {Function File} {} stackplot (@var{X}, @var{ab})
## Builds a stacked line plot of time-domain FID data or frequency-
## domain spectral data. If time-domain data is plotted, @var{ab}
## must hold the time values for the abcissa. If frequency-domain
## values are to be plotted, @var{ab} must have chemical shifts.
## The @var{ab} vector is optional.
## @end deftypefn

function h = stackplot (a, b)
  % check the number of input arguments.
  if (nargin == 1 && ismatrix(a))
    % build a false abscissa.
    X = a;
    ab = [1 : columns(X)]';
  elseif (nargin == 2 && ismatrix(a) && isvector(b))
    % use the given abscissa.
    X = a;
    ab = b;
  else
    % unknown. throw an exception.
    print_usage();
  end

  % get the size of the data matrix.
  [N, K] = size(X);

  % calculate the separation between spectra.
  sep = mean(range(real(X), 2));

  % plot the rows of X.
  figure();
  hold on;
  htmp = plot(ab, linspace(0, N * sep, N)' * ones(1, K) + real(X));
  hold off;

  % see if the plot handle should be returned.
  if (nargout == 1)
    % yes. return the handle.
    h = htmp;
  end
end

