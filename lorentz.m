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
## @anchor{lorentz}
## @deftypefn {Function File} {@var{s} =} lorentz (@var{ppm}, @var{par})
## Simulates a perfectly phased complex Lorentzian peak over an abscissa
## @var{ppm} according to the parameter vector @var{par}. If @var{par}
## is a K-by-3 matrix, then K peaks will be simulated and summed.
##
## The expected contents of each row of @var{par} are as follows:
##
## @code{par(:,1)}: Chemical shift, in @var{ppm} units. @*
## @code{par(:,2)}: Linewidth, in @var{ppm} units. @*
## @code{par(:,3)}: Amplitude, in absolute units.
## @end deftypefn

function s = lorentz (ppm, par)
  % check if the number of expected arguments was passed.
  if (nargin != 2 || nargout != 1 || !isvector(ppm))
    % print the usage statement.
    print_usage();
  end

  % see if the data type of the parameters is a vector or a matrix.
  if (isvector(par))
    % extract the parameters from the vector.
    omega = par(1);
    lambda = par(3);
    amplitude = par(2);

    % ensure the linewidth is greater than zero.
    if (lambda <= 0)
      lambda = 1e-3;
    end

    % generate the output lorentzian peak vector.
    s = (amplitude * lambda) ./ (lambda + i .* (ppm - omega));
  elseif (ismatrix(par))
    % initialize the output spectrum.
    s = zeros(size(ppm));

    % loop through the peaks in the parameter matrix.
    for k = 1 : rows(par)
      % add the current peak to the output spectrum.
      s += lorentz(ppm, par(k,:));
    end
  else
    % the parameter argument doesn't match our expectations.
    error('lorentz: parameters "par" must be a vector or matrix');
  end
end

