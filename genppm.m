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
## @anchor{genppm}
## @deftypefn {Function File} {@var{ppm} =} genppm (@var{n}, @var{parms})
## Uses spectral parameters for number of real data points (@var{n}), spectral
## width and carrier offset in ppm to build a vector of chemical shifts for NMR
## spectra.
##
## This function returns only a single dimension axis at a time. Thus,
## @var{parms} must be a structure.
## @end deftypefn

function ppm = genppm (n, parms)
  % check for a minimum argument count.
  if (nargin != 2 || nargout != 1)
    % print the usage statement.
    print_usage();
  end

  % check the data type of n.
  if (!isscalar(n))
    % hmm... maybe the user accidentally passed a matrix or vector.
    if (isvector(n))
      % use the length of the vector as n.
      n = length(n);
    elseif (ismatrix(n))
      % use the columns of the matrix as n.
      n = columns(n);
    else
      % nope... throw an exception, but don't let on that we did all these
      % generous checks. we mustn't encourage stupidity...
      error('genppm: data point count must be a scalar value');
    end
  end

  % check the parameter structure.
  if (!isstruct(parms) || ...
      !isfield(parms, 'sw') || !isfield(parms.sw, 'ppm') || ...
      !isfield(parms, 'car') || !isfield(parms.car, 'ppm'))
    % invalid structure. throw an exception.
    error('genppm: invalid parameter structure');
  end

  % pull the spectral width and carrier offset from the parameter structure.
  sw = parms.sw.ppm;
  car = parms.car.ppm;

  % generate the output chemical shift vector.
  ppm = flipud(linspace(car - sw / 2, car + sw / 2, n)');
end

