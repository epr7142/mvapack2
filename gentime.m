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
## @anchor{gentime}
## @deftypefn {Function File} {@var{t} =} gentime (@var{n}, @var{parms})
## Uses spectral parameters for number of real data points (@var{n}), spectral
## width (@var{sw}) and carrier offset in ppm (@var{car}) to build a vector of
## chemical shifts for NMR spectra.
##
## This function returns only a single dimension axis at a time. Thus,
## @var{parms} must be a structure.
## @end deftypefn

function t = gentime (n, parms)
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
      % nope... throw an exception.
      error('gentime: data point count must be a scalar value');
    end
  end

  % check the parameter structure.
  if (!isstruct(parms) || !isfield(parms, 'aq'))
    % invalid structure. throw an exception.
    error('gentime: invalid parameter structure');
  end

  % generate the output time vector.
  t = linspace(0, parms.aq, n)';
end

