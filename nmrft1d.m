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
## @anchor{nmrft1d}
## @deftypefn {Function File} {[@var{s}, @var{ppm}, @var{hz}] =} nmrft1d (@var{fid}, @var{parms}, @var{doshift})
## Performs Fourier transformation and shifting to produce a 1D NMR spectrum
## or a 1D NMR spectral data matrix. One-dimensional data in @var{fid} may
## either be a column vector or a data matrix where each free induction
## decay is arranged as a row in the matrix.
##
## Instead of using this function directly, it is recommended that you use
## @ref{nmrft}.
## @end deftypefn

function [s, ppm, hz] = nmrft1d (fid, parms, doshift)
  % check if the number of expected arguments was passed.
  if (nargin != 3 || !any(nargout == [1 : 3]))
    % print the usage statement.
    print_usage();
  end

  % initialize the number of data points.
  npts = 0;

  % check the type of the input data.
  if (isvector(fid))
    % fourier transform the input vector.
    s = fft(fid);

    % shift the input vector, if so requested.
    if (doshift == true)
      v = round(rows(fid) / 2);
      s = circshift(s, v);
    end

    % use the length as the number of points.
    npts = length(s);
  elseif (ismatrix(fid))
    % fourier transform the rows of the input matrix.
    s = fft(fid, [], 2);

    % shift the rows of the input matrix, if so requested.
    if (doshift == true)
      v = round(columns(fid) / 2);
      s = shift(s, v, 2);
    end

    % use the columns as the number of points.
    npts = columns(s);
  else
    % throw an exception.
    error('nmrft1d: input data must be a vector or matrix');
  end

  % build the chemical shift abscissa vector.
  ppm = genppm(npts, parms);

  % return the hertz abscissa vector.
  hz = 0.5 .* parms.sw.hz .* linspace(-1, 1, length(ppm))';
end

