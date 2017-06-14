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
## @anchor{nmrft2d}
## @deftypefn {Function File} {[@var{s}, @var{ppm}, @var{hz}] =} nmrft2d (@var{fid}, @var{parms}, @var{doshift})
## Performs Fourier transformation and shifting to produce a 2D NMR spectral
## matrix or a 2D NMR spectral cell array. Two-dimensional data must be
## arranged with slices of the direct-dimension along the rows.
##
## Instead of using this function directly, it is recommended that you use
## @ref{nmrft}.
## @end deftypefn

function [s, ppm, hz] = nmrft2d (fid, parms, doshift)
  % check if the number of expected arguments was passed.
  if (nargin != 3 || !any(nargout == [1 : 3]))
    % print the usage statement.
    print_usage();
  end

  % initialize the number of data points.
  npts = 0;

  % check the type of the input data.
  if (ismatrix(fid))
    % fourier transform the rows of the input matrix.
    % (this takes care of the direct dimension)
    s = fft(fid, [], 2);

    % shift the rows of the input matrix, if so requested.
    if (doshift == true)
      v = round(columns(s) / 2);
      s = shift(s, v, 2);
    end

    % de-interlace the half-transformed matrix.
    [A, B] = states(s);

    % fourier transform the columns of the states matrices.
    A = fft(A);
    B = fft(B);

    % shift the columns of the states matrices, if so requested.
    if (doshift == true)
      v = round(rows(A) / 2);
      A = shift(A, v);
      B = shift(B, v);
    end

    % re-interlace the fully transformed matrix.
    s = states(A, B);

    % use the matrix to obtain the number of points.
    npts = [columns(s), round(rows(s) / 2)];
  elseif (iscell(fid))
    % initialize the output cell array.
    s = cell(size(fid));
    
    % fourier transform each matrix in the cell array.
    for idx = 1 : length(fid)
      % fourier transform the matrix.
      s{idx} = nmrft2d(fid{idx}, parms, doshift);
    end

    % use the columns as the number of points.
    npts = [columns(s{1}), round(rows(s{1}) / 2)];
  else
    % throw an exception.
    error('nmrft2d: input data must be a matrix or cell array');
  end

  % build the abscissa vectors.
  ppm = cell(2, 1);
  hz = cell(2, 1);
  for idx = 1 : 2
    % build both chemical shift and hertz vectors.
    ppm{idx} = genppm(npts(idx), parms{idx});
    hz{idx} = 0.5 .* parms{idx}.sw.hz .* linspace(1, -1, npts(idx))';
  end

  % flip the direct dimension ppm abscissa, for some reason. (?)
  ppm{1} = flipud(ppm{1});
end

