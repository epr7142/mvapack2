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
## @anchor{nusft}
## @deftypefn {Function File} {@var{s} =} nusft (@var{fid}, @var{t})
## @deftypefnx {Function File} {[@var{s}, @var{ppm}] =} nusft (@var{fid}, @var{t}, @var{parms})
## @deftypefnx {Function File} {[@var{s}, @var{ppm}, @var{hz}] =} nusft (@var{fid}, @var{t}, @var{parms})
## Performs Fourier transformation and shifting to produce an NMR spectrum.
## The data in @var{fid} may either be a column vector or a data matrix
## where each free induction decay is arranged as a row in the matrix.
##
## This function differs from standard NMR Fourier transformation
## (@xref{nmrft}) solely because it uses a non-uniform discrete Fourier
## transform (NDFT) instead of the classical fast Fourier transform (FFT)
## to compute the spectrum, thus allowing the input samples to be arbitrarily
## spaced in time. However, this method of computing the spectrum suffers from
## serious drawbacks, and should not be used in any seriousness.
##
## If a parameter structure is passed as a second argument, a second output
## value will be produced which contains the chemical shift abscissa vector
## that is associated with @var{s}. Optionally, in this case, a third output
## value will be produced which contains the centered abscissa vector in
## hertz units, without the carrier offset applied (@var{hz}).
## @end deftypefn

function [s, ppm, hz] = nusft (fid, t, parms)
  % check if the number of expected arguments was passed.
  if (!any(nargin == [2 : 3]) || !any(nargout == [1 : 3]))
    % print the usage statement.
    print_usage();
  end

  % initialize the number of data points.
  npts = 0;

  % check the type of the input data.
  if (isvector(fid))
    % fourier transform and shift the input vector.
    v = floor(length(fid) / 2);
    [s, hztmp] = ndft(fid, t);
    s = circshift(s, v);

    % use the length as the number of points.
    npts = length(s);
  elseif (ismatrix(fid))
    % fourier transform and shift the rows of the input matrix.
    v = floor(columns(fid) / 2)
    [s, hztmp] = ndft(fid', t);
    s = shift(s', v, 2);

    % use the columns as the number of points.
    npts = columns(s);
  else
    % throw an exception.
    error('nusft: input data must be a vector or matrix');
  end

  % see if a second structure argument was supplied.
  if (nargin >= 2 && !isempty(parms))
    % ensure the second argument is a structure.
    if (!isstruct(parms))
      % throw an exception.
      error('nusft: parms must be a structure');
    end

    % ensure a second return value is requested. we don't want this
    % abscissa vector to fly off into the ether.
    if (nargout < 2)
      % throw an exception.
      error('nusft: chemical shift abscissa requested but not stored');
    end

    % shift the hertz vector over.
    hztmp -= range(hztmp) / 2;

    % build the chemical shift abscissa vector.
    ppm = flipud((hztmp + parms.car.hz) / parms.obs);

    % check if a third return value is requested.
    if (nargout >= 3)
      % yes. return the hertz abscissa vector.
      hz = hztmp;
    end
  end
end

