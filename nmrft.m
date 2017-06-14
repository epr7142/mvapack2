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
## @anchor{nmrft}
## @deftypefn {Function File} {@var{s} =} nmrft (@var{fid}, @var{parms})
## @deftypefnx {Function File} {[@var{s}, @var{ppm}] =} nmrft (@var{fid}, @var{parms})
## @deftypefnx {Function File} {[@var{s}, @var{ppm}, @var{hz}] =} nmrft (@var{fid}, @var{parms})
## @deftypefnx {Function File} {@dots{} =} nmrft (@var{fid}, @var{parms}, @var{doshift})
## Performs Fourier transformation and shifting to produce an NMR spectrum.
##
## In the one-dimensional case, data in @var{fid} may either be a column
## vector or a data matrix where each free induction decay is arranged as
## a row in the matrix.
##
## In the two-dimensional case, data in @var{fid} may either be a data matrix
## where each direct-dimension slice is along the rows, or a cell array that
## contains multiple matrices, each having direct-dimension slices along
## its rows.
##
## A parameter structure (or array) must be passed as a second argument.
## Optionally, a second output value may be produced which contains the
## chemical shift abscissa vector(s) associated with @var{s} (@var{ppm}. Also,
## a third output value may be produced which contains the centered abscissa
## vector in hertz units, without the carrier offset applied (@var{hz}).
##
## A final optional argument, @var{doshift}, may be passed to enable or disable
## the default behavior of shifting the Fourier-transformed data by half the
## number of data points. The default behavior is to shift the data, and you
## should really never have a reason not to.
## @end deftypefn

function [s, ppm, hz] = nmrft (fid, parms, doshift)
  % check if the number of expected arguments was passed.
  if (!any(nargin == [2 : 3]) || !any(nargout == [1 : 3]))
    % print the usage statement.
    print_usage();
  end

  % check if a third argument was provided.
  if (nargin >= 3 && !isempty(doshift))
    % ensure the data type is correct.
    if (!isbool(doshift))
      % invalid type. throw an exception.
      error('nmrft: shifting option must be a boolean');
    end
  else
    % set the default value.
    doshift = true;
  end

  % get the dimensions in the dataset.
  nd = nmrdims(fid, parms);

  % check whether we are fourier-transforming one- or two-dimensional data.
  if (nd == 1)
    % see if we need to reduce the parms cell array to a struct.
    if (iscell(parms))
      parms = parms{1};
    end

    % run a one-dimensional fourier transform.
    [s, ppm, hz] = nmrft1d(fid, parms, doshift);
  elseif (nd == 2)
    % run a two-dimensional fourier transform.
    [s, ppm, hz] = nmrft2d(fid, parms, doshift);
  else
    % throw an exception.
    error('nmrft: input data must be one- or two-dimensional');
  end
end

