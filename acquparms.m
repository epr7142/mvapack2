## Copyright (C) 2014 University of Nebraska Board of Regents.
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
## @anchor{acquparms}
## @deftypefn {Function File} {@var{p} =} acquparms (@var{dirname})
## Reads values from the key-value pairs found in Bruker `acqus' files and
## Varian/Agilent `procpar' files. The type of acquisition parameters to
## parse is determined automatically.
##
## The output @var{p} is a cell array that contains a structure for each
## dimension present in the spectral data. Each element (dimension) in @var{p}
## contains information relevant to reconstructing time- and frequency-domain
## axes (along that dimension) for the spectral intensities:
##
## Universal: @*
## @code{p.td}: number of complex points. @*
## @code{p.obs}: transmitter base frequency. @*
## @code{p.car}: carrier offset frequency. @*
## @code{p.sw}: spectral width.
## @code{p.dim}: dimension.
##
## Bruker only: @*
## @code{p.decim}: decimation ratio. @*
## @code{p.grpdly}: estimated group delay. @*
## @code{p.dspfvs}: dsp firmware version. @*
## @end deftypefn

function p = acquparms (dirname)
  % check for proper arguments.
  if (nargin != 1 || nargout != 1 || !ischar(dirname))
    % improper arguments. print the usage statement.
    print_usage();
  end

  % check if any readable parameter files exist.
  if (fexist([dirname, '/acqus']))
    % yes (bruker). does a second dimension parameter file exist?
    if (fexist([dirname, '/acqu2s']))
      % yes. return a cell array of two parameter structures.
      p = cell(2, 1);
      p{1} = acquparms_bruker(dirname, 'acqus', 1);
      p{2} = acquparms_bruker(dirname, 'acqu2s', 2);
    else
      % no. return a single parameter structure.
      p = {};
      p{1} = acquparms_bruker(dirname, 'acqus', 1);
    end
  elseif (fexist([dirname, '/procpar']))
    % yes (agilent). return a single parameter structure.
    p = {};
    p{1} = acquparms_agilent(dirname, 'procpar', 1);
  else
    % we could not find a parms file. output an error.
    error('acquparms: failed to find acqus or procpar in "%s"', dirname);
  end

  % loop through all dimensions in the parameter array.
  for idx = 1 : length(p)
    % see if the parameters were returned as a structure or a cell array.
    % convert the spectral parameters into ppm units.
    p{idx}.sw.ppm = p{idx}.sw.hz / p{idx}.obs;
    p{idx}.car.ppm = p{idx}.car.hz / p{idx}.obs;

    % calculate the acquisition time of the experiment.
    p{idx}.aq = p{idx}.td / (2 * p{idx}.sw.hz);
  end
end

