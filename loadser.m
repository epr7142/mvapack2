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
## @anchor{loadser}
## @deftypefn {Function File} {@var{ser} =} loadser (@var{filename}, @var{parms}, @var{doswap})
## Loads a Bruker or Agilent ser file based on given parameters.
##
## This function requires that nmrPipe be installed on the system and
## the environment variables required to run nmrPipe are set up.
## @end deftypefn

function ser = loadser (filename, parms, doswap)
  % check for proper arguments.
  if (nargin != 3 || nargout != 1 || ...
      !ischar(filename) || !isbool(doswap) || !iscell(parms))
    % improper arguments. print the usage statement.
    print_usage();
  end

  % ensure the requested file actually exists.
  if (!fexist(filename))
    % file not found. throw an exception.
    error('loadser: file "%s" does not exist', filename);
  end

  % calculate the number of complex points in each dimension.
  n = cellfun(@(x) x.td / 2, parms);

  % check if the swap argument was chosen.
  if (doswap == true)
    doswap = 'no';
  else
    doswap = '';
  end

  % check the type of the source data and parameters.
  if (strcmp(parms{1}.type, 'bruker'))
    % check if the dsp firmware is unknown.
    if (parms{1}.dspfvs == 0)
      % use the default value.
      parms{1}.dspfvs = 12;
    end

    % build a bruk2pipe command string.
    cmd = ['bruk2pipe -in "', filename, '"', ...
           ' -bad 0.0 -', doswap, 'aswap -DMX ', ...
           '-decim ', num2str(parms{1}.decim), ' ', ...
           '-dspfvs ', num2str(parms{1}.dspfvs), ' ', ...
           '-grpdly ', num2str(parms{1}.grpdly), ' ', ...
           '-xN ', num2str(parms{1}.td), ' ', ...
           '-yN ', num2str(parms{2}.td), ' ', ...
           '-xT ', num2str(n(1)), ' ', ...
           '-yT ', num2str(n(2)), ' ', ...
           '-xMODE DQD ', ...
           '-yMODE Echo-AntiEcho ', ...
           '-xSW ', num2str(parms{1}.sw.hz), ' ', ...
           '-ySW ', num2str(parms{2}.sw.hz), ' ', ...
           '-xOBS ', num2str(parms{1}.obs), ' ', ...
           '-yOBS ', num2str(parms{2}.obs), ' ', ...
           '-xCAR ', num2str(parms{1}.car.ppm), ' ', ...
           '-yCAR ', num2str(parms{2}.car.ppm), ' ', ...
           '-xLAB ', parms{1}.nuc, ' ', ...
           '-yLAB ', parms{2}.nuc, ' ', ...
           '-ndim 2 -aq2D States -ov ', ...
           '| pipe2txt.tcl'];
  elseif (strcmp(parms{1}.type, 'agilent'))
    % build a var2pipe command string.
    cmd = ['var2pipe -in "', filename, '"', ...
           ' -', doswap, 'aswap ', ...
           '-xN ', num2str(parms{1}.td), ' ', ...
           '-yN ', num2str(parms{2}.td), ' ', ...
           '-xT ', num2str(n(1)), ' ', ...
           '-yT ', num2str(n(2)), ' ', ...
           '-xMODE Complex ', ...
           '-yMODE Echo-AntiEcho ', ...
           '-xSW ', num2str(parms{1}.sw.hz), ' ', ...
           '-ySW ', num2str(parms{2}.sw.hz), ' ', ...
           '-xOBS ', num2str(parms{1}.obs), ' ', ...
           '-yOBS ', num2str(parms{2}.obs), ' ', ...
           '-xCAR ', num2str(parms{1}.car.ppm), ' ', ...
           '-yCAR ', num2str(parms{2}.car.ppm), ' ', ...
           '-xLAB ', parms{1}.nuc, ' ', ...
           '-yLAB ', parms{2}.nuc, ' ', ...
           '-ndim 2 -aq2D States -ov ', ...
           '| pipe2txt.tcl'];
  else
    % output an error.
    error('loadser: unsupported ser data type');
  end

  % execute the command to load in the data.
  [status, output] = system(cmd);

  % check if the command returned successfully.
  if (status == 0)
    % get the output string.
    mtx = str2num(output);

    % build the complex matrix.
    ser = zeros(n(2), n(1));
    for idx = 1 : rows(mtx)
      % get the indices.
      a = round(mtx(idx, 2));
      b = round(mtx(idx, 1));

      % set the value at the indices.
      x = mtx(idx, 3) + i * mtx(idx, 4);
      ser(a,b) = x;
    end
  else
    % throw an exception.
    error('loadser: failed to extract ser from "%s"', filename);
  end
end

