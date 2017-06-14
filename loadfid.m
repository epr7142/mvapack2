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
## @anchor{loadfid}
## @deftypefn {Function File} {@var{fid} =} loadfid (@var{filename}, @var{parms}, @var{doswap})
## Loads a Bruker or Agilent fid file based on given parameters.
##
## This function requires that nmrPipe be installed on the system and
## the environment variables required to run nmrPipe are set up.
## @end deftypefn

function fid = loadfid (filename, parms, doswap)
  % check for proper arguments.
  if (nargin != 3 || nargout != 1 || !ischar(filename) || ...
      !isbool(doswap) || !(isstruct(parms) || iscell(parms)))
    % improper arguments. print the usage statement.
    print_usage();
  end

  % ensure the requested file actually exists.
  if (!fexist(filename))
    % file not found. throw an exception.
    error('loadfid: file "%s" does not exist', filename);
  end

  % determine the data type of the parms.
  if (iscell(parms))
    % use only the first dimension of the parms.
    parms = parms{1};
  end

  % calculate the number of complex points.
  n = parms.td / 2;

  % check if the swap argument was chosen.
  if (doswap == true)
    doswap = 'no';
  else
    doswap = '';
  end

  % check the type of the source data and parameters.
  if (strcmp(parms.type, 'bruker'))
    % check if the dsp firmware is unknown.
    if (parms.dspfvs == 0)
      % use the default value.
      parms.dspfvs = 12;
    end

    % build a bruk2pipe command string.
    cmd = ['bruk2pipe -in "', filename, '"', ...
           ' -bad 0.0 -', doswap, 'aswap -DMX ', ...
           '-decim ', num2str(parms.decim), ' ', ...
           '-dspfvs ', num2str(parms.dspfvs), ' ', ...
           '-grpdly ', num2str(parms.grpdly), ' ', ...
           '-xN ', num2str(parms.td), ' ', ...
           '-xT ', num2str(n), ' ', ...
           '-xMODE DQD ', ...
           '-xSW ', num2str(parms.sw.hz), ' ', ...
           '-xOBS ', num2str(parms.obs), ' ', ...
           '-xCAR ', num2str(parms.car.ppm), ' ', ...
           '-xLAB 1H -ndim 1 -ov ', ...
           '| pipe2txt.tcl'];
  elseif (strcmp(parms.type, 'agilent'))
    % build a var2pipe command string.
    cmd = ['var2pipe -in "', filename, '"', ...
           ' -', doswap, 'aswap ', ...
           '-xN ', num2str(parms.td), ' ', ...
           '-xT ', num2str(n), ' ', ...
           '-xMODE Complex ', ...
           '-xSW ', num2str(parms.sw.hz), ' ', ...
           '-xOBS ', num2str(parms.obs), ' ', ...
           '-xCAR ', num2str(parms.car.ppm), ' ', ...
           '-xLAB 1H -ndim 1 -ov ', ...
           '| pipe2txt.tcl'];
  else
    % output an error.
    error('loadfid: unsupported fid data type');
  end

  % execute the command to load in the data.
  [status, output] = system(cmd);

  % check if the command returned successfully.
  if (status == 0)
    % get the output string and build the complex vector.
    mtx = str2num(output);
    fid = mtx(:,2) + i .* mtx(:,3);

    % check if the fid needs zeros appended to it.
    if (rows(fid) < n)
      % yeah, it's too short. append zeros.
      suffix = zeros(n - rows(fid), 1) + i .* zeros(n - rows(fid), 1);
      fid = [fid; suffix];
    end
  else
    % throw an exception.
    error('loadfid: failed to extract fid from "%s"', filename);
  end
end

