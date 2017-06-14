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
## @anchor{dmxcorr}
## @deftypefn {Function File} {@var{Fcorr} =} dmxcorr (@var{F}, @var{parms})
## Corrects Bruker DMX-format for group delay issues. This function is usually
## called automatically (by @ref{loaddmx}) without the user having to, but
## may have to be used if zero filling or apodizing must be applied before
## correction.
## @end deftypefn

function Fcorr = dmxcorr (F, parms)
  % check for proper arguments.
  if (nargin != 2 || !(ismatrix(F) || isvector(F) || iscell(F)) || ...
      !(isstruct(parms) || iscell(parms)))
    % improper arguments. print the usage statement.
    print_usage();
  end

  % declare a lookup table. yes, it's a barbaric hack. don't judge me.
  kv = [[   2, 44.7500, 46.0000, 46.311]; ...
        [   3, 33.5000, 36.5000, 36.530]; ...
        [   4, 66.6250, 48.0000, 47.870]; ...
        [   6, 59.0833, 50.1667, 50.229]; ...
        [   8, 68.5625, 53.2500, 53.289]; ...
        [  12, 60.3750, 69.5000, 69.551]; ...
        [  16, 69.5313, 72.2500, 71.600]; ...
        [  24, 61.0208, 70.1667, 70.184]; ...
        [  32, 70.0156, 72.7500, 72.138]; ...
        [  48, 61.3438, 70.5000, 70.528]; ...
        [  64, 70.2578, 73.0000, 72.348]; ...
        [  96, 61.5052, 70.6667, 70.700]; ...
        [ 128, 70.3789, 72.5000, 72.524]; ...
        [ 192, 61.5859, 71.3333,  0.000]; ...
        [ 256, 70.4395, 72.2500,  0.000]; ...
        [ 384, 61.6263, 71.6667,  0.000]; ...
        [ 512, 70.4697, 72.1250,  0.000]; ...
        [ 768, 61.6465, 71.8333,  0.000]; ...
        [1024, 70.4849, 72.0625,  0.000]; ...
        [1536, 61.6566, 71.9167,  0.000]; ...
        [2048, 70.4924, 72.0313,  0.000]];

  % determine the data type of the parms.
  if (iscell(parms))
    % use only the first dimension of the parms.
    parms = parms{1};
  end

  % determine the dsp firmware version.
  kvi = find(kv(:,1) == parms.decim);
  kvj = parms.dspfvs - 8;
  k = floor(kv(kvi,kvj));

  % determine the data type of the fid.
  if (isvector(F))
    % shift the vector.
    Fcorr = shift(F, -k);
  elseif (ismatrix(F))
    % shift the rows of the matrix.
    Fcorr = shift(F, -k, 2);
  elseif (iscell(F))
    % initialize the output array.
    Fcorr = F;

    % loop through the entries of the array.
    for idx = 1 : length(F)
      % shift the rows of the current array entry.
      Fcorr{idx} = shift(F{idx}, -k, 2);
    end
  else
    % invalid type. throw an exception.
    error('dmxcorr: invalid fid data type. must be a vector, matrix or cell');
  end
end

