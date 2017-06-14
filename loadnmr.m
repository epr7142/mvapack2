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
## @anchor{loadnmr}
## @deftypefn {Function File} {@var{f} =} loadnmr (@var{dirname})
## @deftypefnx {Function File} {[@var{f}, @var{parms}] =} loadnmr (@var{dirname})
## @deftypefnx {Function File} {[@var{f}, @var{parms}, @var{t}] =} loadnmr (@var{dirname})
## @deftypefnx {Function File} {[@var{f}, @var{parms}, @var{t}] =} loadnmr (@var{dirname}, @var{doswap})
## @deftypefnx {Function File} {@var{F} =} loadbruker (@var{dirnames})
## @deftypefnx {Function File} {[@var{F}, @var{parms}] =} loadnmr (@var{dirnames})
## @deftypefnx {Function File} {[@var{F}, @var{parms}, @var{t}] =} loadnmr (@var{dirnames})
## @deftypefnx {Function File} {[@var{F}, @var{parms}, @var{t}] =} loadnmr (@var{dirnames}, @var{doswap})
## Loads one or more Bruker or Agilent fid/ser files, automatically extracting
## parameters and automatically determining whether to parse Bruker or Agilent
## data.
##
## Some older data needs to be byte-swapped when it is loaded in. To do this,
## pass @var{doswap} as @code{true} to this function. The default behavior is
## to skip the byte swap.
##
## The parameters can be optionally returned if a second return value is
## requested. If a third optional return value is requested, the
## time abscissa (@var{t}) will be returned.
##
## @strong{NOTE:} Extreme care must be taken to ensure that the acquisition
## parameters of all experimental data specified in @var{dirnames} are
## totally identical! Only one parameter structure (@var{parms}) will be
## returned, so it is assumed that the first experiment holds parameters
## that are representative of the entire dataset.
## @end deftypefn

function [F, parms, t] = loadnmr (dirnames, doswap)
  % check for the minimum argument count.
  if (!any(nargin == [1 : 2]) || !any(nargout == [1 : 3]))
    % this won't work. throw an error.
    print_usage();
  end

  % check the type of the directory input argument.
  if (ischar(dirnames))
    % turn the single string into a cell array...
    dirnames = {dirnames};
  elseif (iscellstr(dirnames))
    % make sure at least one directory was passed.
    if (length(dirnames) == 0)
      % no. throw an exception.
      error('loadnmr: zero input directories supplied');
    end
  else
    % invalid input type. throw an exception.
    error('loadnmr: dirnames must be a string or a string cell array');
  end

  % get the acquisition parameter structure. we have to assume that the
  % parameters are identical across all spectra. if not, incorrect data
  % will ensue.
  p = acquparms(dirnames{1});

  % initialize the output fid data matrix and spectrum count.
  n = length(dirnames);
  F = [];

  % check if the user supplied a byte-swap argument.
  if (nargin < 2)
    % no. default to false.
    doswap = false;
  end

  % check if we are loading 'fid' or 'ser' files.
  if (length(p) > 1 && fexist([dirnames{1}, '/ser']))
    % two-dimensional.
    nd = 2;

    % allocate the cell array.
    F = cell(n, 1);

    % loop through the directories.
    for i = 1 : n
      % build the ser filename string.
      sername = [dirnames{i}, '/ser'];

      % use loadser to pull in the ser data.
      ser = loadser(sername, p, doswap);

      % store the parsed ser to the output cell array.
      F{i} = ser;
    end

    % get the number of data points.
    npts = size(F{1}');

    % see if the output ser should simply be a matrix.
    if (n == 1)
      % yes. squash the cell array down to a matrix.
      F = F{1};
    end
  elseif (fexist([dirnames{1}, '/fid']))
    % one-dimensional.
    nd = 1;

    % ensure only the first dimension of the parms array survives.
    p(2:end) = [];

    % loop through the directories.
    for i = 1 : n
      % build the fid filename string.
      fidname = [dirnames{i}, '/fid'];

      % use loadfid to pull in the fid data.
      fid = loadfid(fidname, p{1}, doswap);

      % append the parsed fid to the output matrix as a row.
      F = [F; fid'];
    end

    % get the number of data points.
    npts = columns(F);

    % see if the output fid matrix is to be reshaped.
    if (n == 1)
      % yes. reshape the fid into a column vector.
      F = reshape(F, npts, 1);
    end
  else
    % throw an exception.
    error('loadnmr: input data must be either one- or two-dimensional');
  end

  % see if the user requested parameters to be returned.
  if (nargout >= 2)
    % yes. return the parameter structure.
    parms = p;
  end

  % see if the user requested the time axis to be returned.
  if (nargout >= 3)
    % yes. build and return time values.
    if (nd == 1)
      % build a column vector.
      t = gentime(p{1}.td / 2, p{1});
    else
      % build a cell array.
      t = cell(nd, 1);

      % loop through each dimension.
      for idx = 1 : nd
        % check if the dimension is nonuniformly sampled.
        if (isfield(p{idx}, 'nus') && ...
            isfield(p{idx}, 'tdnus') && ...
            isfield(p{idx}, 'sched') && ...
            p{idx}.nus == true)
          % rebuild and subsample the time vector.
          tunif = gentime(p{idx}.tdnus / 2, p{idx});
          t{idx} = tunif(p{idx}.sched + 1);
        else
          % build a uniform column vector.
          t{idx} = gentime(p{idx}.td / 2, p{idx});
        end
      end
    end
  end
end

