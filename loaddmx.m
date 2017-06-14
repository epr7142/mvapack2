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
## @anchor{loaddmx}
## @deftypefn {Function File} {@var{f} =} loaddmx (@var{dirname})
## @deftypefnx {Function File} {[@var{f}, @var{parms}] =} loaddmx (@var{dirname})
## @deftypefnx {Function File} {[@var{f}, @var{parms}, @var{t}] =} loaddmx (@var{dirname})
## @deftypefnx {Function File} {[@var{f}, @var{parms}, @var{t}] =} loaddmx (@var{dirname}, @var{correct})
## @deftypefnx {Function File} {@var{F} =} loaddmx (@var{dirnames})
## @deftypefnx {Function File} {[@var{F}, @var{parms}] =} loaddmx (@var{dirnames})
## @deftypefnx {Function File} {[@var{F}, @var{parms}, @var{t}] =} loaddmx (@var{dirnames})
## @deftypefnx {Function File} {[@var{F}, @var{parms}, @var{t}] =} loaddmx (@var{dirnames}, @var{correct})
## Loads one or more 1D Bruker DMX-format fid files, automatically extracting
## parameters. The parameters can be optionally returned if a second return
## value (@var{parms}) is requested. A third optional return value, @var{t},
## can be requested that will contain the time-domain abscissa.
##
## An optional second input argument @var{correct} may be supplied to enable
## or disable group delay correction. The default behavior is to correct for
## group delay, but if you plan on zero-filling or apodizing, you need to
## postpone the correction until after the filling operation. @var{F} may then
## be corrected explicitly for group delay (@xref{dmxcorr}).
## @end deftypefn

function [F, parms, t] = loaddmx (dirnames, correct)
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
      error('loaddmx: zero input directories supplied');
    end
  else
    % invalid input type. throw an exception.
    error('loaddmx: dirnames must be a string or a string cell array');
  end

  % see if the user supplied a correction argument.
  if (nargin < 2 || isempty(correct) || !isbool(correct))
    % no. default to true.
    correct = true;
  end

  % get the acquisition parameter structure.
  p = acquparms(dirnames{1});

  % ensure only the first dimension of the parms array survives. this will
  % force all functions to understand this as one-dimensional data.
  p(2:end) = [];

  % initialize the output fid data matrix and spectrum count.
  n = length(dirnames);
  F = [];

  % loop through the directories.
  for i = 1 : n
    % build the fid filename string.
    fidname = [dirnames{i}, '/fid'];

    % try to load the file.
    fh = fopen(fidname, 'rb');

    % does the file exist?
    if (fh == -1)
      % file not found. throw an exception.
      error('loaddmx: "%s" does not exist', fidname);
    end

    % close the file.
    fclose(fh);

    % access the external c++ code to pull in the binary fid data.
    fid = __loaddmx(fidname, p{1}.td);

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

  % see if the output matrix is to be corrected.
  if (correct == true)
    % yes. correct the output fid data.
    F = dmxcorr(F, p);
  end

  % see if the user requested parameters to be returned.
  if (nargout >= 2)
    % yes. return the parameter structure.
    parms = p;
  end

  % see if the user requested the time axis to be returned.
  if (nargout >= 3)
    % yes. build and return a time vector.
    t = gentime(p{1}.td / 2, p{1});
  end
end

