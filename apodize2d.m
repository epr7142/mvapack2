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
## @anchor{apodize2d}
## @deftypefn {Function File} {@var{wfid} =} apodize2d (@var{fid}, @var{parms}, @var{fn}, @var{opts})
## Performs apodization of a two-dimensional time-domain NMR free-induction
## decay in order to alleviate truncation artifacts that can arise from Fourier
## transformation.
##
## Instead of using this function directly, it is recommended that you use
## @ref{apodize}.
## @end deftypefn

function wfid = apodize2d (fid, parms, fn, opts)
  % check for the minimum number of arguments.
  if (!any(nargin == [2 : 4]) || nargout != 1)
    % print the usage statement.
    print_usage();
  end

  % determine the apodization options of each dimension.
  if (!iscell(opts) || length(opts) != 2)
    % duplicate the options into a 2x1 cell array, as expected later.
    oldopts = opts;
    opts = cell(2, 1);
    opts{1} = oldopts;
    opts{2} = oldopts;
  end

  % get the time axes and apodization windows.
  t = cell(2, 1);
  w = cell(2, 1);
  for idx = 1 : 2
    % compute the axis and the window coefficients.
    t{idx} = gentime(parms{idx}.td / 2, parms{idx});
    w{idx} = fn(t{idx}, opts{idx});
  end

  % check the type of the input data.
  if (ismatrix(fid))
    % apply the direct-dimension apodization.
    wfid = (ones(rows(fid), 1) * w{1}') .* fid;

    % de-interlace, apodize, and re-interlace the input matrix.
    [A, B] = states(wfid);
    A = (w{2} * ones(1, columns(A))) .* A;
    B = (w{2} * ones(1, columns(B))) .* B;
    wfid = states(A, B);
  elseif (iscell(fid))
    % initialize the output cell array.
    wfid = cell(size(fid));
    
    % apodize each matrix in the cell array.
    for idx = 1 : length(fid)
      % apodize the matrix.
      wfid{idx} = apodize2d(fid{idx}, parms, fn, opts);
    end
  else
    % throw an exception.
    error('apodize2d: fid data must be a matrix or a cell array');
  end
end

