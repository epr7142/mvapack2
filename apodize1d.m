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
## @anchor{apodize1d}
## @deftypefn {Function File} {@var{wfid} =} apodize1d (@var{fid}, @var{parms}, @var{fn}, @var{opts})
## Performs apodization of a one-dimensional time-domain NMR free-induction
## decay in order to alleviate truncation artifacts that can arise from Fourier
## transformation.
##
## Instead of using this function directly, it is recommended that you use
## @ref{apodize}.
## @end deftypefn

function wfid = apodize1d (fid, parms, fn, opts)
  % check for the minimum number of arguments.
  if (!any(nargin == [2 : 4]) || nargout != 1)
    % print the usage statement.
    print_usage();
  end

  % get the time axis.
  t = gentime(parms.td / 2, parms);

  % compute the apodization window.
  w = fn(t, opts);

  % check the type of the input data.
  if (isvector(fid))
    % apply the window to the vector input.
    wfid = fid .* w;
  elseif (ismatrix(fid))
    % apply the window to the matrix input.
    wfid = fid * diag(w);
  else
    % throw an exception.
    error('apodize1d: fid data must be a vector or a scalar');
  end
end

