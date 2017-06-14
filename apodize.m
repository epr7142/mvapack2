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
## @anchor{apodize}
## @deftypefn {Function File} {@var{wfid} =} apodize (@var{fid}, @var{parms})
## @deftypefnx {Function File} {@var{wfid} =} apodize (@var{fid}, @var{parms}, @var{fn})
## @deftypefnx {Function File} {@var{wfid} =} apodize (@var{fid}, @var{parms}, @var{fn}, @var{opts})
## Performs apodization of a time-domain NMR free-induction decay in order to
## alleviate truncation artifacts that can arise from Fourier transformation.
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
## The third argument, @var{fn}, is used to specify the type of apodization
## to use (@ref{expwindow}, @ref{gausswindow}, @ref{sinewindow}). The
## default is an exponential window function. An optional fourth argument
## (@var{opts}) may be passed that holds the parameters needed by the
## apodization function @var{fn}.
## @end deftypefn

function wfid = apodize (fid, parms, fn, opts)
  % check if the number of expected arguments was passed.
  if (!any(nargin == [2 : 4]) || nargout != 1)
    % print the usage statement.
    print_usage();
  end

  % see if an apodization function was specified.
  if (nargin >= 3 && !isempty(fn))
    % check the type of the apodization function.
    if (!is_function_handle(fn))
      % invalid. throw an exception.
      error('apodize: apodization method must be a valid function handle');
    end
  else
    % use the default function.
    fn = @expwindow;
  end

  % see if apodization options were passed.
  if (nargin >= 4 && !isempty(opts))
    % check the type of the options.
    if (!isscalar(opts) && !isvector(opts) && !isstruct(opts))
      % invalid. throw an exception.
      error('apodize: the specified apodization options are invalid');
    end
  else
    % use a default option set.
    opts = [];
  end

  % get the dimensions in the dataset.
  nd = nmrdims(fid, parms);

  % check whether we are fourier-transforming one- or two-dimensional data.
  if (nd == 1)
    % see if we need to reduce the parms cell array to a struct.
    if (iscell(parms))
      parms = parms{1};
    end

    % run a one-dimensional apodization.
    wfid = apodize1d(fid, parms, fn, opts);
  elseif (nd == 2)
    % run a two-dimensional apodization.
    wfid = apodize2d(fid, parms, fn, opts);
  else
    % throw an exception.
    error('apodize: input data must be one- or two-dimensional');
  end
end

