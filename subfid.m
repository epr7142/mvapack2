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
## @anchor{subfid}
## @deftypefn {Function File} {[@var{fsub}, @var{tsub}, @var{dfbw}, @var{D}] =} subfid (@var{f}, @var{t}, @var{parms}, @var{Fmin}, @var{Fmax})
## Extracts a band of frequencies (@math{[F_{min},F_{max}]}) from @var{f} into
## @var{fsub}. This function can operate on multiple free induction decays at
## once if @var{f} is supplied as a data matrix. The function also accepts
## two-dimensional time-domain data in the form of matrices or cell arrays.
##
## For one-dimensional data, @var{Fmin} and @var{Fmax} must be scalar values.
## For two-dimensional data, they must be two-vectors containing frequencies
## for each dimension.
##
## Procedurally, the extraction is as follows: modulate @var{f} such that the
## desired frequency band lies at zero frequency, FIR filter @var{f} to avoid
## aliasing, and decimate the filtered form of @var{f}. The extracted band
## will be returned in (@var{fsub},@var{tsub}), and the decimation ratio will
## be returned in @var{D}.
## @end deftypefn

function [fsub, tsub, dfbw, D] = subfid (f, t, parms, Fmin, Fmax)
  % check if the number of expected arguments was passed.
  if (nargin != 5 || nargout != 4)
    % print the usage statement.
    print_usage();
  end

  % get the dimensions in the dataset.
  nd = nmrdims(f, parms);

  % check whether we are filtering one- or two-dimensional data.
  if (nd == 1)
    % see if we need to reduce the parms cell array to a struct.
    if (iscell(parms))
      parms = parms{1};
    end

    % apply a one-dimensional filter.
    [fsub, tsub, dfbw, D] = subfid1d(f, t, parms, Fmin, Fmax);
  elseif (nd == 2)
    % apply a two-dimensional filter.
    [fsub, tsub, dfbw, D] = subfid2d(f, t, parms, Fmin, Fmax);
  else
    % throw an exception.
    error('subfid: input data must be one- or two-dimensional');
  end
end

