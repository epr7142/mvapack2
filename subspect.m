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
## @anchor{subspect}
## @deftypefn {Function File} {[@var{ssub}, @var{absub}] =} subspect (@var{s}, @var{ab}, @var{parms}, @var{Fmin}, @var{Fmax})
## Extracts a band of frequencies (@math{[F_{min},F_{max}]}) from @var{s} into
## @var{ssub}. This function can operate on multiple spectra at once if @var{s}
## is supplied as a data matrix. The function also accepts two-dimensional
## spectral data in the form of matrices or cell arrays.
##
## The values in @var{Fmin} and @var{Fmax} should correspond to those in
## @var{ab}. For one-dimensional data, @var{Fmin} and @var{Fmax} must be
## scalar values. For two-dimensional data, they must be two-vectors
## containing frequencies for each dimension.
##
## Unlike @ref{subfid}, which performs its extraction in the time domain,
## this function has it easy: it just crops out the subspectrum of interest
## from the input spectral data.
## @end deftypefn

function [ssub, absub] = subspect (s, ab, parms, Fmin, Fmax)
  % check if the number of expected arguments was passed.
  if (nargin != 5 || nargout != 2)
    % print the usage statement.
    print_usage();
  end

  % get the dimensions in the dataset.
  nd = nmrdims(s, parms);

  % check whether we are cropping one- or two-dimensional data.
  if (nd == 1)
    % see if we need to reduce the parms cell array to a struct.
    if (iscell(parms))
      parms = parms{1};
    end

    % apply a one-dimensional crop.
    [ssub, absub] = subspect1d(realnmr(s, parms), ab, parms, Fmin, Fmax);
  elseif (nd == 2)
    % apply a two-dimensional filter.
    [ssub, absub] = subspect2d(realnmr(s, parms), ab, parms, Fmin, Fmax);
  else
    % throw an exception.
    error('subspect: input data must be one- or two-dimensional');
  end
end

