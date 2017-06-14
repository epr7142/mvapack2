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
## @anchor{rmnoise}
## @deftypefn {Function File} {[@var{Xrm}, @var{abrm}] =} rmnoise (@var{X}, @var{ab}, @var{idx})
## @deftypefnx {Function File} {[@var{Xrm}, @var{abrm}] =} rmnoise (@var{X}, @var{ab}, @var{idx}, @var{nstd})
## @deftypefnx {Function File} {[@var{Xrm}, @var{abrm}, @var{idxrm}] =} rmnoise (@var{X}, @var{ab}, @var{idx})
## @deftypefnx {Function File} {[@var{Xrm}, @var{abrm}, @var{idxrm}] =} rmnoise (@var{X}, @var{ab}, @var{idx}, @var{nstd})
## Uses a relative standard deviation to distinguish between signal and noise
## in a binned NMR spectral dataset. Bins identified as noise will be removed
## in the output data matrix @var{Xrm} and abscissa @var{abrm}. The @var{idx}
## variable is used to defined a noise region.
##
## An optional argument @var{nstd} may be supplied (default: 2) to set how many
## standard deviations from the mean noise value the threshold will be.
##
## @strong{IMPORTANT:} This routine is @emph{not} designed to be used with
## full-resolution data matrices, @emph{especially} when binning is a
## subsequent step. Such actions are almost surely guaranteed to yield
## strange results.
## @end deftypefn

function [Xrm, abrm, idxrm] = rmnoise (X, ab, idx, nstd)
  % check for proper arguments.
  if (!any(nargin == [3 : 4]) || !any(nargout == [2 : 3]) || ...
      !ismatrix(X) || !isvector(ab) || !isvector(idx))
    % improper arguments. print the usage statement.
    print_usage();
  end

  % see if a fourth argument was supplied.
  if (nargin < 4 || isempty(nstd))
    % use the default value: two standard deviations.
    nstd = 2;
  else
    % ensure the argument is a scalar.
    if (!isscalar(nstd) || !isreal(nstd))
      % not scalar. throw an exception.
      error('rmnoise: nstd argument must be a real scalar');
    end
  end

  % build a standard normal variate copy of the data matrix.
  Z = snv(X);

  % calculate statistics over the variables.
  Zmax = max(Z);
  Zstd = std(Z);
  Zmean = mean(Z);

  % build the threshold rsd value from the noise indices.
  rsd = abs(Zstd ./ Zmean);
  rsdTh = mean(rsd(idx)) + nstd * std(rsd(idx));

  % identify the variables that have rsd values lower than the threshold
  % and also have negative z-values.
  idxrm_rsd = find(rsd <= rsdTh);
  idxrm_max = find(Zmax <= 0);
  idxrm = intersect(idxrm_rsd, idxrm_max);

  % remove the identified variables.
  [Xrm, abrm] = rmvar(X, ab, idxrm);
end

