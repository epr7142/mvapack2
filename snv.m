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
## @anchor{snv}
## @deftypefn {Function File} {@var{Xn} =} snv (@var{X})
## @deftypefnx {Function File} {[@var{Xn}, @var{mu}] =} snv (@var{X})
## @deftypefnx {Function File} {[@var{Xn}, @var{mu}, @var{s}] =} snv (@var{X})
## Normalize the observations of a data matrix using standard normal variate
## normalization (SNV). Row centering is also performed in the process. The
## variables used to center and scale @var{X} may optionally
## be returned.
## @end deftypefn

function [Xnew, mu, s] = snv (X)
  % check for proper arguments.
  if (nargin != 1 || !ismatrix(X))
    % improper arguments. print the usage statement.
    print_usage();
  end

  % return the normalized output matrix.
  Xnew = diag(1 ./ std(X, 0, 2)) * center(X, 2);

  % see if the center was requested.
  if (nargout >= 2)
    % yes. return the center.
    mu = mean(X, 2);
  end

  % see if the scale was requested.
  if (nargout >= 3)
    % yes. return the scale.
    s = std(X, 0, 2);
  end
end

