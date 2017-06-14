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
## @anchor{ransy1d}
## @deftypefn {Function File} {@var{R} =} ransy1d (@var{X}, @var{k})
## Use the Ratio Analysis SpectroscopY method (RANSY) to extract single
## compound spectra from spectra of complex mixtures, defined here:
##
## @quotation
## S. Wei, et. al. `Ratio Analysis Nuclear Magnetic Resonance Spectroscopy
## for Selective Metabolite Identification in Complex Samples'. Analytical
## Chemistry 2011(83): 7616-7623.
## @end quotation
## @end deftypefn

function R = ransy1d (X, k)
  % check for proper arguments.
  if (nargin != 2 || nargout != 1 || !ismatrix(X) || !isscalar(k))
    % improper arguments. print the usage statement.
    print_usage();
  end

  % get the dimensions of the input matrix.
  [n, m] = size(X);

  % run the ransy calculation.
  D = X ./ (X(:,k) * ones(1, m));
  R = mean(D) ./ std(D);

  % remove infinities from the ransy matrix.
  R(k) = max(R(R != inf));
end

