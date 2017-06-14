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
## @anchor{stocsy}
## @deftypefn {Function File} {@var{C} =} stocsy (@var{X}, @var{ab})
## Use the Statistical Total Correlation SpectroscopY (STOCSY) method to
## generate a contour map of correlations between spectral variables in
## a 1D NMR dataset, defined here:
##
## @quotation
## O. Cloarec, et. al. `Statistical Total Correlation Spectroscopy: An
## Exploratory Approach for Latent Biomarker Identification from Metabolic
## 1H NMR Data Sets'. Analytical Chemistry 2005(77): 1282-1289.
## @end quotation
## @end deftypefn

function C = stocsy (X, ab, n)
  % check for proper arguments.
  if (nargin < 2 || !ismatrix(X) || !isvector(ab))
    % improper arguments. print the usage statement.
    print_usage();
  end

  % see if a level count was requested.
  if (nargin >= 3 && !isempty(n))
    % yes, check the type of the level argument.
    if (isscalar(n))
      % build a list of levels than spans the range.
      lev = sort((linspace(min(vec(X)), max(vec(X)), n).^2) ./ (rows(X) - 1));
    elseif (isvector(n))
      % use the list of levels requested.
      lev = unique(sort(n));
    else
      % invalid type.
      error('stocsy: level argument must be a scalar or vector type');
    end
  else
    % the default number of contours is ten.
    lev = sort((linspace(min(vec(X)), max(vec(X)), 10).^2) ./ (rows(X) - 1));
  end

  % turn the levels vector into a column vector.
  lev = reshape(lev, length(lev), 1);
  ab = reshape(ab, length(ab), 1);

  % run the stocsy calculation.
  C = __stocsy(X, ab, lev);
end

