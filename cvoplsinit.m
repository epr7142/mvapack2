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
## @anchor{cvoplsinit}
## @deftypefn {Function File} {@var{CV} =} cvoplsinit (@var{mdl})
## Initializes and returns a cell array used for internal cross-validation
## of OPLS models.
## @end deftypefn

function CV = cvoplsinit (mdl)
  % initialize the cross-validation cell array.
  CV = cell(mdl.cv.niter, 1);

  % loop through the number of cross-validation iterations.
  for i = 1 : mdl.cv.niter
    % initialize the cross-validation group indices.
    CV{i}.idx = cvindices(mdl.N, mdl.cv.ngroup);

    % split up the PLS E and F matrices for cross-validation, which are
    % initialized to X and Y at the time that this function will be called.
    [CV{i}.Xt, CV{i}.Xv] = cvsplit(mdl.E, CV{i}.idx);
    [CV{i}.Yt, CV{i}.Yv] = cvsplit(mdl.F, CV{i}.idx);

    % initialize the X-weights, loadings and Y-weights cell arrays.
    CV{i}.Wp = cell(mdl.cv.ngroup, 1);
    CV{i}.Wo = cell(mdl.cv.ngroup, 1);
    CV{i}.Pp = cell(mdl.cv.ngroup, 1);
    CV{i}.Po = cell(mdl.cv.ngroup, 1);
    CV{i}.C = cell(mdl.cv.ngroup, 1);

    % initialize the estimated Y matrix.
    CV{i}.Yh = zeros(size(mdl.F));
  end
end

