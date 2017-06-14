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
## @anchor{cvpcainit}
## @deftypefn {Function File} {@var{CV} =} cvpcainit (@var{mdl})
## Initializes and returns a cell array used for internal cross-validation
## of PCA models.
## @end deftypefn

function CV = cvpcainit (mdl)
  % initialize the cross-validation cell array.
  CV = cell(mdl.cv.niter, 1);

  % loop through the number of cross-validation iterations.
  for i = 1 : mdl.cv.niter
    % initialize the cross-validation group indices.
    CV{i}.nidx = cvindices(mdl.N, mdl.cv.ngroup);
    CV{i}.kidx = cvindices(mdl.K, mdl.cv.ngroup);

    % allocate arrays for the group index sets.
    CV{i}.nset = cell(mdl.cv.ngroup, 1);
    CV{i}.kset = cell(mdl.cv.ngroup, 1);

    % allocate arrays for the group masks.
    CV{i}.nmask = cell(mdl.cv.ngroup, 1);
    CV{i}.kmask = cell(mdl.cv.ngroup, 1);

    % loop through the number of cross-validation groups.
    for g = 1 : mdl.cv.ngroup
      % compute the index set for the group.
      CV{i}.nset{g} = find(CV{i}.nidx == g);
      CV{i}.kset{g} = find(CV{i}.kidx == g);

      % compute the masks for the index set.
      CV{i}.nmask{g} = (CV{i}.nidx == g);
      CV{i}.kmask{g} = (CV{i}.kidx == g);
    end
  end
end

