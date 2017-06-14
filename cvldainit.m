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
## @anchor{cvldainit}
## @deftypefn {Function File} {@var{CV} =} cvldainit (@var{mdl})
## Initializes and returns a cell array used for internal cross-validation
## of LDA models.
## @end deftypefn

function CV = cvldainit (mdl)
  % initialize the cross-validation cell array.
  CV = cell(mdl.cv.niter, 1);

  % loop through the number of cross-validation iterations.
  for i = 1 : mdl.cv.niter
    % run a check to ensure that no class is singled out in a partitioning.
    ok = false;
    while (!ok)
      % be really optimistic.
      ok = true;

      % initialize the cross-validation group indices.
      CV{i}.idx = cvindices(mdl.N, mdl.cv.ngroup);

      % split up the PLS E and F matrices for cross-validation, which are
      % initialized to X and Y at the time that this function will be called.
      [CV{i}.Xt, CV{i}.Xv] = cvsplit(mdl.E, CV{i}.idx);
      [CV{i}.Yt, CV{i}.Yv] = cvsplit(mdl.F, CV{i}.idx);

      % initialize the X-weights, loadings and Y-weights cell arrays.
      CV{i}.W = cell(mdl.cv.ngroup, 1);
      CV{i}.P = cell(mdl.cv.ngroup, 1);
      CV{i}.C = cell(mdl.cv.ngroup, 1);

      % initialize the estimated Y matrix.
      CV{i}.Yh = zeros(size(mdl.F));

      % check that no training set has a single class entry.
      for g = 1 : mdl.cv.ngroup
        % count the members of each class in the training partition.
        Yt = CV{i}.Yt{g};
        nt = arrayfun(@(m) length(find(Yt(:,m) == 1)), [1 : columns(Yt)]);

        % check if we have an issue.
        if (any(nt <= 1))
          % yeah, we have to remake this partitioning.
          ok = false;
          break;
        end
      end
    end
  end
end

