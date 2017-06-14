## Copyright (C) 2014 University of Nebraska Board of Regents.
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
## @anchor{cvplsscores}
## @deftypefn {Function File} {@var{T} = } cvplsscores (@var{mdl})
## @deftypefnx {Function File} {@var{T} = } cvplsscores (@var{mdl}, @var{n})
## Returns the cross-validated scores from a PLS model in the form of an
## array. Each element of the array contains the reconstructed scores
## from a monte carlo leave-n-out cross validation. The length of the
## array equals the number of cross-validation iterations.
## @end deftypefn

function T = cvplsscores (mdl, n)
  % check for proper arguments.
  if (nargin < 2)
    % improper arguments. print the usage statement.
    print_usage();
  end

  % initialize the scores cell arrays.
  T = cell(mdl.N, 1);
  Ti = cell(mdl.cv.niter, 1);

  % loop for every cross-validation iteration.
  for i = 1 : mdl.cv.niter
    % re-initialize a cell array for every group in the iteration.
    clear Tg;
    Tg = cell(mdl.cv.ngroup, 1);

    % compute scores for every group in the iteration.
    for g = 1 : mdl.cv.ngroup
      % compute T, the scores of the current validation set.
      Tg{g} = mdl.CV{i}.Xv{g} * mdl.CV{i}.Ws{g};
    end

    % join the validation-set scores into a single matrix for the
    % current iteration.
    Ti{i} = cvjoin([], Tg, mdl.CV{i}.idx);
  end

  % loop for every observation.
  for obs = 1 : mdl.N
    % initialize a matrix for the current observation.
    T{obs} = [];

    % loop for every cross-validation iteration.
    for i = 1 : mdl.cv.niter
      % append the current iteration observation to the final matrix.
      T{obs} = [T{obs}; Ti{i}(obs, 1 : min([n; mdl.A]))];
    end
  end
end
