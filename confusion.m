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
## @anchor{confusion}
## @deftypefn {Function File} {@var{D} = } confusion (@var{mdl})
## @deftypefnx {Function File} {[@var{D}, @var{hits}] = } confusion (@var{mdl})
## @deftypefnx {Function File} {[@var{D}, @var{hits}, @var{misses}] = } confusion (@var{mdl})
## Returns the results of cross-validating a PLS-DA model in the format of
## a confusion matrix.
## @end deftypefn

function [D, hits, misses] = confusion (mdl)
  % check for proper arguments.
  if (nargin < 1)
    % improper arguments. print the usage statement.
    print_usage();
  end

  % make sure the model argument is a structure.
  if (!isstruct(mdl))
    % it isn't! throw an exception.
    error('confusion: model argument must be a structure type');
  end

  % check that the model has a cross-validation structure.
  if (!isfield(mdl, 'cv') || !isfield(mdl, 'CV') || isempty(mdl.CV))
    % throw an error.
    error('confusion: model has no cross-validation information');
  end

  % initialize the confusion matrix.
  D = zeros(mdl.M, mdl.M);

  % loop over all cross-validation iterations.
  for i = 1 : mdl.cv.niter
    % loop over all cross-validation groups.
    for g = 1 : mdl.cv.ngroup
      % backscale the predicted and true validation-set response matrices.
      Ybs = backscale(mdl.CV{i}.Yh{g}, mdl.mean.Y, mdl.scale.Y);
      Y = backscale(mdl.CV{i}.Yv{g}, mdl.mean.Y, mdl.scale.Y);

      % classify the validation-set observations based on predicted response.
      [dmax, Mmax] = max(Ybs, [], 2);
      Yhat = eye(mdl.M)(Mmax,:);

      % loop over the classified observations.
      for n = 1 : rows(Y)
        % get the confusion matrix indices to increment.
        a = find(Y(n,:));
        b = find(Yhat(n,:));

        % increment the appropriate matrix element.
        D(a,b) += 1;
      end
    end
  end

  % check if the number of hits was requested.
  if (nargout >= 2)
    % return the number of correct classifications (hits).
    hits = sum(diag(D));
  end

  % check if the number of misses was requested.
  if (nargout >= 3)
    % return the number of incorrect classifications (misses).
    misses = sum(vec(D - diag(diag(D))));
  end
end
