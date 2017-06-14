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
## @anchor{plsroc}
## @deftypefn {Function File} {[@var{x}, @var{y}] = } plsroc (@var{mdl})
## @deftypefnx {Function File} {[@var{x}, @var{y}, @var{auc}] = } plsroc (@var{mdl})
## Returns the results of cross-validating a two-class PLS-DA model in the
## format of a receiver operating characteristic (ROC) curve.
## @end deftypefn

function [x, y, auc] = plsroc (mdl)
  % check for proper arguments.
  if (nargin < 1 || nargout < 2)
    % improper arguments. print the usage statement.
    print_usage();
  end

  % set the number of density axis values.
  nth = 256;

  % make sure the model argument is a structure.
  if (!isstruct(mdl))
    % it isn't! throw an exception.
    error('plsroc: model argument must be a structure type');
  end

  % check that the model has a cross-validation structure.
  if (!isfield(mdl, 'cv') || !isfield(mdl, 'CV') || isempty(mdl.CV))
    % throw an error.
    error('plsroc: model has no cross-validation information');
  end

  % check that the model contains two classes.
  if (mdl.M != 2)
    % throw an error.
    error('plsroc: only two-class pls-da models are supported');
  end

  % initialize the concatenated responses.
  Yh = [];
  Yv = [];

  % loop over all cross-validation iterations.
  for i = 1 : mdl.cv.niter
    % extract the true and predicted validation-set response matrices.
    Yh = [Yh; cvjoin([], mdl.CV{i}.Yh, mdl.CV{i}.idx)];
    Yv = [Yv; cvjoin([], mdl.CV{i}.Yv, mdl.CV{i}.idx)];
  end

  % backscale the extracted matrices.
  Yh = backscale(Yh, mdl.mean.Y, mdl.scale.Y);
  Yv = backscale(Yv, mdl.mean.Y, mdl.scale.Y);

  % grab the first column of the matrices.
  Yh = Yh(:,1);
  Yv = Yv(:,1);

  % initialize the output vectors.
  x = ones(nth, 1);
  y = zeros(nth, 1);

  % build a vector of threshold values.
  th = linspace(min(Yh), max(Yh), nth + 2)(2:end-1)';

  % loop over the threshold values.
  for i = 1 : nth
    % get the threshold value.
    t = th(i);

    % threshold the predicted y-matrix.
    Yth = (Yh > t);

    % get the number of positives and negatives.
    P = sum(Yth == 1);
    N = sum(Yth == 0);

    % get the number of true positives.
    TP = sum((Yth == Yv) .* (Yth == 1));

    % get the number of true negatives.
    TN = sum((Yth == Yv) .* (Yth == 0));

    % compute the true positive rate.
    if (P != 0)
      % only compute for nonzero P. if there are no positives, there are also
      % no true positives, so the result will remain zero.
      y(i) = TP / P;
    end

    % compute the false positive rate.
    if (N != 0)
      % only compute for nonzero N. if there are no negatives, there are also
      % no true negatives, so the result will remain one.
      x(i) = 1 - TN / N;
    end
  end

  % ensure the curve begins at (0,0) and ends at (1,1).
  x = [0; x; 1];
  y = [0; y; 1];

  % check if the user requested the area under the curve.
  if (nargout >= 3)
    % compute the area.
    auc = trapz(x, y);
  end
end
