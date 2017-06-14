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
## @anchor{cvlda}
## @deftypefn {Function File} {[@var{Qsq}, @var{Qstd}, @var{CV}] =} cvlda (@var{mdl}, @var{Amax})
## Performs internal cross-validation of an LDA model and returns
## a @math{Q^2} value for inferring model reliability. The @var{CV}
## cell array returned is a modified version of that found in the
## passed model.
## @end deftypefn

function [Qsq, Qstd, CV] = cvlda (mdl, Amax)
  % initialize the output structures.
  Q = zeros(mdl.cv.niter, 1);
  Qsq = [];
  Qstd = [];
  CV = mdl.CV;

  % loop through the number of cross-validation iterations.
  for i = 1 : mdl.cv.niter
    % loop through the number of cross-validation groups.
    for g = 1 : mdl.cv.ngroup
      % calculate a new decomposition for cross-validation.
      [P, D, U] = ldacomp(CV{i}.Xt{g}, CV{i}.Yt{g});

      % store the decomposition.
      CV{i}.P{g} = P;
      CV{i}.D{g} = D;
      CV{i}.U{g} = U;
    end
  end

  % loop through all possible component counts.
  for a = 1 : Amax
    % compose a set of dimension-limited parameters.
    Pq = P(:,1:a);

    % loop through the number of cross-validation iterations.
    for i = 1 : mdl.cv.niter
      % loop through the number of cross-validation groups.
      for g = 1 : mdl.cv.ngroup
        % estimate the left-out responses.
        Yhat = ldaclassify(Pq, U, CV{i}.Xv{g});
        CV{i}.Yh(find(CV{i}.idx == g), :) = Yhat;
      end

      % calculate the predictive ability statistic.
      Q(i) = 1 - ssratio(mdl.Y - CV{i}.Yh, mdl.Y);
    end

    % compute the standard deviation of Qsq values over all iterations.
    Qstd = [Qstd; std(Q)];

    % determine how to report the Qsq values.
    if (a == 1)
      % this is the first component. return the value verbatim.
      Qsq = [Qsq; mean(Q)];
    else
      % this is a subsequent component. return the value with the previous
      % value subtracted out.
      Qsq = [Qsq; mean(Q) - sum(Qsq)];
    end
  end
end

