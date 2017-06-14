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
## @anchor{cvpca}
## @deftypefn {Function File} {[@var{Qsq}, @var{Qstd}, @var{CV}] =} cvpca (@var{mdl}, @var{t}, @var{p})
## Performs internal cross-validation of a PCA model and returns
## a @math{Q^2} value for inferring model reliability.
## @end deftypefn

function [Qsq, Qstd, CV] = cvpca (mdl, t, p)
  % initialize the output structures.
  Q = zeros(mdl.cv.niter, 1);
  CV = mdl.CV;

  % get simpler handles on a few things.
  N = mdl.N;
  K = mdl.K;

  % loop through the number of cross-validation iterations.
  for i = 1 : mdl.cv.niter
    % initialize the estimated data matrix.
    Ehat = zeros(size(mdl.E));

    % loop through the number of cross-validation groups.
    for g = 1 : mdl.cv.ngroup
      % row-reduce the scores.
      tr = t;
      tr(CV{i}.nset{g}) = [];

      % extract the row-reduced submatrix.
      Er = mdl.E;
      Er(CV{i}.nset{g}, :) = [];

      % calculate a new component loading.
      [jnk, pcv, iter] = __pcacomp(Er, tr);

      % run a parity check on the loadings.
      if (p' * -pcv > p' * pcv)
        pcv = -pcv;
      end

      % loop again through the groups.
      for h = 1 : mdl.cv.ngroup
        % get the number of left-out columns.
        Kout = length(CV{i}.kset{h});

        % extract the column-reduced submatrix.
        Ec = mdl.E;
        Ec(:, CV{i}.kset{h}) = [];

        % calculate a new component score.
        [tcv, jnk, iter] = __pcacomp(Ec, t);

        % adjust the component score.
        tcv = tcv .* sqrt(K / (K - Kout));

        % run a parity check on the scores.
        if (t' * -tcv > t' * tcv)
          tcv = -tcv;
        end

        % store the estimated matrix elements.
        Ehat += (CV{i}.nmask{g} .* tcv) * (CV{i}.kmask{h} .* pcv)';
      end

      % clear the temporary submatrices.
      clear Er Ec;
    end

    % calculate Q, one minus the ratio of the error and true sum of squares.
    Q(i) = 1 - ssratio(mdl.X - Ehat, mdl.X);

    % clear the estimated data matrix.
    clear Ehat;
  end

  % compute the mean and standard deviation of the iterations' Qsq values.
  Qsq = mean(Q);
  Qstd = std(Q);
end

