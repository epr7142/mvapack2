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
## @anchor{j2}
## @deftypefn {Function File} {@var{j2v} =} j2 (@var{mdl})
## @deftypefnx {Function File} {@var{j2v} =} j2 (@var{mdl}, @var{Y})
## Computes the @math{J_2} clustering statistic (ratio of the determinants
## of entire dataset covariance to each cluster covariance) on the
## computed scores of a PCA, PLS or OPLS model. PCA models need an
## accompanying @var{Y}-matrix to define classes. It may either be
## passed as a second argument to this function or added to the
## model (@xref{addclasses}).
##
## For more information on the @math{J_2} metric, see:
## K. Koutroumbas, S. Theodoridis. `Pattern Recognition'. Elsevier Press,
## Amsterdam, 2006.
## @quotation
## 
## @end quotation
## @end deftypefn

function j2v = j2 (mdl, Y)
  % check if the class matrix is provided.
  if (nargin < 2 || isempty(Y))
    % not provided. check if it should have been provided.
    if (!isfield(mdl, 'Y'))
      % yes. throw an exception.
      error('j2: class matrix (Y) required for unsupervised data');
    else
      % no. use the classes embedded in the model.
      Y = backscaleclasses(mdl);
    end
  end

  % extract the scores from the model.
  T = scores(mdl);

  % get the number of classes.
  M = columns(Y);

  % initialize the output j2 vector.
  j2v = [];

  % calculate the spread of the whole dataset.
  det_all = det(cov(mdl.T));

  % loop through the classes.
  for m = 1 : M
    % get the scores for the current class.
    idx = classidx(Y, m);

    % calculate the spread of the current class.
    det_grp = det(cov(mdl.T(idx,:)));

    % calculate the J2 statistic for the current class.
    j2v = [j2v; det_all / det_grp];
  end
end

