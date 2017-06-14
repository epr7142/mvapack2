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
## @anchor{classes}
## @deftypefn {Function File} {@var{Y} =} classes ([@var{n1}, @var{n2}, @dots{}, @var{nM}])
## Creates a Y-matrix suitable for discriminant analysis, where class
## membership is denoted by a 1 in the corresponding class column. One
## vector argument is required, where the number of elements in the vector
## equals the number of classes and each element in the vector gives
## the number of observations in that class.
## @end deftypefn

function Y = classes (counts)
  % make sure the classes function was given a proper argument.
  if (nargin != 1 || nargout != 1 || !isvector(counts))
    % print the function usage.
    print_usage();
  end

  % check that the counts argument has at least two entries.
  if (length(counts) < 2)
    % throw an exception.
    error('classes: input counts vector must have at least two entries');
  end

  % the sum of the class sizes should equal the total number of observations.
  N = sum(counts);

  % the number of class sizes should equal the number of classes.
  M = length(counts);

  % build an initial matrix of zeros for the class matrix.
  Y = zeros(N, M);

  % build a matrix of observation indices where each class
  % begins (column 1) and ends (column 2).
  k = [cumsum(counts) - counts + 1; cumsum(counts)]';

  % loop through the classes.
  for i = 1 : M
    % extract the relevant row of the index matrix.
    krow = k(i,:);

    % insert ones into the relevant fields of the class matrix (Y).
    Y(krow(1):krow(2),i) = ones(counts(i), 1);
  end
end

