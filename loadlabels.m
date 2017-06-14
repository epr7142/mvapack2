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
## @anchor{loadlabels}
## @deftypefn {Function File} {@var{labels} =} loadlabels (@var{fname})
## @deftypefnx {Function File} {[@var{labels}, @var{indices}] =} loadlabels (@var{fname})
## @deftypefnx {Function File} {[@var{labels}, @var{indices}, @var{Y}] =} loadlabels (@var{fname})
## Loads in class label assignments for a dataset from a text file, where the
## n'th line in the text file contains the class label of the n'th observation.
##
## An optional second return value (@var{indices}) can be requested that
## contains the row indices that will bring the associated data matrix
## into sync with the class membership matrix, like so:
## @code{X = X(indices,:);}
##
## An optional third return value (@var{Y}) can be requested that contains the
## binary class membership matrix used during PLS.
## @end deftypefn

function [labels, indices, Y] = loadlabels (fname)
  % check for proper arguments.
  if (nargin != 1 || !any(nargout == [1 : 3]) || !ischar(fname))
    % improper arguments. print the usage statement.
    print_usage();
  end

  % open the input file.
  fh = fopen(fname, "rb");

  % make sure the file was opened successfully.
  if (fh == -1)
    % we failed. return an error.
    error('loadlabels: file not found: "%s"', fname);
  end

  % initialize the strings and labels arrays.
  strings = {};
  labels = {};

  % loop through the file until we've reached the end.
  while (!feof(fh))
    % get the next line of the file.
    str = fgetl(fh);

    % append the string to the strings array.
    strings = cat(2, strings, str);
  end

  % build the labels array as the unique set of strings.
  [jnk, idx] = sort(strings);
  labels = sort(unique(strings));

  % close the input file.
  fclose(fh);

  % see if indices were requested.
  if (nargout >= 2)
    % return the indices.
    indices = idx;
  end

  % see if classes were requested.
  if (nargout >= 3)
    % initialize a count vector.
    counts = zeros(1, length(labels));

    % loop through the labels.
    for i = 1 : length(labels)
      % store the count for the current label.
      counts(i) = sum(arrayfun(@(x,y) strcmp(x,y), strings, labels(i)));
    end

    % return the classes.
    Y = classes(counts);
  end
end

