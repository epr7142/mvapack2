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
## @anchor{rmroi}
## @deftypefn {Function File} {@var{roirm} =} rmroi (@var{roi}, @var{rmzones})
## Removes any regions of interest from @var{roi} that overlaps the regions of
## interest specified in @var{rmzones}.
## @end deftypefn

function roirm = rmroi (roi, rmzones)
  % check for proper arguments.
  if (nargin != 2 || nargout != 1 || ...
      !ismatrix(roi) || !ismatrix(rmzones))
    % improper arguments. print the usage statement.
    print_usage();
  end

  % ensure the roi matrix and the zone matrix match.
  if (columns(roi) != columns(rmzones))
    % throw an exception.
    error('rmroi: region of interest column counts do not match');
  end

  % initialize the list of regions to remove.
  idx = [];

  % check whether the roi matrix is 1d or 2d.
  if (columns(roi) == 2)
    % pull out convenience vectors from the roi matrix.
    L = min(roi, [], 2);
    U = max(roi, [], 2);

    % loop through the removal zones.
    for z = 1 : rows(rmzones)
      % pull out the zone boundaries.
      zL = min(rmzones(z,:));
      zU = max(rmzones(z,:));

      % find all overlapping regions of interest.
      iadd = (L >= zL) .* (L <= zU) + ...
             (U >= zL) .* (U <= zU) + ...
             (L <= zL) .* (U >= zU);

      % add the identified regions.
      idx = [idx; find(iadd)];
    end
  elseif (columns(roi) == 4)
    % pull out convenience vectors from the roi matrix.
    L1 = roi(:, 1);
    U1 = roi(:, 2);
    L2 = roi(:, 3);
    U2 = roi(:, 4);

    % loop through the removal zones.
    for z = 1 : rows(rmzones)
      % pull out the zone boundaries.
      zL1 = rmzones(z, 1);
      zU1 = rmzones(z, 2);
      zL2 = rmzones(z, 3);
      zU2 = rmzones(z, 4);

      % find all overlapping regions of interest in the first dim.
      iadd1 = (L1 >= zL1) .* (L1 <= zU1) + ...
              (U1 >= zL1) .* (U1 <= zU1) + ...
              (L1 <= zL1) .* (U1 >= zU1);

      % find all overlapping regions of interest in the second dim.
      iadd2 = (L2 >= zL2) .* (L2 <= zU2) + ...
              (U2 >= zL2) .* (U2 <= zU2) + ...
              (L2 <= zL2) .* (U2 >= zU2);

      % compute the final overlap vector.
      iadd = iadd1 .* iadd2;

      % add the identified regions.
      idx = [idx; find(iadd)];
    end
  else
    % throw an exception.
    error('rmroi: regions of interest must be one- or two-dimensional');
  end

  % get a list without duplicated regions.
  idx = unique(idx);

  % remove the overlapped regions from the matrix.
  roirm = roi;
  roirm(sort(idx), :) = [];
end

