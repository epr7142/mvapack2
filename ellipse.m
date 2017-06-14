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
## @anchor{ellipse}
## @deftypefn {Function File} {@var{E} =} ellipse (@var{X}, @var{correlated})
## Assuming the input data rows in @var{X} are normally distributed,
## calculates the @code{alpha = 0.05} confidence ellipse around the
## points (rows) in @var{X}. A second optional argument, @var{correlated},
## may be passed to use either a diagonal or full covariance matrix. The
## default behavior is to use a full (correlated) matrix.
## @end deftypefn

function E = ellipse (X, correlated = true)
  % get the dimensions of the input matrix.
  [N, K] = size(X);

  % check if the dimensions are expected to have correlated values.
  if (correlated == true)
    % yes. run a full covariance analysis.
    C = cov(X);
  else
    % no. calculate only the variances.
    C = diag(var(X));
  end

  % calculate the eigenvalue-sorted eigendecomposition of the covariances.
  [V, lambda] = eigsort(C);

  % store the critical value of the chi-squared distribution at 0.95 for
  % K degrees of freedom.
  F = chi2inv(0.95, K);

  % store the number of ellipse dimensions.
  E.dims = K;

  % calculate the center of the ellipse.
  E.center = mean(X)';

  % calculate the radius of the ellipse. this assumes an underlying normal
  % distribution and sets the ellipse radius at two standard deviations from
  % the mean value.
  E.radii = sqrt(lambda .* F);

  % store the rotation of the ellipse.
  E.rots = V;

  % see if a two-dimensional ellipse was requested.
  if (K == 2)
    % build a parameter vector.
    t = [0 : pi / 128 : 2 * pi];

    % build a circle matrix. :)
    cst = [cos(t); sin(t)];

    % build a matrix of offsets (the centering matrix).
    cents = E.center * ones(size(t));

    % calculate the final cartesian coordinates of the ellipse points.
    E.xy = (cents + E.rots * diag(E.radii) * cst)';
  elseif (K == 3)
    % no, but a three-dimensional one was. build the mesh.
    [E.x, E.y, E.z] = ellipsoid(0, 0, 0, 1, 1, 1, 10);

    % scale the mesh coordinates.
    E.x = E.x .* E.radii(1);
    E.y = E.y .* E.radii(2);
    E.z = E.z .* E.radii(3);

    % rotate the mesh coordinates into temporary matrices.
    xtmp = arrayfun(@(x, y, z) V(1,:) * [x; y; z], E.x, E.y, E.z);
    ytmp = arrayfun(@(x, y, z) V(2,:) * [x; y; z], E.x, E.y, E.z);
    ztmp = arrayfun(@(x, y, z) V(3,:) * [x; y; z], E.x, E.y, E.z);

    % replace the temporary values with the rotated values.
    E.x = xtmp;
    E.y = ytmp;
    E.z = ztmp;

    % translate the mesh coordinates.
    E.x += E.center(1);
    E.y += E.center(2);
    E.z += E.center(3);
  end
end

