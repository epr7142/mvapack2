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
## @anchor{cvsplit}
## @deftypefn {Function File} {[@var{Xt}, @var{Xv}] =} cvsplit (@var{X}, @var{idx})
## Creates a one-dimensional array of datasets from a single dataset,
## where the array of indices @var{idx} was generated using the appropriate
## functions (@xref{cvindices}).
## @end deftypefn

function [Xt, Xv] = cvsplit (X, idx)
  % initialize the output cell arrays.
  Xt = {};
  Xv = {};

  % loop through the cross-validation indices.
  for g = min(idx) : max(idx)
    % build the cross-validation sub-matrix for the current CV group.
    Xt{g} = X(idx != g,:);
    Xv{g} = X(idx == g,:);
  end
end

