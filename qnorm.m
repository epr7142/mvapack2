## Copyright (C) 2017 University of Nebraska Board of Regents.
## Written by Thao Vu and Eli Riekeberg <eriekeberg@huskers.unl.edu>.
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

%Want to remove numbers from final functions. Succinctness.
## -*- texinfo -*-
## @anchor{qnorm1}
## @deftypefn {Function File} {@var{Xn} =} qnorm1 (@var{X})}
## Normalize the observations of a data matrix using quantile normalization (qnorm).
## @end deftypefn

function [Xnorm] = qnorm1(matrix);

# Can we add argument checking?


# check the size of matrix
dim = size(matrix);
nrow = dim(1);
ncol = dim(2);
# sort the input matrix column-wise
Xsort = sort(matrix);

# create unit vector to do projection on
D = []; 
# values of each element in D
D_value = 1/sqrt(ncol);
D(1:ncol) = D_value;

#Computing the scaling factor

scale_factor = Xsort * D';

# projection on D

proj = scale_factor * D;

#for loop to get the rankings

for i = 1:ncol
  sorted_col = sort(matrix(:,i));
  ranking = lookup(sorted_col, matrix(:,i));
  rankings(:,i) = ranking;
endfor

# return the normalized matrix 
Xnorm = proj(rankings);

endfunction
