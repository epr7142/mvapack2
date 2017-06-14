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
## @anchor{nmrist}
## @deftypefn {Function File} {[@var{recfid}] =} nmrist (@var{fid}, @var{parms})
## @deftypefnx {Function File} {[@var{recfid}] =} nmrist (@var{fid}, @var{parms}, @var{phc})
## Performs Iterative Soft Thresholding (IST) reconstruction of a nonuniformly
## sampled 2D NMR time domain matrix or a 2D NMR time domain cell array.
## Two-dimensional data must be arranged with slices of the direct-dimension
## along the rows.
##
## The product of this operation is a uniformly sampled data structure that
## may be processed in the same way as any other 2D data.
##
## An optional third argument, @var{phc}, may be passed as a two-element vector
## that supplies manual phase correction values for the direct dimension. If
## @var{phc} is not provided or is empty, automatic phase correction will be
## applied during reconstruction.
## @end deftypefn

function recfid = nmrist (fid, parms, phc)
  % check if the number of expected arguments was passed.
  if (!any(nargin == [2 : 3]) || nargout != 1)
    % print the usage statement.
    print_usage();
  end

  % check the parameter structure.
  if (!iscell(parms) || length(parms) != 2)
    % throw an exception.
    error('nmrist: parameters must be a two-element cell array');
  end

  % check the phase correction vector.
  if (nargin < 3 || isempty(phc))
    % set an empty phase correction vector.
    phc = [];
  end

  % check the type of the input data.
  if (ismatrix(fid))
    % fourier transform the rows of the input matrix.
    % (this takes care of the direct dimension)
    s = fft(fid, [], 2);

    % shift the rows of the input matrix.
    v = round(columns(s) / 2);
    s = shift(s, v, 2);

    % check if a two-element phase correction vector was provided.
    if (isempty(phc) || length(phc) != 2)
      % run automatic phase correction of the first cosine-modulated slice.
      [junk, phc0, phc1] = autophase1d(s(2,:)', @simplex_entropy);
      phc = [phc0, phc1];
    end

    % phase-correct the direct dimension.
    s = phase(s, parms, [phc(1), 0], [phc(2), 0]);

    % de-interlace the half-transformed matrix.
    [Anus, Bnus] = states(s);

    % initialize the internally zero-filled uniform states matrices.
    A = zeros(parms{2}.tdnus / 2, columns(s));
    B = zeros(parms{2}.tdnus / 2, columns(s));

    % fill the zero matrices with sampled values.
    A(parms{2}.sched + 1, :) = Anus;
    B(parms{2}.sched + 1, :) = Bnus;

    % reconstruct the columns of the states matrices.
    A = ist(A, parms{2}.sched + 1);
    B = ist(B, parms{2}.sched + 1);

    % re-interlace the fully reconstructed matrix.
    s = states(A, B);

    % de-shift the rows of the input matrix and inverse ft.
    s = shift(s, -v, 2);
    recfid = ifft(s, [], 2);
  elseif (iscell(fid))
    % initialize the output cell array.
    recfid = cell(size(fid));
    
    % fourier transform each matrix in the cell array.
    for idx = 1 : length(fid)
      % reconstruct the matrix.
      recfid{idx} = nmrist(fid{idx}, parms, phc);
    end
  else
    % throw an exception.
    error('nmrist: input data must be a matrix or cell array');
  end
end

