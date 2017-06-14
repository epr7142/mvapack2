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
## @anchor{states}
## @deftypefn {Function File} {[@var{A}, @var{B}] =} states (@var{X})
## @deftypefnx {Function File} {@var{A} =} states (@var{X})
## @deftypefnx {Function File} {@var{X} =} states (@var{A}, @var{B})
## De-interlaces States/Haberkorn/Ruben cosine- and sine-modulated rows
## of a complex matrix @var{X} into the complex matrices @var{A} and @var{B}.
## @tex
##
## The de-interlacing procedure is as follows:
## $$ X_{cos} \leftarrow X_{2:2:N,\cdot} $$
## $$ X_{sin} \leftarrow X_{1:2:N-1,\cdot} $$
## $$ A \leftarrow Re\{X_{cos}\} + i \cdot Re\{X_{sin}\} $$
## $$ B \leftarrow Im\{X_{cos}\} + i \cdot Im\{X_{sin}\} $$
## @end tex
##
## Alternatively, this function can re-interlace the matrices @var{A} and
## @var{B} to yield a new original matrix @var{X}.
## @tex
##
## The re-interlacing procedure is as follows:
## $$ X_{cos} \leftarrow Re\{A\} + i \cdot Re\{B\} $$
## $$ X_{sin} \leftarrow Im\{A\} + i \cdot Im\{B\} $$
## $$ X_{2:2:N,\cdot} \leftarrow X_{cos} $$
## $$ X_{1:2:N-1,\cdot} \leftarrow X_{sin} $$
## @end tex
## @end deftypefn

function [out1, out2] = states (in1, in2)
  % check for proper arguments.
  if (!any(nargin == [1, 2]))
    % improper arguments. print the usage statement.
    print_usage();
  end

  % check if we are de-interlacing or re-interlacing.
  if (nargin == 1 && any(nargout == [1, 2]))
    % get the input argument.
    X = in1;

    % extract the sine and cosine submatrices.
    Xcos = X(2 : 2 : end, :);
    Xsin = X(1 : 2 : end - 1, :);

    % build the output matrices.
    A = complex(real(Xcos), real(Xsin));
    B = complex(imag(Xcos), imag(Xsin));

    % return the output arguments.
    out1 = A;
    out2 = B;
  elseif (nargin == 2 && nargout == 1)
    % get the input arguments.
    A = in1;
    B = in2;

    % build the sine and cosine submatrices.
    Xcos = complex(real(A), real(B));
    Xsin = complex(imag(A), imag(B));

    % build the output matrix.
    X = zeros(rows(A) + rows(B), columns(A));
    X(2 : 2 : end, :) = Xcos;
    X(1 : 2 : end - 1, :) = Xsin;

    % return the output argument.
    out1 = X;
  else
    % invalid combo. throw an exception.
    error('states: invalid input/output combination');
  end
end

