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
## @anchor{peakpick}
## @deftypefn {Function File} {@var{P} =} peakpick (@var{s}, @var{ab})
## Picks peaks in a 1D spectrum vector or 2D spectrum matrix. The input
## spectrum @var{s} and abscissa @var{ab} are both required.
##
## The one-dimensional peak-picking algorithm is an implementation of:
##
## @quotation
## Du et. al., `Improved peak detection in mass spectrum by incorporating
## continuous wavelet transform-based pattern matching', Bioinformatics, 2006.
## @end quotation
##
## The two-dimensional peak-picking algorithm is not yet implemented.
## @end deftypefn

function P = peakpick (s, ab)
  % check if the number of expected arguments was passed.
  if (!any(nargin == [1 : 2]) || nargout != 1)
    % print the usage statement.
    print_usage();
  end

  % check if the spectrum contains complex values.
  if (iscomplex(s))
    % use only the real portion of the complex spectrum.
    s = real(s);
  end

  % check if an abscissa was provided.
  if (nargin >= 2 && !isempty(ab))
    % check the vector.
    if (!isvector(ab) || !isreal(ab) || length(ab) != length(s))
      % invalid abiscssa.
      error('peakpick: abscissa must be a real vector');
    end
  else
    % build a fake one.
    ab = linspace(-1, 1, length(s))';
  end

  % check the type of the input data.
  if (isvector(s))
    % estimate the noise floor of the spectrum.
    [nm, ns] = __estnoise1d(s);
    nf = nm + 2 * ns;

    % calculate the cwt of the spectrum.
    C = cwt(s, [], [], [], [], false);

    % calculate the fourier derivative coefficients.
    vd = -2 .* pi .* i .* linspace(-1, 1, rows(C))';
    vdd = vd .^ 2;

    % calculate the derivatives of the cwt.
    Cd = (vd * ones(1, columns(C))) .* C;
    Cdd = (vdd * ones(1, columns(C))) .* C;

    % inverse fourier transform all the cwt matrices.
    C = real(ifft(C));
    Cd = real(ifft(Cd));
    Cdd = real(ifft(Cdd));

    % run the compiled subroutine for one-dimensional peak picking.
    P = __peakpick1d(s, ab, nf, C, Cd, Cdd);
  elseif (ismatrix(s))
    % FIXME: functionality not yet implemented.
    error('peakpick: functionality (2D peak picking) not implemented');
  else
    % throw an exception.
    error('peakpick: input data must be a vector or matrix');
  end
end

