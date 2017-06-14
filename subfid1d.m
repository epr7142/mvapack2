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
## @anchor{subfid1d}
## @deftypefn {Function File} {[@var{fsub}, @var{tsub}, @var{dfbw}, @var{D}] =} subfid1d (@var{f}, @var{t}, @var{parms}, @var{Fmin}, @var{Fmax})
## Extracts a band of frequencies (@math{[F_{min},F_{max}]}) from @var{f} into
## @var{fsub}.
##
## It is highly recommended that you use @ref{subfid} instead of calling this
## function directly.
## @end deftypefn

function [fsub, tsub, dfbw, D] = subfid1d (f, t, parms, Fmin, Fmax)
  % check if the number of expected arguments was passed.
  if (nargin != 5 || nargout != 4)
    % print the usage statement.
    print_usage();
  end

  % check the decay data argument.
  if (!(isvector(f) || ismatrix(f)) || !iscomplex(f))
    % invalid type. throw an exception.
    error('subfid1d: decay data must be a complex vector or a matrix.');
  end

  % check the time axis argument.
  if (!isvector(t) || !isreal(t))
    % invalid type. throw an exception.
    error('subfid1d: time axis must be a real vector.');
  end

  % check the parameter argument.
  if (!isstruct(parms))
    % invalid type. throw an exception.
    error('subfid1d: parameters must be a structure');
  end

  % check the frequency boundary arguments.
  if (!isscalar(Fmin) || !isscalar(Fmax) || !isreal(Fmin) || !isreal(Fmax))
    % invalid type. throw an exception.
    error('subfid1d: frequency bounds must be real scalar values.');
  end

  % initialize some default options.
  % @expand: fractional amount to expand the FIR filter bandwidth.
  % @N: number of taps to use for the FIR filter.
  expand = 0.1;
  N = 512;

  % calculate the window center frequency and bandwidth.
  F0 = (Fmin + Fmax) / 2;
  Fw = abs(Fmax - Fmin);

  % calculate the digital filter bandwidth.
  dfbw = (1 + expand) * Fw;
  sw = parms.sw.hz;

  % calculate the decimation ratio.
  D = round(sw / dfbw);

  % build the modulation values.
  Amod = reshape(exp(-i .* 2 .* pi .* F0 .* t), length(t), 1);

  % build the filter coefficients.
  ht = 0.5 .* t(N) .* linspace(-1, 1, N)';
  h = sinc(dfbw .* ht) .* blackman(N);

  % apply the extraction operations.
  if (isvector(f))
    % modulate the single decay.
    fmod = f .* Amod;

    % filter the modulated decay.
    ffir = shift(filter(h, 1, fmod), -N / 2);

    % decimate the filtered, modulated decay.
    fsub = ffir([1 : D : end]);
    tsub = t([1 : D : end]);
  elseif (ismatrix(f))
    % modulate the set of decays.
    fmod = f .* (ones(rows(f), 1) * Amod');

    % filter the modulated decays.
    ffir = zeros(size(fmod));
    for i = 1 : rows(ffir)
      % filter the currently indexed decay.
      ffir(i,:) = shift(filter(h, 1, fmod(i,:)')', -N / 2);
    end

    % decimate the filtered, modulated decays.
    fsub = ffir(:, [1 : D : end]);
    tsub = t([1 : D : end]);
  else
    % throw an exception.
    error('subfid1d: fid data must be a vector or a matrix');
  end
end

