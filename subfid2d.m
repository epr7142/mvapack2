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
## @anchor{subfid2d}
## @deftypefn {Function File} {[@var{fsub}, @var{tsub}, @var{dfbw}, @var{D}] =} subfid2d (@var{f}, @var{t}, @var{parms}, @var{Fmin}, @var{Fmax})
## Extracts a band of frequencies (@math{[F_{min},F_{max}]}) from @var{f} into
## @var{fsub}.
##
## It is highly recommended that you use @ref{subfid} instead of calling this
## function directly.
## @end deftypefn

function [fsub, tsub, dfbw, D] = subfid2d (f, t, parms, Fmin, Fmax)
  % check if the number of expected arguments was passed.
  if (nargin != 5 || nargout != 4)
    % print the usage statement.
    print_usage();
  end

  % check the decay data argument.
  if (!(ismatrix(f) && iscomplex(f)) && !iscell(f))
    % invalid type. throw an exception.
    error('subfid2d: decay data must be a complex matrix or a cell array');
  end

  % check the time axis argument.
  if (!iscell(t) || length(t) != 2 || ...
      !isvector(t{1}) || !isreal(t{1}) || ...
      !isvector(t{2}) || !isreal(t{2}))
    % invalid type. throw an exception.
    error('subfid2d: time axis must be a cell array of two real vectors');
  end

  % check the parameter argument.
  if (!iscell(parms) || length(parms) < 2)
    % invalid type. throw an exception.
    error('subfid2d: parameters must be a two-element cell array');
  end

  % check the frequency boundary arguments.
  if (!isvector(Fmin) || !isreal(Fmin) || length(Fmin) != 2 || ...
      !isvector(Fmax) || !isreal(Fmax) || length(Fmax) != 2)
    % invalid type. throw an exception.
    error('subfid2d: frequency bounds must be real two-vectors');
  end

  % reshape the frequency boundaries to column vectors.
  Fmin = reshape(Fmin, 2, 1);
  Fmax = reshape(Fmax, 2, 1);

  % initialize some default options.
  % @expand: fractional amount to expand the FIR filter bandwidth.
  % @N: number of taps to use for the FIR filter.
  expand = 0.1;
  N = [512; 64];

  % calculate the window center frequency and bandwidth.
  F0 = (Fmin + Fmax) ./ 2;
  Fw = abs(Fmax - Fmin);

  % calculate the digital filter bandwidth.
  dfbw = (1 + expand) .* Fw;

  % grab the spectral widths.
  sw = cellfun(@(x) x.sw.hz, parms);
  sw = sw(1:2);

  % calculate the decimation ratio. this step ensures that a minimum
  % number of complex points in each dimension (16) is always kept, so
  % decimation will never drop 'all but a single point'.
  D = min(round(sw ./ dfbw), size(f')' ./ 16);

  % initialize the modulation and filtering coefficients.
  Amod = cell(2, 1);
  h = cell(2, 1);

  % loop for each dimension.
  for idx = 1 : 2
    % get the current time vector.
    ti = t{idx};

    % build the modulation values.
    Amod{idx} = reshape(exp(-i .* 2 .* pi .* F0(idx) .* ti), length(ti), 1);

    % build the filter coefficients.
    ht = 0.5 .* ti(N(idx)) .* linspace(-1, 1, N(idx))';
    h{idx} = sinc(dfbw(idx) .* ht) .* blackman(N(idx));
  end

  % apply the extraction operations.
  if (ismatrix(f))
    % initialize the time vectors.
    tsub = cell(2, 1);
    t1 = t{1};
    t2 = t{2};

    % modulate along the first dimension.
    fmod = f .* (ones(rows(f), 1) * Amod{1}');

    % filter the modulated decays along the first dimension.
    ffir = zeros(size(fmod));
    for idx = 1 : rows(ffir)
      % filter the currently indexed decay.
      ffir(idx,:) = shift(filter(h{1}, 1, fmod(idx,:)')', -N(1) / 2);
    end

    % decimate the filtered, modulated decays.
    fsub = ffir(:, [1 : D(1) : end]);
    tsub{1} = t1([1 : D(1) : end]);

    % extract the two interleaved submatrices for second dimension filtering.
    [P, Q] = states(fsub);

    % modulate along the second dimension.
    Pmod = P .* (Amod{2} * ones(1, columns(P)));
    Qmod = Q .* (Amod{2} * ones(1, columns(Q)));

    % filter the modulated decays along the second dimension.
    Pfir = zeros(size(Pmod));
    Qfir = zeros(size(Qmod));
    for idx = 1 : columns(Pfir)
      % filter the currently indexed decay.
      Pfir(:,idx) = shift(filter(h{2}, 1, Pmod(:,idx)), -N(2) / 2);
      Qfir(:,idx) = shift(filter(h{2}, 1, Qmod(:,idx)), -N(2) / 2);
    end

    % decimate the filtered, modulated decays.
    Psub = Pfir([1 : D(2) : end], :);
    Qsub = Qfir([1 : D(2) : end], :);
    tsub{2} = t2([1 : D(2) : end]);

    % re-interleave the filtered submatrices.
    fsub = states(Psub, Qsub);
  elseif (iscell(f))
    % initialize the output cell array.
    fsub = cell(length(f), 1);

    % filter each matrix in the cell array.
    for idx = 1 : length(f)
      % filter the matrix.
      fsub{idx} = subfid2d(f{idx}, t, parms, Fmin, Fmax);
    end
  else
    % throw an exception.
    error('subfid2d: fid data must be a matrix or a cell array');
  end
end

