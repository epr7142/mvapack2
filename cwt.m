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
## @anchor{cwt}
## @deftypefn {Function File} {@var{C} =} cwt (@var{x})
## @deftypefnx {Function File} {@var{C} =} cwt (@var{x}, @var{w})
## @deftypefnx {Function File} {@var{C} =} cwt (@var{x}, @var{w}, @var{voices})
## @deftypefnx {Function File} {@var{C} =} cwt (@var{x}, @var{w}, @var{voices}, @var{octave})
## @deftypefnx {Function File} {@var{C} =} cwt (@var{x}, @var{w}, @var{voices}, @var{octave}, @var{scale})
## @deftypefnx {Function File} {@var{C} =} cwt (@var{x}, @var{w}, @var{voices}, @var{octave}, @var{scale}, @var{doifft})
## Performs continuous wavelet transformation of a vector signal to produce a
## time-frequency transform (CWT) matrix, @var{C}. The default wavelet is
## the Mexican hat (@xref{wavelet_sombrero}), but can be set to another
## function by passing @var{w} with a function handle of the form:
##
## @code{function psi = wavelet_function (t)
##   ...
## end}
##
## The optional argument @var{voices} sets the number of voices per octave,
## and defaults to 8. The optional argument @var{octave} sets the initial
## octave, and defaults to 4. The optional argument @var{scale} sets the
## initial scale, and defaults to 8. You can tweak these numbers to your
## satisfaction.
##
## The final optional argument @var{doifft} sets whether the function should
## perform a final inverse fourier transform. This is useful if you wish to
## calculate derivatives in frequency space. By default, @var{doifft} is set
## to @code{true}.
## @end deftypefn

function C = cwt (x, w, voices, octave, scale, doifft)
  % check if the number of expected arguments was passed.
  if (!any(nargin == [1 : 6]) || nargout != 1 || !isvector(x))
    % print the usage statement.
    print_usage();
  end

  % check for the wavelet argument.
  if (nargin >= 2 && !isempty(w))
    % check the type of the argument.
    if (!is_function_handle(w))
      % invalid type. throw an exception.
      error('cwt: wavelet argument must be a valid function handle');
    end
  else
    % use the default wavelet: sombrero.
    w = @wavelet_sombrero;
  end

  % check for the voices argument.
  if (nargin >= 3 && !isempty(voices))
    % check the type of the argument.
    if (!isscalar(voices) || !isreal(voices))
      % invalid type. throw an exception.
      error('cwt: voices argument must be a real scalar');
    end
  else
    % use the default value: 8.
    voices = 8;
  end

  % check for the octave argument.
  if (nargin >= 4 && !isempty(octave))
    % check the type of the argument.
    if (!isscalar(octave) || !isreal(octave))
      % invalid type. throw an exception.
      error('cwt: octave argument must be a real scalar');
    end
  else
    % use the default value: 4.
    octave = 4;
  end

  % check for the scale argument.
  if (nargin >= 5 && !isempty(scale))
    % check the type of the argument.
    if (!isscalar(scale) || !isreal(scale))
      % invalid type. throw an exception.
      error('cwt: scale argument must be a real scalar');
    end
  else
    % use the default value: 8.
    scale = 8;
  end

  % check for the doifft argument.
  if (nargin >= 6 && !isempty(doifft))
    % check the type of the argument.
    if (!isbool(doifft))
      % invalid type. throw an exception.
      error('cwt: doifft argument must be a boolean');
    end
  else
    % use the default value: true.
    doifft = true;
  end

  % get the length of the signal, to the closest power of two.
  n = 2 ^ nextpow2(length(x));

  % calculate the number of octaves and scales.
  octaves = log2(n) - octave;
  scales = voices * octaves;

  % reshape, pad and fourier transform the input signal.
  x = reshape(postpad(x, n), n, 1);
  xhat = fft(x);

  % initialize the column counter and initial scale.
  sc = scale;
  j = 1;

  % initialize the output matrix.
  C = zeros(n, scales);

  % loop through the octaves.
  for oc = 1 : octaves
    % and loop through the voices.
    for vc = 1 : voices
      % calculate the scale of the current voice.
      voice = sc * 2 ^ (vc / voices);

      % calculate the scaled, fourier transformed daughter wavelet.
      t = (2 * pi / voice) .* linspace(-n / 2, n / 2, n)';
      what = abs(fft(w(t) ./ sqrt(voice)));

      % multiply the frequency-space wavelet and signal.
      chat = xhat .* what;

      % store the current wavelet transform.
      C(:,j) = chat;

      % increment the column counter.
      j++;
    end

    % go to the next scale.
    sc *= 2;
  end

  % should we perform the inverse fourier transform?
  if (doifft == true)
    % yes, perform the transformation.
    C = ifft(C);
  end
end

