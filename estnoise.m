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
## @anchor{estnoise}
## @deftypefn {Function File} {@var{nfloor} =} estnoise (@var{s}, @var{parms})
## @deftypefnx {Function File} {[@var{mu}, @var{sigma}] =} estnoise (@var{s}, @var{parms})
## Roughly estimates the mean and variance of the spectral baseline in a
## one-dimensional spectrum vector or a two-dimensional spectrum matrix.
## If only a single return value is requested, then the noise floor,
## defined as the sum of the mean and two times the standard deviation,
## will be reported.
## @end deftypefn

function [out1, out2] = estnoise (s, parms)
  % check if the number of expected arguments was passed.
  if (nargin != 2 || !any(nargout == [1 : 2]))
    % print the usage statement.
    print_usage();
  end

  % check if the spectrum contains complex values.
  if (iscomplex(s))
    % use only the real portion of the complex spectrum.
    s = real(s);
  end

  % get the number of dimensions in the data.
  nd = nmrdims(s, parms);

  % estimate noise based on dimensionality.
  if (nd == 1)
    % check the type of the data argument.
    if (isvector(s))
      % get the noise parameters of the spectrum.
      [mu, sigma] = __estnoise1d(s);
    elseif (ismatrix(s))
      % initialize the parameter vectors.
      mu = zeros(rows(s), 1);
      sigma = zeros(rows(s), 1);

      % loop through the spectra in the dataset.
      for idx = 1 : rows(s)
        % get the noise parameters of the current spectrum.
        [mu(idx), sigma(idx)] = __estnoise1d(s(idx,:)');
      end
    else
      % invalid type.
      error('estnoise: one-dimensional data must be a vector or a matrix');
    end
  elseif (nd == 2)
    % check the type of the data argument.
    if (ismatrix(s))
      % get the noise parameters of the spectrum.
      [mu, sigma] = __estnoise2d(s);
    elseif (iscell(s))
      % initialize the parameter vectors.
      mu = zeros(length(s), 1);
      sigma = zeros(length(s), 1);

      % loop through the spectra in the dataset.
      for idx = 1 : length(s)
        % get the noise parameters of the current spectrum.
        [mu(idx), sigma(idx)] = __estnoise2d(s{idx});
      end
    else
      % invalid type.
      error('estnoise: two-dimensional data must be a matrix or an array');
    end
  else
    % invalid dimensionality. throw an exception.
    error('estnoise: input data must be one- or two-dimensional');
  end

  % check how many return values were requested.
  if (nargout == 1)
    % return only the estimated noise floor.
    out1 = mu + 2 .* sigma;
  elseif (nargout == 2)
    % return the mean and standard deviation.
    out1 = mu;
    out2 = sigma;
  end
end

