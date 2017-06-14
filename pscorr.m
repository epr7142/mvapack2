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
## @anchor{pscorr}
## @deftypefn {Function File} {@var{Xn} =} pscorr (@var{X})
## @deftypefnx {Function File} {[@var{Xn}, @var{b}] =} pscorr (@var{X})
## @deftypefnx {Function File} {@var{Xn} =} pscorr (@var{X}, @var{r})
## @deftypefnx {Function File} {[@var{Xn}, @var{b}] =} pscorr (@var{X}, @var{r})
## Use simultaneous phase-scatter correction to normalize spectra into better
## alignment with each other, row-wise at least. If the optional @var{r}
## reference observation is not provided, the mean of the observations in
## @var{X} will be used. Normalization factors may be requested through the
## second optional return value @var{b}. More information here:
##
## @quotation
## B. Worley, R. Powers, `Simultaneous Phase and Scatter Correction in NMR
## Datasets', Chemometr. Intell. Lab. Syst., 2013, submitted.
## @end quotation
## @end deftypefn

function [Xnew, b] = pscorr (X, r)
  % check for proper arguments.
  if (!any(nargin == [1 : 2]) || !any(nargout == [1 : 2]) || !ismatrix(X))
    % improper arguments. print the usage statement.
    print_usage();
  end

  % ensure the input data matrix is complex.
  if (!iscomplex(X))
    % throw an exception.
    error('pscorr: complex data is required for processing');
  end

  % check if the data is actually a vector.
  if (isvector(X))
    % reshape it into a row vector.
    X = reshape(X, 1, length(X));
  end

  % check if the reference spectrum was provided.
  if (nargin >= 2 && !isempty(r) && isvector(r))
    % yes. ensure it's properly shaped.
    r = reshape(r, 1, length(r));
  else
    % no. use the default value, the mean.
    r = mean(X);
  end

  % calculate the centered spectrum.
  r = real(r);
  rbar = mean(r);
  rc = r - rbar;

  % initialize the output matrix.
  Xnew = X;
  Xc = center(X, 2);

  % check if the second output argument was requested.
  if (nargout >= 2)
    % yes. initialize it.
    b = ones(rows(X), 1);
  end

  % loop through the observations of the data matrix.
  for i = 1 : rows(X)
    % extract the current observation.
    x = [real(Xc(i,:)'), imag(Xc(i,:)')];

    % set up the arguments to leasqr().
    pin = [0; 0; 1];
    stol = 1e-4;
    niter = 100;
    wt = ones(size(rc'));
    dp = 1e-3 .* ones(size(pin));

    % calculate the optimal correction factors.
    [f, pout] = leasqr(x, rc', pin, @pscorr_func, ...
                       stol, niter, wt, dp, 'dfdp', {});

    % store the corrected spectrum using the optimal factors.
    Xnew(i,:) = (pscorr_build(x, pout) + rbar)';

    % check if the normalization factors were requested.
    if (nargout >= 2)
      % yes. store the factor.
      b(i) = 1 / pout(3);
    end
  end

  % check if the data is actually a vector.
  if (isvector(Xnew))
    % reshape it back into a column vector.
    Xnew = reshape(Xnew, length(Xnew), 1);
  end
end

