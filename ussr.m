## Copyright (C) 2014 University of Nebraska Board of Regents.
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
## @deftypefn {Function File} {@var{R} =} ussr (@var{F}, @var{P}, @var{ppm})
## @deftypefnx {Function File} {@var{R} =} ussr (@var{F}, @var{P}, @var{ppm}, @var{alpha})
## Perform Uncomplicated Statistical Spectral Remodeling (USSR) on
## two matrices @var{F} and @var{P}, where spectra are paired in the
## two matrices, arranged row-wise. The @var{P} matrix is remodeled
## to yield @var{R} according to the procedure in:
##
## B. Worley et. al., 'Uncomplicated Statistical Spectral Remodeling',
## J. Biomol. NMR., in preparation.
##
## The chemical shift values on the abscissas of all spectra are assumed
## to be identical to within the digital resolution of the experiment, so
## only a single abscissa vector @var{ppm}.
##
## An optional fourth input variable @var{alpha} may be set to define the level
## of confidence required to keep a signal in the reconstructed spectrum.
## @end deftypefn

function [R, mu, sigma, scale, T] = ussr (F, P, ppm, alpha)
  % get the size of the datasets.
  N = rows(F);
  K = columns(F);

  % ensure the dataset sizes all match each other.
  if (N != rows(P) || K != columns(P) || K != length(ppm))
    % one of the sizes has a mismatch. output an error.
    error('size mismatch identified in input variables');
  end

  % see if a critical alpha value was passed.
  if (nargin < 4 || isempty(alpha))
    % no value was passed. assume 0.05.
    alpha = 0.05;
  end

  % build a gaussian kernel for smoothing stuff and things.
  gkern = fspecial('gaussian', [1 7]);

  % initialize the peak indicator matrix for the free spectra.
  Fpeaks = zeros(N, K);

  % compute the noise mean and standard deviation of each free spectrum.
  [Rmean, Rstd] = estnoise(F, cell(1, 1));

  % loop through the spectra to calculate the peak indicators.
  for i = 1 : N
    % calculate a threshold in order to identify true peaks.
    Fpeaks(i,:) = 1 .* (F(i,:) > Rmean(i) + 5 * Rstd(i));

    % smooth based on the gaussian kernel.
    Fpeaks(i,:) = imfilter(Fpeaks(i,:), gkern);
  end

  % calculate the differences matrix.
  D = P - F;

  % estimate the true baseline mean from the median.
  mu = median(D);

  % estimate the true baseline deviation from median absolute deviation.
  sigma = median(abs(D - ones(N, 1) * mu)) ./ norminv(0.75, 0, 1);

  % initialize the matrix of t-values.
  T = zeros(N, K);

  % initialize the matrix of remodeled spectra.
  R = zeros(N, K);

  % initialize the mean scale factors.
  scale = ones(N, 1);

  % loop through the spectra.
  for i = 1 : N
%    % loop a while to better fit the baseline into the spectrum.
%    for j = 1 : 20
%      % calculate the t- and p-values for the currently indexed spectrum.
      t = (P(i,:) - scale(i) .* mu) ./ (sigma ./ sqrt(N));
      p = 1.0 - tcdf(t, N - 1);
%
%      % calculate the scale factor from the p-values.
%      scale(i) = sum(p .* P(i,:) .* mu) / sum(p .* mu .* mu);
%    end

    % calculate the weights from the p-values and spectral intensities.
    w = ones(1, K) .* (p < alpha / K) .* (P(i,:) > mu);

    % loop a few times to smooth out the weights.
    for j = 1 : 7
      % smooth based on the gaussian kernel.
      w = imfilter(w, gkern);
    end

    % estimate the noise level of the spectrum using far-upfield intensity.
    estnoi = normrnd(0, 0.2 .* std(F(i,find(ppm < -0.1))), 1, K);

    % calculate the remodeled spectrum.
    w = w .* Fpeaks(i,:);
    smu = scale(i) .* mu;
    R(i,:) = w .* max(P(i,:) - smu, estnoi) + (1 - w) .* estnoi;

    % store the t-value row.
    T(i,:) = t;
  end
end

