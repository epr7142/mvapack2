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
## @anchor{decompose1d}
## @deftypefn {Function File} {@var{T} =} decompose1d (@var{f}, @var{t}, @var{parms})
## @deftypefnx {Function File} {@var{T} =} decompose1d (@var{f}, @var{t}, @var{parms}, @var{roi})
## @deftypefnx {Function File} {@var{T} =} decompose1d (@var{f}, @var{t}, @var{parms}, @var{roi}, @var{minbw})
## @deftypefnx {Function File} {@var{T} =} decompose1d (@var{f}, @var{t}, @var{parms}, @var{roi}, @var{minbw}, @var{fitopts})
## Performs Complete Reduction to Amplitude and Frequency Table (CRAFT)
## analysis of a 1D time-domain NMR free induction decay (FID) in @var{f},
## with a time abscissa in @var{t} and parameters in @var{parms}.
##
## It is highly recommended that you use @ref{decompose} instead of calling
## this function directly.
## @end deftypefn

function T = decompose1d (f, t, parms, roi, minbw, fitopts)
  % check if the number of expected arguments was passed.
  if (nargin != 6 || nargout != 1)
    % print the usage statement.
    print_usage();
  end

  % enforce one-dimensional parameters.
  if (iscell(parms))
    parms = parms{1};
  end

  % check if an ROI argument was provided.
  if (nargin >= 4 && !isempty(roi))
    % check the type of the ROI argument.
    if (!(ismatrix(roi) && columns(roi) == 2) && !is_function_handle(roi))
      % invalid argument type. throw an exception.
      error('decompose1d: roi argument must be a function handle or a matrix');
    end
  else
    % use the default automatic ROI generation routine.
    roi = @roibin;
  end

  % check if a minimum bandwidth was provided.
  if (nargin >= 5 && !isempty(minbw))
    % check the type of the argument.
    if (!isscalar(minbw) && !isreal(minbw))
      % invalid argument type. throw an exception.
      error('decompose1d: minbw argument must be a real scalar');
    end
  else
    % use the default minimum bandwidth.
    minbw = parms.sw.hz / 100;
  end

  % check if an options structure was passed.
  if (nargin >= 6 && !isempty(fitopts))
    % check the type.
    if (!isstruct(fitopts))
      % invalid type.
      error('decompose1d: fit options argument must be a structure');
    end
  else
    % pass an empty options structure to the fitting routine.
    fitopts = [];
  end

  % calculate the data matrix sizes.
  if (isvector(f))
    % single observation.
    N = 1;
    K = length(f);
    T = [];
  elseif (ismatrix(f))
    % multiple observations.
    [N, K] = size(f);
    T = cell(N, 1);
  else
    % invalid data type. throw an exception.
    error('decompose1d: input data must be a vector or matrix');
  end

  % run quick and dirty automated processing of the time-domain data.
  [s, ppm, hz] = nmrft(f, parms);
  s = real(autophase(s, parms));

  % see if we need to automatically build regions of interest.
  if (is_function_handle(roi))
    % yes. build them now.
    R = roi(s, hz, parms, minbw);
    roi = R;
  end

  % get the number of ROIs.
  nroi = rows(roi);

  % initialize the extracted sub-FID information.
  fsub = cell(nroi, 1);
  tsub = cell(nroi, 1);
  dfbw = cell(nroi, 1);
  D = cell(nroi, 1);

  % loop through all the regions of interest.
  for idx = 1 : nroi
    % get the frequency boundaries of the region.
    Fmin = roi(idx, 1);
    Fmax = roi(idx, 2);

    % extract the frequency band of interest from the total FID.
    [fs, ts, bs, ds] = subfid(f, t, parms, Fmin, Fmax);

    % store the extracted sub-FID information.
    fsub{idx} = fs;
    tsub{idx} = ts;
    dfbw{idx} = bs;
    D{idx} = ds;
  end

  % loop again through the regions of interest.
  for idx = 1 : nroi
    % pull the ROI-based values into more easily digestable variables.
    fsubi = fsub{idx};
    tsubi = tsub{idx};
    dfbwi = dfbw{idx};
    Di = D{idx};

    % check how many observations we must process.
    if (N == 1)
      % estimate peak parameters from the extracted sub-FID.
      Ti = rjmcmc1d(fsubi, tsubi, fitopts);

      % shift the frequencies back into baseband.
      omega = Ti(:,2);
      omega = ifelse(omega > dfbwi / 2, omega - dfbwi, omega);
      omega += mean(roi(idx,:));
      Ti(:,2) = omega;

      % append the signals to the output table.
      T = [T; Ti];
    else
      % loop for each extracted sub-FID.
      for n = 1 : N
        % estimate peak parameters from the extracted sub-FID.
        Ti = rjmcmc1d(fsubi(n,:), tsubi, fitopts);

        % shift the frequencies back into baseband.
        omega = Ti(:,2);
        omega = ifelse(omega > dfbwi / 2, dfbwi - omega, omega);
        omega += mean(roi(idx,:));
        Ti(:,2) = omega;

        % append the signals to the output table.
        T{n} = [T{n}; Ti];
      end
    end
  end
end

