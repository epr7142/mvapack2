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
## @anchor{rjmcmc2d}
## @deftypefn {Function File} {@var{P} =} rjmcmc2d (@var{y}, @var{t})
## @deftypefnx {Function File} {@var{P} =} rjmcmc2d (@var{y}, @var{t}, @var{opts})
## Execute a Reversible Jump Markov Chain Monte Carlo (RJ-MCMC) simulation
## in order to obtain estimates of the number of signals in a 2D time-domain
## free induction decay @var{y}, and the corresponding parameters of those
## signals. A two-element time cell array @var{t} is required to be paired
## with @var{y}, the input time-domain matrix.
##
## A third (optional) argument may be passed to specify options for the fitting
## algorithm. The argument, @var{opts}, must be a structure and is expected to
## contain any of the following option fields:
##
## The option field @var{iters} may be passed as a two-element vector,
## whose first element is the number of desired burn-in iterations, and whose
## second element is the number of desired simulation iterations. @var{iters}
## may also contain a third element that indicates the modulus for which
## iteration to save during the simulation (@emph{i.e.} only save every
## hundredth iteration).
##
## The option field @var{lambda} may be passed to indicate the expected
## number of signals the model will contain. The default value is 5. The
## option field @var{delta} may be passed to indicate the expected signal
## to noise ratio of the data. The default value is 10. Further option fields
## @var{omega0} and @var{rho0} may be provided to initialize the simulation at
## more sane parameter values (frequency and decay rate, respectively).
##
## The RJ-MCMC algorithm is an extension of @ref{rjmcmc1d} to two-dimensional
## time-domain data matrices. See @ref{rjmcmc1d} for more information and
## literature references.
## @end deftypefn

function P = rjmcmc2d (y, t, opts)
  % check the arguments.
  if (nargin != 3 || nargout != 1)
    % invalid arguments.
    print_usage();
  end

  % check the type of the time array.
  if (!iscell(t) || length(t) != 2 || ...
      !isvector(t{1}) || !isreal(t{1}) || ...
      !isvector(t{2}) || !isreal(t{2}))
    % invalid type. throw an exception.
    error('rjmcmc2d: time axis must be an array of two real vectors');
  end

  % check the type of the time data.
  if (!ismatrix(y) || !iscomplex(y))
    % invalid type. throw an exception.
    error('rjmcmc2d: time-domain data must be a complex matrix');
  end

  % get the lengths of the time vectors.
  tlen = cellfun(@(x) length(x), t);

  % reshape the time axis vectors into column format.
  for idx = 1 : 2
    t{idx} = reshape(t{idx}, tlen(idx), 1);
  end

  % ensure a parameter array exists in the options structure.
  if (!isfield(opts, 'parms'))
    % throw an exception.
    error('rjmcmc2d: parameters are required');
  end

  % grab the parameter array.
  parms = opts.parms;

  % check if an options structure was passed.
  if (nargin >= 3 && !isempty(opts) && isstruct(opts))
    % check for iteration parameters.
    if (isfield(opts, 'iters'))
      % use the option field.
      iters = opts.iters;
    end

    % check for an expected signal count.
    if (isfield(opts, 'lambda'))
      % use the option field.
      lambda = opts.lambda;
    end

    % check for an expected signal to noise ratio.
    if (isfield(opts, 'delta'))
      % use the option field.
      delta = opts.delta;
    end

    % check for an explicit frequency axis.
    if (isfield(opts, 'w'))
      % use the option field.
      w = opts.w;
    end

    % check for initial frequency values.
    if (isfield(opts, 'omega0'))
      % use the option field.
      omega0 = opts.omega0;
    end

    % check for initial decay rate values.
    if (isfield(opts, 'rho0'))
      % use the option field.
      rho0 = opts.rho0;
    end

    % check for threshold values.
    if (isfield(opts, 'thresh'))
      % use the option field.
      thresh = opts.thresh;
    end

    % check for a verbosity flag.
    if (isfield(opts, 'verbose'))
      % use the option field.
      verbose = opts.verbose;
    end
  end

  % check if iteration parameters were passed.
  if (exist('iters'))
    % check the type of the iteration argument.
    if (isvector(iters) && any(length(iters) == [2, 3]))
      % use the specified iteration counts.
      n_burn = iters(1);
      n_sim = iters(2);

      % see if a third count was passed.
      if (length(iters) == 3)
        % yes. use the specified modulus.
        n_mod = iters(3);
      else
        % no. save every iteration.
        n_mod = 1;
      end
    else
      % invalid type. throw an exception.
      error('rjmcmc2d: iteration argument must be a 2/3-element vector');
    end
  else
    % use the default iteration counts.
    n_burn = 1000;
    n_sim = 4000;
    n_mod = 1;
  end

  % check that the modulus evenly divides the simulation iteration count.
  if (mod(n_sim, n_mod) != 0)
    % this can't be allowed to happen. throw an exception.
    error('rjmcmc2d: iteration storage modulus must evenly divide count');
  end

  % check if an expected signal count was passed.
  if (exist('lambda'))
    % check the argument type.
    if (!isscalar(lambda) || !isreal(lambda) || lambda < 1)
      % invalid argument. throw an exception.
      error('rjmcmc2d: lambda must be a real positive scalar');
    end
  else
    % use the default value.
    lambda = 5;
  end

  % check if an expected signal to noise ratio was passed.
  if (exist('delta'))
    % check the argument type.
    if (!isscalar(delta) || !isreal(delta) || delta <= 0)
      % invalid argument. throw an exception.
      error('rjmcmc2d: delta must be a real positive scalar');
    end
  else
    % use the default value.
    delta = 10;
  end

  % see if a frequency axis was passed.
  if (exist('w'))
    % check the type of the frequency axis.
    if (!iscell(w) || length(w) != 2 || ...
        !isvector(w{1}) || !isreal(w{1}) || ...
        !isvector(w{2}) || !isreal(w{2}))
      % invalid type. throw an exception.
      error('rjmcmc2d: frequency axis must be an array of two real vectors');
    end
  else
    % build uniform frequency axes (as a best guess).
    w = cell(2, 1);
    for idx = 1 : 2
      w_max = 1 / min(diff(t{idx}));
      w{idx} = linspace(0, w_max, 2 ^ nextpow2(tlen(idx)))';
    end
  end

  % check if an initial decay vector was passed.
  if (exist('rho0'))
    % check the type of the initial decay rates matrix.
    if (!ismatrix(rho0) || !isreal(rho0) || columns(rho0) != 2)
      % invalid type. throw an exception.
      error('rjmcmc2d: decay rates matrix must be a real two-column matrix');
    end

    % check that the number of decays matches the expected signal count.
    if (rows(rho0) != rows(omega0))
      % throw an exception.
      error('rjmcmc2d: decay rate and frequency vector lengths must match.');
    end
  else
    % initialize the decay rate parameters.
    rho0 = 1.0 .* ones(lambda, 2);
  end

  % check if threshold values were passed.
  if (exist('thresh'))
    % check the type of the threshold argument.
    if (!isvector(thresh) || !isreal(thresh) || length(thresh) != 2)
      % throw an exception.
      error('rjmcmc2d: thresholds must be a two-element real vector.');
    end
  else
    % initialize the default thresholds.
    thresh = [2, pi];
  end

  % check if a verbosity flag was passed.
  if (exist('verbose'))
    % check the type of the verbosity argument.
    if (!isbool(verbose))
      % throw an exception.
      error('rjmcmc1d: verbosity flag must be a boolean');
    end
  else
    % initialize the default verbosity.
    verbose = false;
  end

  % initialize the amplitude parameters.
  a0 = complex(zeros(lambda, 1), zeros(lambda, 1));

  % apodize the FID in the direct dimension, just to remove any signal that
  % appears to have a negative decay rate. such signals induce massive
  % artifacts during transformation into frequency space.
  y = (ones(rows(y), 1) * sinewindow(t{1})') .* y;

  % calculate the frequency transform of the signal.
  Y = nmrft(y, parms, false);

  % perform an initial automatic phase correction of the signal.
  [Y, phi0, phi1] = autophase2d(Y, @simplex_entropy);

  % whitening makes no distinction between positive and negative peaks, but
  % the RJ-MCMC algorithm needs all peaks to be positive. make it so.
  if (max(vec(real(states(Y)))) < -min(vec(real(states(Y)))))
    % apply a 180-degree phase correction.
    phi0(1) += 180;
  end

  % check if an initial frequency vector was passed.
  if (exist('omega0'))
    % check the type of the initial frequencies matrix.
    if (!ismatrix(omega0) || !isreal(omega0) || columns(omega0) != 2)
      % invalid type. throw an exception.
      error('rjmcmc2d: frequencies matrix must be a real two-column matrix');
    end

    % check that the number of frequencies matches the expected signal count.
    if (rows(omega0) != lambda)
      % update the lambda parameter.
      lambda = rows(omega0);
    end
  else
    % randomly initialize the frequency parameters.
    omega0 = __rjmcmc2d_freqrnd(w, states(Y), lambda);
  end

  % use the phase correction gained from the spectrum to phase the signal
  % in the time domain. this ensures mostly absorptive peaks in the real
  % spectrum, which allows the sampler to discard phase angles during
  % its random walk (and avoid adding dispersive signals).
  y = phase2d(y, phi0, [0; 0]);

  % extract and vectorize the cosine- and sine-modulated slices from the
  % time-domain data.
  ycos = y(2 : 2 : end, :);
  ysin = y(1 : 2 : end - 1, :);

  % run the compiled RJ-MCMC subroutine.
  [omega1, omega2, rho1, rho2, a, K, ...
   omega1_map, omega2_map, rho1_map, rho2_map, a_map] = ...
    __rjmcmc2d(t, ycos, ysin, w, states(Y), ...
               n_burn, n_sim, n_mod, lambda, delta, ...
               omega0, rho0, a0, verbose);

  % sort the MAP parameters by first-dimension frequency.
  [omega1_map, idx] = sort(omega1_map);
  omega2_map = omega2_map(idx);
  rho1_map = rho1_map(idx);
  rho2_map = rho2_map(idx);
  a_map = a_map(idx);

  % finalize the output peak table.
  P = [a_map, omega1_map, omega2_map, rho1_map, rho2_map, ...
       zeros(rows(a_map), 2)];
end

