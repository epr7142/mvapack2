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
## @anchor{rjmcmc1d}
## @deftypefn {Function File} {@var{P} =} rjmcmc1d (@var{y}, @var{t})
## @deftypefnx {Function File} {@var{P} =} rjmcmc1d (@var{y}, @var{t}, @var{opts})
## Execute a Reversible Jump Markov Chain Monte Carlo (RJ-MCMC) simulation
## in order to obtain estimates of the number of signals in a time-domain
## free induction decay @var{y}, and the corresponding parameters of those
## signals. A time axis @var{t} is required to be paired with @var{y}.
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
## number of signals the model will contain. The default value is 20. The
## option field @var{delta} may be passed to indicate the expected signal
## to noise ratio of the data. The default value is 10. Further option fields
## @var{omega0} and @var{rho0} may be provided to initialize the simulation at
## more sane parameter values (frequency and decay rate, respectively).
##
## The RJ-MCMC algorithm is implemented according to information reported in:
##
## @quotation
## Andrieu C., `Joint Bayesian Model Selection and Estimation of Noisy
## Sinusoids via Reversible Jump MCMC', IEEE Transactions on Signal
## Processing, 1999.
## @end quotation
##
## @quotation
## Rubtsov D., Griffin J., `Time-domain Bayesian detection and estimation of
## noisy damped sinusoidal signals applied to NMR spectroscopy', Journal of
## Magnetic Resonance, 2007.
## @end quotation
##
## @quotation
## Roodaki A., `Note on the computation of the Metropolis-Hastings ratio for
## Birth-or-Death moves in trans-dimensional MCMC algorithms for signal
## decomposition problems', arXiv:1111.6245v2, 2012.
## @end quotation
## @end deftypefn

function P = rjmcmc1d (y, t, opts)
  % check the arguments.
  if (!any(nargin == [2 : 3]) || nargout != 1)
    % invalid arguments.
    print_usage();
  end

  % check the type of the time axis.
  if (!isvector(t) || !isreal(t))
    % invalid type. throw an exception.
    error('rjmcmc1d: time axis must be a real vector');
  end

  % check the type of the time data.
  if (!isvector(y) || !iscomplex(y))
    % invalid type. throw an exception.
    error('rjmcmc1d: time data must be a complex vector');
  end

  % reshape the time axis vector.
  t = reshape(t, length(t), 1);

  % reshape the time data vector.
  y = reshape(y, length(y), 1);

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
      error('rjmcmc1d: iteration argument must be a 2/3-element vector');
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
    error('rjmcmc1d: iteration storage modulus must evenly divide count');
  end

  % check if an expected signal count was passed.
  if (exist('lambda'))
    % check the argument type.
    if (!isscalar(lambda) || !isreal(lambda) || lambda < 1)
      % invalid argument. throw an exception.
      error('rjmcmc1d: lambda must be a real positive scalar');
    end
  else
    % use the default value.
    lambda = 20;
  end

  % check if an expected signal to noise ratio was passed.
  if (exist('delta'))
    % check the argument type.
    if (!isscalar(delta) || !isreal(delta) || delta <= 0)
      % invalid argument. throw an exception.
      error('rjmcmc1d: delta must be a real positive scalar');
    end
  else
    % use the default value.
    delta = 10;
  end

  % see if a frequency axis was passed.
  if (exist('w'))
    % check the type of the frequency axis.
    if (!isvector(w) || !isreal(w))
      % invalid type. throw an exception.
      error('rjmcmc1d: frequency axis must be a real vector');
    end
  else
    % build a uniform frequency axis (as a best guess).
    w = linspace(0, 1 / min(diff(t)), 2 ^ nextpow2(length(t)))';
  end

  % check if an initial frequency vector was passed.
  if (exist('omega0'))
    % check that the number of frequencies matches the expected signal count.
    if (length(omega0) != lambda)
      % update the lambda parameter.
      lambda = length(omega0);
    end
  else
    % randomly initialize the frequency parameters.
    omega0 = unifrnd(0, w(end), lambda, 1);
  end

  % check if an initial decay vector was passed.
  if (exist('rho0'))
    % check that the number of decays matches the expected signal count.
    if (length(rho0) != length(omega0))
      % throw an exception.
      error('rjmcmc1d: decay rate and frequency vector lengths must match');
    end
  else
    % initialize the decay rate parameters.
    rho0 = 1e-1 .* ones(lambda, 1);
  end

  % check if threshold values were passed.
  if (exist('thresh'))
    % check the type of the threshold argument.
    if (!isvector(thresh) || !isreal(thresh) || length(thresh) != 2)
      % throw an exception.
      error('rjmcmc1d: thresholds must be a two-element real vector');
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

  % apodize the FID, just to remove any signal that appears to have a negative
  % decay rate. such signals induce massive artifacts during transformation
  % into frequency space.
  y = sinewindow(t) .* y;

  % calculate the frequency transform of the signal. the NDFT is used, simply
  % because it accepts an output frequency axis at which the values of the
  % spectrum are to be computed. yes, i know it's slower.
  Y = ndft(y, t, w);

  % perform an initial automatic phase correction of the signal.
  [Y, phi0, phi1] = autophase1d(Y, @simplex_entropy);

  % use the phase correction gained from the spectrum to phase the signal
  % in the time domain. this ensures mostly absorptive peaks in the real
  % spectrum, which allows the sampler to discard phase angles during
  % its random walk (and avoid adding dispersive signals).
  y = phase1d(y, phi0, 0);

  % extract the absolute value of the spectrum.
  Yabs = abs(Y);

  % run the compiled RJ-MCMC subroutine.
  [omega, rho, a, K, omega_map, rho_map, a_map] = ...
    __rjmcmc1d(t, y, w, Yabs, n_burn, n_sim, n_mod, ...
               lambda, delta, omega0, rho0, a0, ...
               verbose);

  % sort the MAP parameters by frequency.
  [omega_map, idx] = sort(omega_map);
  rho_map = rho_map(idx);
  a_map = a_map(idx);

  % calculate the magnitude and phase from the complex amplitude. it should
  % really only be a bit off from the 'phase-free' amplitudes we receive from
  % the main subroutine.
  D = __rjmcmc1d_basis(t, omega_map, rho_map);
  a_fit = D \ y;
  phi_map = (180 / pi) .* imag(log(a_fit ./ abs(a_fit)));

  % finalize the output peak table.
  P = [a_map, omega_map, rho_map, phi_map];
end

