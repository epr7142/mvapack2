## Copyright (C) 2014 University of Nebraska Board of Regents.
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
## @anchor{ist}
## @deftypefn {Function File} {@var{y} =} ist (@var{x}, @var{sched})
## @deftypefnx {Function File} {@var{y} =} ist (@var{x}, @var{sched}, @var{opts})
## Perform Iterative Soft Thresholding (IST) reconstruction of a vector or
## matrix @var{x} based on the known points in @var{sched}.
##
## The IST algorithm is performed along the first non-singleton dimension of
## the array. Thus if @var{x} is a matrix, @code{ist (x, sched)} performs IST
## for each column of @var{x}. The known point indices in @var{sched} must be
## in vector form.
##
## An optional argument @var{opts} may be provided to specify advanced options.
## For more information on IST options, see @ref{ist_options}.
## @end deftypefn

function y = ist (x, sched, opts)
  % check the arguments.
  if (!any(nargin == [2 : 3]) || nargout != 1)
    % invalid arguments.
    print_usage();
  end

  % check the schedule vector.
  if (isempty(sched) || !isvector(sched) || !isreal(sched))
    % invalid type. throw an exception.
    error('ist: schedule is expected to be a real vector');
  end

  % check the boundaries of the schedule vector.
  if (min(sched) < 1 || max(sched) > rows(x))
    % invalid bounds. throw an exception.
    error('ist: schedule indices are out of bounds');
  end

  % check if options were provided.
  if (nargin < 3 || isempty(opts) || !isstruct(opts))
    % call in the default options.
    opts = ist_options();
  end

  % extract the parameters.
  tau = opts.tau;
  N = opts.iters;

  % initialize the threshold.
  lambda = max(abs(fft(x)));

  % build a vector of sampling weights.
  w = zeros(size(x));
  w(sched) = ones(size(sched));

  % check the type of the input data.
  if (isvector(x))
    % get the length of the input data.
    m = length(x);

    % initialize the output vectors.
    Y = zeros(2 * m, 1);
    y = zeros(m, 1);

    % loop for the number of iterations.
    for n = 1 : N
      % compute the residual over the measured data points.
      r = w .* (x - y);

      % sum the current residual into the output f.d. vector.
      Y += fft([r; zeros(m, 1)]);

      % threshold into the output frequency domain vector.
      Yth = Y .* (abs(Y) > lambda);
      Y = Yth .* (1 - lambda ./ abs(Y));

      % inverse fourier transform the thresholded spectrum.
      y = ifft(Y)(1 : m);

      % update the threshold.
      lambda *= (1 - tau);
    end
  elseif (ismatrix(x))
    % initialize the output matrix.
    y = zeros(size(x));

    % loop through the data columns.
    for n = 1 : columns(x)
      % perform IST on the current column.
      y(:, n) = ist(x(:, n), sched, opts);
    end
  end
end
