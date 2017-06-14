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
## @anchor{sinewindow}
## @deftypefn {Function File} {@var{wfid} =} sinewindow (@var{t}, @var{opts})
## Calculate window coefficients for sinusoidal windowing of a signal.
##
## The parameters of the window are specified as fields in the @var{opts}
## structure. The parameter @var{offset} is the relative offset of the sine
## bell. The parameter @var{ending} is the relative ending of the sine bell.
## Both @var{offset} and @var{ending} span (@math{[0, 1]}). The final
## parameter @var{exponent} specifies the order of the sinusoid.
##
## The default parameter values are 0.5, 0 and 2, resulting in a squared
## cosine window.
## @tex
##
## The equation applied to the free induction decay is as follows:
## $$ wfid(t) = fid(t) \cdot \sin \left [ \pi (offset) + \pi (ending - offset) 
##    \left ( t \over t_{max} \right ) \right ] $$
## @end tex
##
## Instead of using this function directly, it is recommended that you use
## @ref{apodize}.
## @end deftypefn

function w = sinewindow (t, opts)
  % check for the minimum number of arguments.
  if (!any(nargin == [1 : 2]) || nargout != 1)
    % print the usage statement.
    print_usage();
  end

  % check for a passed option set.
  if (nargin >= 2 && !isempty(opts))
    % yes. check the type of the options.
    if (!isstruct(opts))
      % invalid argument. throw an exception.
      error('sinewindow: window options must be a structure data type');
    end

    % is there an offset option?
    if (isfield(opts, 'offset'))
      % yes. use it.
      offset = opts.offset;
    end

    % is there an ending option?
    if (isfield(opts, 'ending'))
      % yes. use it.
      ending = opts.ending;
    end

    % is there an exponent option?
    if (isfield(opts, 'exponent'))
      % yes. use it.
      exponent = opts.exponent;
    end
  end

  % check for an offset parameter.
  if (exist('offset'))
    % yes. check it.
    if (!isscalar(offset) || !isreal(offset) || offset < 0 || offset > 1)
      % invalid. throw an exception.
      error('sinewindow: invalid offset parameter.');
    end
  else
    % no. use the default.
    offset = 0.5;
  end

  % check for an ending parameter.
  if (exist('ending'))
    % yes. check it.
    if (!isscalar(ending) || !isreal(ending) || ending < 0 || ending > 1)
      % invalid. throw an exception.
      error('sinewindow: invalid ending parameter.');
    end
  else
    % no. use the default.
    ending = 1.0;
  end

  % check for an exponent parameter.
  if (exist('exponent'))
    % yes. check it.
    if (!isscalar(exponent) || !isreal(exponent) || exponent < 1)
      % invalid. throw an exception.
      error('sinewindow: invalid exponent parameter.');
    end
  else
    % no. use the default.
    exponent = 2.0;
  end

  % compute the window coefficients.
  w = sin(pi * offset + pi .* (ending - offset) .* (t - min(t)) ./ range(t));
  w = w .^ exponent;
end

