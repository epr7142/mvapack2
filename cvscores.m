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
## @anchor{cvscores}
## @deftypefn {Function File} {@var{T} = } cvscores (@var{mdl})
## @deftypefnx {Function File} {@var{T} = } cvscores (@var{mdl}, @var{n})
## Returns the cross-validated scores from a PLS model in the form of an
## array. Each element of the array contains the reconstructed scores
## from a monte carlo leave-n-out cross validation. The length of the
## array equals the number of cross-validation iterations.
## @end deftypefn

function T = cvscores (mdl, n)
  % check for proper arguments.
  if (nargin < 1)
    % improper arguments. print the usage statement.
    print_usage();
  end

  % check if a count was passed.
  if (nargin < 2)
    % set the default value.
    n = [];
  end

  % make sure the model argument is a structure.
  if (!isstruct(mdl))
    % it isn't! throw an exception.
    error('cvscores: model argument must be a structure type');
  end

  % check that the model has a type.
  if (!isfield(mdl, 'type'))
    % throw an error.
    error('cvscores: model has no type');
  end

  % check that the model has a component count.
  if (!isfield(mdl, 'A') || !isscalar(mdl.A))
    % throw an error.
    error('cvscores: model has no component count');
  end

  % check that the model has a cross-validation structure.
  if (!isfield(mdl, 'cv') || !isfield(mdl, 'CV') || isempty(mdl.CV))
    % throw an error.
    error('cvscores: model has no cross-validation information');
  end

  % run based on the model type.
  if (strcmp(mdl.type, 'pls'))
    % check if we need a default component count.
    if (isempty(n))
      n = mdl.A;
    end

    % return the pls cv-scores.
    T = cvplsscores(mdl, n);
  elseif (strcmp(mdl.type, 'opls'))
    % check if we need a default component count.
    if (isempty(n))
      n = mdl.Ap;
    end

    % return the opls cv-scores.
    T = cvoplsscores(mdl, n);
  else
    % unsupported model type.
    error('cvscores: only pls and opls models are supported');
  end
end
