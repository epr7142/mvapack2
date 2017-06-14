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
## @anchor{scoresplot}
## @deftypefn {Function File} {} scoresplot (@var{mdl})
## @deftypefnx {Function File} {} scoresplot (@var{mdl}, @var{d})
## @deftypefnx {Function File} {} scoresplot (@var{mdl}, @var{d}, @var{coloring})
## @deftypefnx {Function File} {} scoresplot (@var{mdl}, @var{d}, @var{coloring}, @var{numbers})
## Builds a scores plot of PCA, PLS or OPLS modeled data. The number of
## components to plot may be supplied as the optional second argument @var{d}.
## An optional matrix @var{coloring} may be passed to define the color scheme
## of the plot. An optional third argument (@var{numbers}) may be set to true
## (default is false) to change the plotted points into observation numbers.
## @end deftypefn

function scoresplot (mdl, d, coloring, numbers)
  % make sure we have at least one argument.
  if (nargin < 1)
    % print the usage statement.
    print_usage();
  end

  % make sure the model argument is a structure.
  if (!isstruct(mdl))
    % it isn't! throw an exception.
    error('scoresplot: model argument must be a structure type');
  end

  % determine whether the data is two- or three-dimensional.
  if (nargin < 2 || isempty(d))
    % auto-detect the dimensionality based on the number of components.
    if (mdl.A == 2)
      % set 2D data.
      d = 2;
    elseif (mdl.A == 3)
      % set 3D data.
      d = 3;
    else
      % yeah, this will probably always happen for opls models...
      error('scoresplot: could not auto-detect number of dimensions');
    end
  end

  % determine if a coloring scheme was passed.
  if (nargin < 3 || isempty(coloring))
    % set the default coloring scheme based on whether the model contains a
    % class matrix.
    if (isfield(mdl, 'Y'))
      % yes class matrix. use class coloring.
      coloring = @clscolors;
    else
      % no class matrix. use observation coloring.
      coloring = @obscolors;
    end
  else
    % make sure the coloring argument is a function handle.
    if (!is_function_handle(coloring))
      % whoops. throw an exception.
      error('scoresplot: the coloring argument must be a function handle');
    end

    % make sure the model contains a class matrix if class coloring was
    % requested.
    if (coloring == @clscolors && !isfield(mdl, 'Y'))
      % the matrix must exist! throw an exception.
      error('scoresplot: clscolors requires that the model contain classes');
    end
  end

  % determine if numbering was requested.
  if (nargin < 4 || isempty(numbers))
    % use the default behavior: don't plot numbers.
    numbers = false;
  else
    % ensure the passed argument is a boolean.
    if (!isbool(numbers))
      % it isn't. throw an exception.
      error('scoresplot: numbers argument must be a boolean');
    end
  end

  % make sure the model contains at least the number of requested components.
  if (mdl.A < d)
    % it doesn't. throw an exception.
    error('scoresplot: model contains an insufficient number of components');
  end

  % see which dimensionality of plotting we'll be doing.
  if (d == 2)
    % run the 2D function.
    scoresplot2(mdl, coloring, numbers);
  elseif (d == 3)
    % run the 3D function.
    scoresplot3(mdl, coloring, numbers);
  else
    error('scoresplot: unable to plot desired dimensionality of data');
  end
end

