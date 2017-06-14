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
## @anchor{splot}
## @deftypefn {Function File} {} splot (@var{mdl})
## @deftypefnx {Function File} {} splot (@var{mdl}, @var{a})
## @deftypefnx {Function File} {} splot (@var{mdl}, @var{a}, @var{numbers})
## @deftypefnx {Function File} {@var{pdata} =} splot (@dots{})
## Builds an S-plot from (O)PLS modeled data to facilitate identification of
## variables that contribute strongly to class distinction. An optional
## second argument @var{a} may be passed to select which predictive component
## is to be analyzed. The default behavior is to plot the first predictive
## component. An optional third argument (@var{numbers}) may be set to true
## (default is false) to change the plotted points into variable numbers.
## @end deftypefn

function pdata = splot (mdl, a, numbers)
  % check for proper arguments.
  if (nargin < 1 || !isstruct(mdl))
    % improper arguments. print the usage statement.
    print_usage();
  end

  % check if a second argument was passed.
  if (nargin < 2 || isempty(a))
    % use the default component index.
    a = 1;
  else
    % ensure the component index is a real scalar.
    if (!isscalar(a) || !isreal(a))
      % throw an exception.
      error('splot: component index argument (a) must be a real scalar');
    end
  end

  % check if numbers were requested.
  if (nargin < 3 || !isbool(numbers))
    % default to plotting points.
    numbers = false;
  end

  % check that the model has a type field.
  if (!isfield(mdl, 'type'))
    % invalid model. throw an exception.
    error('splot: model does not contain a type field');
  end

  % pull out some data from the model based on its type field.
  if (strcmp(mdl.type, 'pls'))
    % extract the scores and component count.
    T = mdl.T;
    A = mdl.A;

    % set the score label.
    tlbl = 't';
  elseif (strcmp(mdl.type, 'opls'))
    % extract the predictive scores and component count.
    T = mdl.Tp;
    A = mdl.Ap;

    % set the score label.
    tlbl = 't_p';
  else
    % invalid model type. throw an exception.
    error('splot: model type must be either pls or opls');
  end

  % check if the model contains the desired component.
  if (A < a)
    % nope. throw an exception.
    error('splot: model does not contain %d components', a);
  end

  % extract the desired score vector from the model.
  t = T(:, a);

  % extract the unscaled data matrix from the model.
  X = mdl.X0;

  % calculate the s-plot values.
  s1 = (t' * X) ./ (mdl.N - 1);
  s2 = s1 ./ (std(t) .* std(X));

  % see if the data was requested.
  if (nargout < 1)
    % initialize the figure.
    figure();
    hold on;
    title(sprintf('S-plot: cov(%s, X) vs corr(%s, X)', tlbl, tlbl));
    xlabel(sprintf('cov(%s, X)', tlbl));
    ylabel(sprintf('corr(%s, X)', tlbl));

    % see if numbers or points are desired.
    if (numbers == true)
      % white out the score value points.
      scatter(s1, s2, [], [1, 1, 1], '^');

      % loop through the loadings.
      for k = 1 : length(s1)
        % print the current variable index.
        text(s1(k), s2(k), num2str(k), 'color', [0, 0, 0]);
      end
    else
      % scatter plot the score values.
      scatter(s1, s2, [], [0, 0, 0], '^');
    end

    % release the figure for plotting.
    hold off;
  else
    % return the generated s-plot data instead of plotting.
    pdata = [[1 : length(s1)]; s1; s2]';
  end
end

