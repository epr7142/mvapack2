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
## @anchor{susplot}
## @deftypefn {Function File} {} susplot (@var{mdl1}, @var{mdl2})
## @deftypefnx {Function File} {} susplot (@var{mdl1}, @var{mdl2}, @var{a})
## @deftypefnx {Function File} {} susplot (@var{mdl1}, @var{mdl2}, @var{a}, @var{numbers})
## @deftypefnx {Function File} {@var{pdata} =} susplot (@dots{})
## Builds a Shared and Unique Structure (SUS)-plot from (O)PLS modeled data
## to facilitate identification of variables that covary and contravary in
## two models. An optional third argument @var{a} may be passed to select
## which predictive component is to be analyzed. The default behavior is to
## plot the first predictive component. An optional fourth argument
## (@var{numbers}) may be set to true (default is false) to change the
## plotted points into variable numbers.
## @end deftypefn

function pdata = susplot (mdl1, mdl2, a, numbers)
  % check for proper arguments.
  if (nargin < 2 || !isstruct(mdl1) || !isstruct(mdl2))
    % improper arguments. print the usage statement.
    print_usage();
  end

  % check if a second argument was passed.
  if (nargin < 3 || isempty(a))
    % use the default component index.
    a = 1;
  else
    % ensure the component index is a real scalar.
    if (!isscalar(a) || !isreal(a))
      % throw an exception.
      error('susplot: component index argument (a) must be a real scalar');
    end
  end

  % check if numbers were requested.
  if (nargin < 4 || !isbool(numbers))
    % default to plotting points.
    numbers = false;
  end

  % check that the model has a type field.
  if (!isfield(mdl1, 'type') || !isfield(mdl2, 'type'))
    % invalid models. throw an exception.
    error('susplot: one or both models lacks a type field');
  end

  % check that the model variable counts match.
  if (mdl1.K != mdl2.K)
    % this won't work. throw an exception.
    error('susplot: model variable counts does not match');
  end

  % pull out some data from the first model based on its type field.
  if (strcmp(mdl1.type, 'pls'))
    % extract the scores and component count.
    T1 = mdl1.T;
    A1 = mdl1.A;

    % set the score label.
    t1lbl = 't';
  elseif (strcmp(mdl1.type, 'opls'))
    % extract the predictive scores and component count.
    T1 = mdl1.Tp;
    A1 = mdl1.Ap;

    % set the score label.
    t1lbl = 't_p';
  else
    % invalid model type. throw an exception.
    error('susplot: first model type must be either pls or opls');
  end

  % pull out some data from the second model based on its type field.
  if (strcmp(mdl2.type, 'pls'))
    % extract the scores and component count.
    T2 = mdl2.T;
    A2 = mdl2.A;

    % set the score label.
    t2lbl = 't';
  elseif (strcmp(mdl2.type, 'opls'))
    % extract the predictive scores and component count.
    T2 = mdl2.Tp;
    A2 = mdl2.Ap;

    % set the score label.
    t2lbl = 't_p';
  else
    % invalid model type. throw an exception.
    error('susplot: second model type must be either pls or opls');
  end

  % check if the first model contains the desired component.
  if (A1 < a)
    % nope. throw an exception.
    error('susplot: first model does not contain %d components', a);
  end

  % check if the second model contains the desired component.
  if (A2 < a)
    % nope. throw an exception.
    error('susplot: second model does not contain %d components', a);
  end

  % extract the desired score vectors from the models.
  t1 = T1(:, a);
  t2 = T2(:, a);

  % extract the unscaled data matrices from the models.
  X1 = mdl1.X0;
  X2 = mdl2.X0;

  % calculate the sus-plot values.
  s1 = (t1' * X1) ./ (mdl1.N - 1) ./ (std(t1) .* std(X1));
  s2 = (t2' * X2) ./ (mdl2.N - 1) ./ (std(t2) .* std(X2));

  % see if the data was requested.
  if (nargout < 1)
    % initialize the figure.
    figure();
    hold on;
    title(sprintf('SUS-plot: corr(%s, X)[1] vs corr(%s, X)[2]', t1lbl, t2lbl));
    xlabel(sprintf('corr(%s, X)[1]', t1lbl));
    ylabel(sprintf('corr(%s, X)[2]', t2lbl));

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

