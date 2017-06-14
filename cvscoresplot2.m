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
## @anchor{cvscoresplot2}
## @deftypefn {Function File} {} cvscoresplot2 (@var{mdl}, @var{coloring}, @var{numbers})
## Builds a two-dimensional scores plot of cross-validated scores from PLS
## or OPLS modeled data. It is recommended that you not use this function
## directly. Use @ref{cvscoresplot} instead and specify two components.
## @end deftypefn

function cvscoresplot2 (mdl, coloring, numbers)
  % get the two required CV score vectors from the model.
  T = cvscores(mdl, 2);

  % initialize the color matrix.
  colors = [0, 0, 0];

  % initialize the class matrix and classy state.
  classy = false;
  Y = [];

  % build the color matrix.
  if (coloring == @clscolors)
    % hell yeah, we're classy.
    classy = true;

    % backscale the embedded class matrix.
    Y = backscaleclasses(mdl);

    % apply the coloring to the class matrix.
    [colors, clsmap] = coloring(Y);
  else
    % apply the coloring based on the scores matrix.
    colors = coloring(scores(mdl, 2));
  end

  % read in the model quality parameters.
  qual = rqdiff(mdl);

  % compute ellipse information on a per-observation basis.
  E = cellfun(@(t) ellipse(t), T, 'UniformOutput', false);

  % begin a new figure.
  figure();
  hold on;

  % set the figure title.
  title('CV scores plot: t_{CV,1} vs t_{CV,2}');

  % see if quality values are available.
  if (rows(qual) >= 2)
    % build the figure axis labels.
    xstr = sprintf('t_{CV,1} (R^2=%.3f, Q^2=%.3f)', qual(1,1), qual(1,2));
    ystr = sprintf('t_{CV,2} (R^2=%.3f, Q^2=%.3f)', qual(2,1), qual(2,2));
  else
    % build the labels without values.
    xstr = 't_{CV,1}';
    ystr = 't_{CV,2}';
  end

  % set the figure axis labels.
  xlabel(xstr);
  ylabel(ystr);

  % see if we should plot class labels.
  if (classy == true && isfield(mdl, 'labels'))
    % loop over each class.
    for m = 1 : mdl.M
      % plot a fake point.
      plot([0, 0], [1, 1], 'color', clsmap(m,:));
    end

    % print the legend labelling.
    legend(mdl.labels);

    % loop again over each class.
    for m = 1 : mdl.M
      % plot a fake point.
      plot([0, 0], [1, 1], 'color', [1, 1, 1], 'linewidth', 4);
    end
  end

  % loop through every observation in the model.
  for n = 1 : mdl.N
    % get the observation mean.
    c1 = E{n}.center(1);
    c2 = E{n}.center(2);

    % get the observation ellipse.
    e1 = E{n}.xy(:,1);
    e2 = E{n}.xy(:,2);

    % check if we're classy.
    if (classy == true)
      % determine class membership for the current observation.
      m = find(Y(n,:));

      % plot a class-colored ellipse.
      plot(e1, e2, 'color', clsmap(m,:), '-');

      % are we plotting numbers or points?
      if (numbers == true)
        % print a class-colored number.
        text(c1, c2, num2str(n), 'color', clsmap(m,:));
      else
        % plot a class-colored point.
        scatter(c1, c2, [], clsmap(m,:), '-');
      end
    else
      % plot a simpler ellipse.
      plot(e1, e2, 'color', colors(n,:), '-');

      % are we plotting numbers or points?
      if (numbers == true)
        % print a number.
        text(c1, c2, num2str(n), 'color', colors(n,:));
      else
        % plot a point.
        scatter(c1, c2, [], colors(n,:), '-');
      end
    end
  end

  % plot the main ellipse.
  E = ellipse(scores(mdl, 2), false);
  plot(E.xy(:,1), E.xy(:,2), 'color', [0, 0, 0], 'linewidth', 2, '-');

  % finish the figure.
  hold off;
end
