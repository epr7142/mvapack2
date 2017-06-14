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
## @anchor{cvscoresplot3}
## @deftypefn {Function File} {} cvscoresplot3 (@var{mdl}, @var{coloring}, @var{numbers})
## Builds a three-dimensional scores plot of cross-validated scores from PLS
## or OPLS modeled data. It is recommended that you not use this function
## directly. Use @ref{cvscoresplot} instead and specify two components.
## @end deftypefn

function cvscoresplot3 (mdl, coloring, numbers)
  % get the three required CV score vectors from the model.
  T = cvscores(mdl, 3);

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
    colors = coloring(scores(mdl, 3));
  end

  % read in the model quality parameters.
  qual = rqdiff(mdl);

  % compute ellipse information on a per-observation basis.
  E = cellfun(@(t) ellipse(t), T, 'UniformOutput', false);

  % begin a new figure.
  figure();
  hold on;

  % set the figure title.
  title('CV scores plot: t_{CV,1} vs t_{CV,2} vs t_{CV,3}');

  % see if quality values are available.
  if (rows(qual) >= 2)
    % build the figure axis labels.
    xstr = sprintf('t_{CV,1} (R^2=%.3f, Q^2=%.3f)', qual(1,1), qual(1,2));
    ystr = sprintf('t_{CV,2} (R^2=%.3f, Q^2=%.3f)', qual(2,1), qual(2,2));
    zstr = sprintf('t_{CV,3} (R^2=%.3f, Q^2=%.3f)', qual(3,1), qual(3,2));
  else
    % build the labels without values.
    xstr = 't_{CV,1}';
    ystr = 't_{CV,2}';
    zstr = 't_{CV,3}';
  end

  % set the figure axis labels.
  xlabel(xstr);
  ylabel(ystr);
  zlabel(zstr);

  % turn off hidden line removal.
  hidden('off');

  % see if we should plot class labels.
  if (classy == true && isfield(mdl, 'labels'))
    % loop over each class.
    for m = 1 : mdl.M
      % plot a fake point.
      plot3([0, 0], [1, 1], [2, 2], 'color', clsmap(m,:));
    end

    % print the legend labelling.
    legend(mdl.labels);

    % loop again over each class.
    for m = 1 : mdl.M
      % plot a fake point.
      plot3([0, 0], [1, 1], [2, 2], 'color', [1, 1, 1], 'linewidth', 4);
    end
  end

  % loop through every observation in the model.
  for n = 1 : mdl.N
    % get the observation mean.
    c1 = E{n}.center(1);
    c2 = E{n}.center(2);
    c3 = E{n}.center(3);

    % get the observation ellipse.
    e1 = E{n}.x;
    e2 = E{n}.y;
    e3 = E{n}.z;

    % check if we're classy.
    if (classy == true)
      % determine class membership for the current observation.
      m = find(Y(n,:));

      % plot a class-colored ellipsoid.
      plot3(e1, e2, e3, 'color', clsmap(m,:), '-');
      plot3(e1', e2', e3', 'color', clsmap(m,:), '-');

      % are we plotting numbers or points?
      if (numbers == true)
        % print a class-colored number.
        text(c1, c2, c3, num2str(n), 'color', clsmap(m,:));
      else
        % plot a class-colored point.
        scatter3(c1, c2, c3, [], clsmap(m,:), '-');
      end
    else
      % plot a simpler ellipse.
      plot3(e1, e2, e3, 'color', colors(n,:), '-');
      plot3(e1', e2', e3', 'color', colors(n,:), '-');

      % are we plotting numbers or points?
      if (numbers == true)
        % print a number.
        text(c1, c2, c3, num2str(n), 'color', colors(n,:));
      else
        % plot a point.
        scatter3(c1, c2, c3, [], colors(n,:), '-');
      end
    end
  end

  % finish the figure.
  hold off;
end
