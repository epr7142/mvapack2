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
## @anchor{scoresplot3}
## @deftypefn {Function File} {} scoresplot3 (@var{mdl}, @var{coloring}, @var{numbers})
## Builds a three-dimensional scores plot of PCA, PLS or OPLS modeled data.
## It is recommended that you not use this function directly. Use
## @ref{scoresplot} instead and specify three components.
## @end deftypefn

function scoresplot3 (mdl, coloring, numbers)
  % get the three required score vectors from the model.
  T = scores(mdl, 3);

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
    colors = coloring(T);
  end

  % read in the model quality parameters.
  qual = rqdiff(mdl);

  % begin a new figure.
  figure();
  hold on;

  % set the figure title.
  title('Scores plot: t_1 vs t_2 vs t_3');

  % see if quality values are available.
  if (rows(qual) >= 3)
    % build the figure axis labels.
    xstr = sprintf('t_1 (R^2=%.3f, Q^2=%.3f)', qual(1,1), qual(1,2));
    ystr = sprintf('t_2 (R^2=%.3f, Q^2=%.3f)', qual(2,1), qual(2,2));
    zstr = sprintf('t_3 (R^2=%.3f, Q^2=%.3f)', qual(3,1), qual(3,2));
  else
    % build the labels without values.
    xstr = 't_1';
    ystr = 't_2';
    zstr = 't_3';
  end

  % set the figure axis labels.
  xlabel(xstr);
  ylabel(ystr);
  zlabel(zstr);

  % turn off hidden line removal.
  hidden('off');

  % see if the figure will be... classeh.
  if (classy == true)
    % see if numbers or points are desired.
    if (numbers == true)
      % plot white points.
      scatter3(T(:,1), T(:,2), T(:,3), [], ones(size(colors)), '^');

      % loop through the observations.
      for i = 1 : rows(T)
        % print the current observation index.
        text(T(i,1), T(i,2), T(i,3), num2str(i), 'color', colors(i,:));
      end
    else
      % yes indeed we have class.
      for m = 1 : columns(Y)
        % extract the scores and ellipse for the current class.
        idx = classidx(Y, m);
        Tcls = T(idx,:);
        Ecls = ellipse(Tcls);

        % plot the ellipsoid for the current class.
        plot3(Ecls.x(:,1), Ecls.y(:,1), Ecls.z(:,1), ...
              'color', clsmap(m,:), '-');
      end

      % see if the model contains a labels field.
      if (isfield(mdl, 'labels'))
        % print the legend labelling.
        legend(mdl.labels);
      end

      % loop again through the classes.
      for m = 1 : columns(Y)
        % extract the scores for the current class.
        idx = classidx(Y, m);
        Tcls = T(idx,:);
        Ecls = ellipse(Tcls);

        % plot the ellipsoid for the current class.
        plot3(Ecls.x, Ecls.y, Ecls.z, 'color', clsmap(m,:), '-');
        plot3(Ecls.x', Ecls.y', Ecls.z', 'color', clsmap(m,:), '-');

        % plot the points for the current class.
        scatter3(Tcls(:,1), Tcls(:,2), Tcls(:,3), [], clsmap(m,:), '-');
      end
    end
  else
    % see if numbers or points are desired.
    if (numbers == true)
      % plot white points.
      scatter3(T(:,1), T(:,2), T(:,3), [], ones(size(colors)), '^');

      % loop through the observations.
      for i = 1 : rows(T)
        % print the current observation index.
        text(T(i,1), T(i,2), T(i,3), num2str(i), 'color', colors(i,:));
      end
    else
      % classless. just plot the points.
      scatter3(T(:,1), T(:,2), T(:,3), [], colors, '^');
    end
  end

  % finish the figure.
  hold off;
end

