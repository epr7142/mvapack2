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
## @anchor{savescores}
## @deftypefn {Function File} {} savescores (@var{mdl}, @var{filename})
## @deftypefnx {Function File} {} savescores (@var{mdl}, @var{filename}, @var{A})
## @deftypefnx {Function File} {} savescores (@var{mdl}, @var{filename}, @var{A}, @var{Y})
## @deftypefnx {Function File} {} savescores (@var{mdl}, @var{filename}, @var{A}, @var{Y}, @var{labels})
## Exports scores from a PCA, PLS or OPLS model to SIMCA-P+ format text.
## @end deftypefn

function savescores (mdl, filename, A, Y, labels)
  % check for proper arguments..
  if (nargin < 2 || !isstruct(mdl) || !ischar(filename))
    % improper arguments. print the usage statement.
    print_usage();
  end

  % check if a third argument wasn't provided.
  if (nargin < 3 || isempty(A))
    % see if the model has a component count.
    if (isfield(mdl, 'A'))
      % use the model's component count.
      A = mdl.A;
    else
      % invalid.
      error('savescores: model contains no component count. must be wrong.');
    end
  end

  % ensure the requested component count is met by the model.
  if (!isfield(mdl, 'A') || mdl.A < A)
    % invalid.
    error('savescores: at least two model dimensions required');
  end

  % check if the class matrix is provided.
  if (nargin < 4 || isempty(Y))
    % not provided. check if it should have been provided.
    if (strcmp(mdl.type, 'pca') == 1 && !isfield(mdl, 'Y'))
      % yes. throw an exception.
      error('savescores: class matrix (Y) required for unsupervised data');
    else
      % no. use the classes embedded in the model.
      Y = backscaleclasses(mdl);
    end
  end

  % check if a fifth argument was provided.
  if (nargin < 5 || isempty(labels))
    % see if the model contains labels.
    if (isfield(mdl, 'labels'))
      % use the labels from the model.
      labels = mdl.labels;
    else
      % invalid.
      error('savescores: model contains no labels');
    end
  end

  % extract the requested number of scores from the model.
  T = scores(mdl, A);

  % open the output file.
  fh = fopen(filename, 'wb');
  fprintf(fh, 'Obs ID (Primary)\tObs ID (Obs. Sec. ID:1)');

  % loop through the components, outputting the header.
  for a = 1 : A
    fprintf(fh, '\tM1.t[%d]', a);
  end
  fprintf(fh, '\n');

  % loop through the observations, outputting the data.
  for i = 1 : mdl.N
    % output the first two values.
    fprintf(fh, '%d\t%s', i, labels{find(Y(i,:))});

    % loop through the components, outputting the values.
    for a = 1 : A
      fprintf(fh, '\t%f', T(i,a));
    end
    fprintf(fh, '\n');
  end

  % close the output file.
  fclose(fh);
end

