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
## @anchor{backscaleplot}
## @deftypefn {Function File} {} backscaleplot (@var{ab}, @var{mdl})
## @deftypefnx {Function File} {} backscaleplot (@var{ab}, @var{mdl}, @var{coloring})
## @deftypefnx {Function File} {@var{pdata} =} backscaleplot (@var{ab}, @var{mdl}, @var{coloring})
## Builds a backscaled loadings plot from PCA, PLS or OPLS modeled data. An
## optional second argument @var{coloring} may be set to enable or disable
## coloring of the backscaled loadings values. The default behavior is to
## generate uncolored plots, because plot coloring can take a while.
##
## Because external plotting programs like @code{gnuplot} are considerably
## quicker at plotting multicolored lines, colored backscaled data ready
## for plotting may be returned into @var{pdata}.
##
## More information on the technique of backscaling model loadings can be
## found in the following two references:
##
## @quotation
## O. Cloarec et. al., `Evaluation of the Orthogonal Projection on Latent
## Structure Model Limitations Caused by Chemical Shift Variability and
## Improved Visualization of Biomarker Changes in 1H NMR Spectroscopic
## Metabonomic Studies', Anal. Chem., 2005(77): 517--526.
## @end quotation
##
## @quotation
## O. Cloarec et. al., `Statistical Total Correlation Spectroscopy: An
## Exploratory Approach for Latent Biomarker Identification from Metabolic
## 1H NMR Data Sets', Anal. Chem., 2005(77): 1282--1289.
## @end quotation
## @end deftypefn

function pdata = backscaleplot (ab, mdl, coloring)
  % use a default number of colormap levels.	
  nmap = 256;

  % check for proper arguments.
  if (!any(nargin == [2 : 3]) || !isstruct(mdl))
    % improper arguments. print the usage statement.
    print_usage();
  end

  % check if the model contains enough components.
  if (mdl.A < 1)
    % nope. throw an exception.
    error('backscaleplot: at least one model dimension required');
  end

  % check if a custom coloring method was selected.
  if (nargin >= 3 && !isempty(coloring))
    % ensure the second argument is a boolean.
    if (!isbool(coloring))
      % nope. throw an exception.
      error('backscaleplot: coloring argument must be a boolean');
    end
  else
    % no. default to no coloring.
    coloring = false;
  end

  % extract the backscaled loadings from the model.
  [B, center, scale] = backscale(mdl);

  % get the size of the loadings matrix.
  [N, K] = size(B);

  % subtract the mean back out from the loadings.
  B -= ones(N, 1) * center;

  % determine if the model is multiblock.
  if (ismultiblock(mdl))
    % build a normalization vector that brings the blocks into the same range
    % of relative intensities.
    k = 1;
    Bnorm = [];
    snorm = [];
    for b = 1 : mdl.B
      % build the normalization for this block of B.
      nval = max(vec(B(:, k : k + mdl.blocks{b}.K - 1)));
      Bnorm = [Bnorm, repmat(nval, 1, mdl.blocks{b}.K)];
      % build the normalization for this block of scale.
      nval = max(scale(k : k + mdl.blocks{b}.K - 1));
      snorm = [snorm, repmat(nval, 1, mdl.blocks{b}.K)];

      % move to the next block.
      k += mdl.blocks{b}.K;
    end

    % normalize the loadings and scale factor by the computed vector.
    B = B .* ((ones(N, 1) * Bnorm) .^ -1);
    scale = scale .* (snorm .^ -1);
  end

  % calculate the separation between rows.
  if (N > 1)
    sep = mean(range(B, 2)) / 2;
  else
    sep = 0;
  end

  % calculate the row-separated B matrix.
  Bsep = linspace(0, N * sep, N)' * ones(1, K) + B;

  % check if an empty abscissa was passed.
  if (isempty(ab) || !isvector(ab))
    % build an 'uninformative' abscissa vector.
    ab = linspace(1, mdl.K, mdl.K)';
  end

  % see if an output argument was requested.
  if (nargout > 0)
    % add the backscaled matrix into the output data matrix.
    pdata = [ab, Bsep'];
  else
    % no output requested. initialize the figure.
    figure();
    hold on;
    title('Loadings plot: backscaled');
    xlabel('p');
  end

  % act based on whether the plot is to be colored.
  if (coloring == true)
    % build the colormap.
    cmap = hsv(floor(1.5 * nmap));
    cmap = cmap(1 : nmap, :);
    cmap = flipud(cmap);

    % calculate the weight values for coloring purposes.
    Smin = min(scale);
    Srange = range(scale);
    w = (scale - Smin) ./ Srange;

    % calculate the color values from the weights.
    cidx = floor(w .* 255 + 1);

    % build the coloring matrix.
    c = cmap(cidx, :);

    % see if an output argument was requested.
    if (nargout > 0)
      % add the coloring values to the output data matrix.
      pdata = [pdata, w'];
    else
      % loop through the variables.
      k = 1;
      len = 1;
      while (k < K)
        % elongate until we find a differently colored variable.
        while (k + len + 1 <= K && cidx(k + len + 1) == cidx(k)), len++; end

        % build the current segment.
        abseg = ab(k : k + len)';
        Bseg = Bsep(:, k : k + len)';

        % plot the data.
        plot(abseg, Bseg, 'color', c(k, :));

        % move onto the next segment.
        k += len;
        len = 1;
      end
    end
  else
    % see if an output argument was requested.
    if (nargout == 0)
      % plot without coloring.
      plot(ab, Bsep);
    end
  end

  % release the figure for plotting only if we made one.
  if (nargout == 0)
    hold off;
  end
end

