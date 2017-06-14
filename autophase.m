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
## @anchor{autophase}
## @deftypefn {Function File} {@var{sp} =} autophase (@var{s}, @var{parms})
## @deftypefnx {Function File} {[@var{sp}, @var{phc0}, @var{phc1}] =} autophase (@var{s}, @var{parms})
## @deftypefnx {Function File} {@var{sp} =} autophase (@var{s}, @var{parms}, @var{objective})
## @deftypefnx {Function File} {[@var{sp}, @var{phc0}, @var{phc1}] =} autophase (@var{s}, @var{parms}, @var{objective})
## Performs automatic phase correction of a one- or two-dimensional NMR
## spectrum or spectral dataset, using a simplex optimization algorithm:
##
## @quotation
## M. Siegel. `The use of the modified simplex method for automatic
## phase correction in Fourier-transform Nuclear Magnetic Resonance
## spectroscopy'. Analytica Chimica Acta, 1981. 133(1981): 103-108.
## @end quotation
##
## For one-dimensional spectra, the algorithm uses an entropy minimization
## objective during optimization (@xref{simplex_entropy}). For two-dimensional
## spectra, a whitening objective is used (@xref{simplex_whiten}).
##
## In the one-dimensional case, data in @var{s} may either be a column
## vector or a data matrix where each spectrum is arranged as
## a row in the matrix.
##
## In the two-dimensional case, data in @var{s} may either be a data matrix
## where each direct-dimension slice is along the rows, or a cell array that
## contains multiple matrices, each having direct-dimension slices along
## its rows.
##
## A parameter structure (or array) must be passed as a second argument.
## The phase correction values @var{phc0} and @var{phc1} may also be
## returned if desired.
## @end deftypefn

function [sp, phc0, phc1] = autophase (s, parms, objective)
  % check for proper arguments.
  if (!any(nargin == [2 : 3]) || !any(nargout == [1 : 3]) || ...
      !(iscell(s) || ismatrix(s) || isvector(s)))
    % improper arguments. print the usage statement.
    print_usage();
  end

  % get the dimensions in the dataset.
  nd = nmrdims(s, parms);

  % check if an objective function was specified.
  if (nargin >= 3 && !isempty(objective))
    % ensure a valid function handle was passed.
    if (!is_function_handle(objective))
      % invalid. throw an exception.
      error('autophase: objective must be a valid function handle');
    end
  else
    % set the default objective based on the data dimensionality.
    if (nd == 1)
      % use the entropy objective.
      objective = @simplex_entropy;
    elseif (nd == 2)
      % use the whitening objective.
      objective = @simplex_whiten;
    end
  end

  % check whether we are phasing one- or two-dimensional data.
  if (nd == 1)
    % run a one-dimensional phase correction.
    [sp, phc0, phc1] = autophase1d(s, objective);
  elseif (nd == 2)
    % run a two-dimensional phase correction.
    [sp, phc0, phc1] = autophase2d(s, objective);
  else
    % throw an exception.
    error('autophase: input data must be one- or two-dimensional');
  end
end

