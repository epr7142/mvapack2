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
## @anchor{responses}
## @deftypefn {Function File} {[@var{Y}, @var{Yfit}, @var{Yhat}] =} responses (@var{mdl})
## Returns the univariate backscale response values from a PLS or OPLS
## model. The outputs @var{Y}, @var{Yfit} and @var{Yhat} contain the input
## response values, output fitted responses and output re-estimated responses.
## @end deftypefn

function [Y, Yfit, Yhat] = responses (mdl)
  % check for proper arguments.
  if (nargin != 1 || !any(nargout == [1 : 3]) || !isstruct(mdl))
    % improper arguments. print the usage statement.
    print_usage();
  end

  % see if the model is of the correct type.
  if (!strcmp(mdl.type, 'pls') && !strcmp(mdl.type, 'opls'))
    % throw an exception.
    error('responses: invalid model type');
  end

  % retrieve the input data matrix from the model.
  X = mdl.X;

  % if the model is an opls model, remove the orthogonal variation.
  if (strcmp(mdl.type, 'opls'))
    % yes, remove the orthogonal components.
    X -= mdl.To * mdl.Po';
  end

  % return the original responses from the model.
  Y = backscale(mdl.Y, mdl.mean.Y, mdl.scale.Y);

  % backscale the responses for output.
  Yfit = backscale(mdl.U * mdl.C', mdl.mean.Y, mdl.scale.Y);
  Yhat = backscale(X * mdl.B, mdl.mean.Y, mdl.scale.Y);
end

