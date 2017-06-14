## Copyright (C) 2015 University of Nebraska Board of Regents.
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
## @anchor{press}
## @deftypefn {Function File} {@var{S} =} press (@var{mdl}, @var{Yhat})
## @deftypefnx {Function File} {[@var{S}, @var{Ylim}] =} press (@var{mdl}, @var{Yhat})
## Calculates the predicted residual sum of squares (PRESS) statistic
## used in determining cross-validation model reliability metrics.
##
## This function is used by @ref{cvpls} and @ref{cvopls}, and automatically
## decides to compute a standard PRESS or a discriminant PRESS based on
## whether the model in @var{mdl} is a discriminant analysis (DA) model
## or not. If PRESSD is computed, the bounded response matrix is returned
## in @var{Ylim}.
## @end deftypefn

function [S, Ylim] = press (mdl, Yhat)
  % check if the model is of the discriminant type.
  if (mdl.isda == true && mdl.M == 1)
    % pull the original unscaled class matrix.
    Y0 = mdl.Y0;
    Y = mdl.Y;

    % backscale the predicted class matrix.
    Y0hat = ones(mdl.N, 1) * mdl.mean.Y + Yhat * diag(mdl.scale.Y);
    Ylim = Yhat;

    % bound the backscaled matrix.
    idx = find((Y0 == 1) .* (Y0hat > 1) + (Y0 == 0) .* (Y0hat < 0));
    Ylim(idx) = Y(idx);

    % compute the bounded residual matrix.
    S = mdl.Y - Ylim;
  else
    % return the straight up difference. do not modify the responses.
    S = mdl.Y - Yhat;
    Ylim = Yhat;
  end
end

