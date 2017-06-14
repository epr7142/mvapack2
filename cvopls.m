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
## @anchor{cvopls}
## @deftypefn {Function File} {[@var{Qsq}, @var{Qstd}, @var{CV}] =} cvopls (@var{mdl}, @var{V}, @var{Ao})
## Performs internal cross-validation of a PLS or OPLS model and returns
## a @math{Q^2} value for inferring model reliability. The @var{CV}
## cell array returned is a modified version of that found in the passed
## model.
## @end deftypefn

function [Qsq, Qstd, CV] = cvopls (mdl, V, Ao)
  % initialize the output structures.
  Q = zeros(mdl.cv.niter, 1);
  CV = mdl.CV;

  % loop through the number of cross-validation iterations.
  for i = 1 : mdl.cv.niter
    % loop through the number of cross-validation groups.
    for g = 1 : mdl.cv.ngroup
      % calculate a new component for cross-validation.
      [w, t, p, u, c, Wo, To, Po, iter] = ...
        oplscomp(CV{i}.Xt{g}, CV{i}.Yt{g}, V, Ao);

      % store the new predictive component in the cross-validation structure.
      CV{i}.Wp{g} = [CV{i}.Wp{g}, w];
      CV{i}.Pp{g} = [CV{i}.Pp{g}, p];
      CV{i}.C{g} = [CV{i}.C{g}, c];

      % store the new orthogonal components in the cross-validation structure.
      CV{i}.Wo{g} = [CV{i}.Wo{g}, Wo];
      CV{i}.Po{g} = [CV{i}.Po{g}, Po];

      % calculate the deflation-free weights and regression coefficients.
      CV{i}.Wps{g} = CV{i}.Wp{g} * inv(CV{i}.Pp{g}' * CV{i}.Wp{g});
      CV{i}.Wos{g} = CV{i}.Wo{g} * inv(CV{i}.Po{g}' * CV{i}.Wo{g});
      CV{i}.B{g} = CV{i}.Wps{g} * CV{i}.C{g}';

      % deflate the data and output matrices.
      CV{i}.Xt{g} = CV{i}.Xt{g} - t * p' - To * Po';
      CV{i}.Yt{g} = CV{i}.Yt{g} - t * c';

      % estimate the left-out responses.
      Yhat = CV{i}.Xv{g} * CV{i}.B{g};
      CV{i}.Yh(find(CV{i}.idx == g), :) = Yhat;
    end

    % calculate the predictive ability statistic.
    [PRESSi, CV{i}.Yh] = press(mdl, CV{i}.Yh);
    Q(i) = 1 - ssratio(PRESSi, mdl.Y);
  end

  % compute the standard deviation of the computed Qsq values.
  Qstd = std(Q);

  % determine how to report the Qsq values.
  if (mdl.A == 0)
    % this is the first component. return the value verbatim.
    Qsq = mean(Q);
  else
    % this is a subsequent component. return the value with the previous
    % value subtracted out.
    Qsq = mean(Q) - sum(mdl.Qsq.comp);
  end
end

