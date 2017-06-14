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
## @anchor{acquparms_agilent}
## @deftypefn {Function File} {@var{p} =} acquparms_agilent (@var{dirname}, @var{fbase}, @var{dim})
## Reads values from the key-value pairs found in the Agilent 'procpar' file
## from a filename provided by the @ref{acquparms} function. It is highly
## recommended that you use @ref{acquparms} instead of this function.
## @end deftypefn

function p = acquparms_agilent (dirname, fbase, dim)
  % check for proper arguments.
  if (nargin != 3 || nargout != 1 || !ischar(dirname) || !ischar(fbase))
    % improper arguments. print the usage statement.
    print_usage();
  end

  % build the acquisition parameter filename.
  fname = [dirname, '/', fbase];

  % set the nmr data source type.
  p.type = 'agilent';
  p.dim = dim;

  % set the default group delay parameter.
  p.grpdly = -1;

  % open the file.
  fh = fopen(fname, 'rb');

  % loop through the file until we've reached the end.
  while (!feof(fh))
    % get the next line of the file and tokenize it based on ' ' delimiters.
    str = fgetl(fh);
    kv = strsplit(str, ' ');

    % see if the tokenized line is of the parameter header type.
    if (length(kv) == 11)
      % trim the key string of extra nonsense characters.
      k = strtrim(kv{1});

      % read in the next line of parameters and split them.
      vv = strsplit(fgetl(fh), ' ');

      % determine which key string was parsed.
      if (strcmp(k, 'np') == 1)
        % original number of data points.
        tdorig = str2num(strtrim(vv{2}));

        % actual number of data points.
        p.td = pow2(prevpow2(tdorig));
      elseif (strcmp(k, 'sfrq') == 1)
        % carrier base frequency (mhz).
        p.obs = str2num(strtrim(vv{2}));
      elseif (strcmp(k, 'sw') == 1)
        % spectral width (hz).
        p.sw.hz = str2num(strtrim(vv{2}));
      elseif (strcmp(k, 'reffrq') == 1)
        % reference frequency (mhz).
        reffrq = str2num(strtrim(vv{2}));
      end
    end
  end

  % close the file.
  fclose(fh);

  % calculate the carrier offset (hz).
  p.car.hz = (p.obs - reffrq) * 1.0e6;
end

