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
## @anchor{acquparms_bruker}
## @deftypefn {Function File} {@var{p} =} acquparms_bruker (@var{dirname}, @var{fbase}, @var{dim})
## Reads values from the key-value pairs found in the Bruker 'acqus' file
## from a filename provided by the @ref{acquparms} function. It is highly
## recommended that you use @ref{acquparms} instead of this function.
## @end deftypefn

function p = acquparms_bruker (dirname, fbase, dim)
  % check for proper arguments.
  if (nargin != 3 || nargout != 1 || !ischar(dirname) || !ischar(fbase))
    % improper arguments. print the usage statement.
    print_usage();
  end

  % build the acquisition parameter filename.
  fname = [dirname, '/', fbase];

  % set the nmr data source type.
  p.type = 'bruker';
  p.dim = dim;

  % set the default group delay parameter.
  p.grpdly = -1;

  % initially assume that the data is uniformly sampled.
  nustd = 0;

  % open the file.
  fh = fopen(fname, 'rb');

  % loop through the file until we've reached the end.
  while (!feof(fh))
    % get the next line of the file and tokenize it based on '=' delimiters.
    str = fgetl(fh);
    kv = strsplit(str, '=');

    % see if the tokenized line is of the key=value type.
    if (length(kv) == 2)
      % trim the key string of extra nonsense characters.
      k = strtrim(kv{1});
      k = k(4:end);

      % determine which key string was parsed.
      if (strcmp(k, 'DECIM') == 1)
        % decimation ratio.
        p.decim = str2num(strtrim(kv{2}));
      elseif (strcmp(k, 'GRPDLY') == 1)
        % group delay.
        p.grpdly = str2num(strtrim(kv{2}));
      elseif (strcmp(k, 'DSPFVS') == 1)
        % dsp firmware version string.
        p.dspfvs = str2num(strtrim(kv{2}));
      elseif (strcmp(k, 'TD') == 1)
        % number of data points.
        p.td = str2num(strtrim(kv{2}));
      elseif (strcmp(k, 'NusTD') == 1)
        % if this value is nonzero and greater than TD, then
        % the data has been nonuniformly sampled.
        nustd = str2num(strtrim(kv{2}));
      elseif (strcmp(k, 'SFO1') == 1)
        % carrier base frequency (mhz).
        p.obs = str2num(strtrim(kv{2}));
      elseif (strcmp(k, 'SW_h') == 1)
        % spectral width (hz).
        p.sw.hz = str2num(strtrim(kv{2}));
      elseif (strcmp(k, 'O1') == 1)
        % carrier offset (hz).
        p.car.hz = str2num(strtrim(kv{2}));
      elseif (strcmp(k, 'NUC1') == 1)
        % nucleus.
        p.nuc = strtrim(kv{2});

        % trim off the brackets from the nucleus string.
        [ss, se] = regexp(p.nuc, '[^<>]+');
        p.nuc = p.nuc(ss:se);
      end
    end
  end

  % close the file.
  fclose(fh);

  % check if 'NusTD' was set in the acqus file.
  if (nustd > 0 && nustd != p.td)
    % store the 'uniform' total data point count.
    p.tdnus = nustd;

    % store the nus status as true.
    p.nus = true;

    % check if a sampling schedule exists.
    if (fexist([dirname, '/nuslist']))
      % load the schedule matrix and extract the relevant column.
      sch = load([dirname, '/nuslist']);
      sch = sch(:, dim - 1);

      % store the schedule vector in the parameter structure.
      p.sched = sch;
    end
  else
    % store the nus status as false.
    p.nus = false;
  end
end

