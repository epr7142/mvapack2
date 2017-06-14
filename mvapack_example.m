% This Octave example script details and documents (in what is hopefully
% painful detail) the processing and modeling of 1D NMR data using the
% MVAPACK software package for GNU Octave. It is highly recommended that
% you use the latest available Octave interpreter and packages, as the
% majority of the MVAPACK code was written and tested on 3.2 and 3.6
% interpreters. OK, enough boring chit-chat. Let's begin.

% Add the MVAPACK utilities to the function search path.
% --
% This tells Octave where to look for the MVAPACK functions you'll be
% calling later on.
addpath('/opt/mvapack');

% Generate a list of directory names to load into the processing environment.
% --
% This function uses the standard Linux wildcard characters to load in all
% desired directories. To wildcard a single character, use the '?' character.
% To wildcard multiple characters, use the '*' character. Remember to include
% the experiment subdirectory in the path glob.
F.dirs = glob('my-experiment-???/1');

% Load in the FID data.
% --
% This function automatically reads in the FID data to Octave matrices using
% the nmrPipe executable as a backend. You may have to use the 'doswap' flag
% in order to get correct data (This will depend on how modern the console
% you collected the data with was, mainly).
% --
% F: a data matrix of (complex) free induction decays, arranged as rows
% parms: a structure containing the parsed spectral parameters of interest
% t: the time-domain data abscissa, useful/required for later functions
[F.data, F.parms, F.t] = loadnmr(F.dirs);

% Define a class matrix to hold the class memberships.
% --
% This function is used when strict class membership Y-matrices are desired.
% All model building routines may alternatively use continuous data for true
% regression, as opposed to the discriminant analysis performed in this
% example code.
cls.Y = classes([10, 10, 10, 9, 9]);
cls.Ypls = classes([10 + 10 + 10 + 9, 9]);

% Define an array of class labels for human-readability.
% --
% The order of appearance of these labels should match the class assignments
% in Y, which should correspond to class assignments of all observations.
cls.labels = {'A', 'B', 'C', 'D', 'E'};

% Alternatively, load in a labels file.
% --
% This function generates labels, indices and Y. The indices vector is a
% vector that can be used to permute the elements of any data matrix into
% alignment with the class matrix. See below:
cls.filename = 'labelsfile.txt';
[cls.labels, cls.indices, cls.Y] = loadlabels(cls.filename);
F.data = F.data(cls.indices,:);

% Apply an exponential apodization to the FID data.
% --
% This function applies an exponential line-broadening window function to
% the FID data, with a broadening value argument in Hertz. The only reason
% that LB is defined prior to the function call is for clarity; you can
% lump the literal value into the function call instead.
% --
% We don't use apodization, due to some sort of dogma. It shouldn't strictly
% be necessary, and is totally moot if the final spectra are subjected to any
% binning procedures.
F.apod.lb = 0.3;
F.apod.fn = @expwindow;
F.data = apodize(F.data, F.parms, F.apod.fn, F.apod.lb);

% Apply a zero-filling operation to the apodized FID data.
% --
% This function expands the FID data length based on power-of-two scaling
% rules. With only the data matrix passed as an argument, a single zero-fill
% will be performed. An optional second parameter may be passed to specify
% the number of desired zero-fills.
F.nzf = 2;
F.data = zerofill(F.data, F.parms, F.nzf);

% Fourier-transform the apodized, zero-filled FID data to spectra.
% --
% This function just performs the appropriate fast fourier transform and
% circular shifting of the input time-domain data to yield proper NMR
% spectra.
% --
% S: the spectral data, with spectra arranged along the rows
[S.data, S.ppm, S.hz] = nmrft(F.data, F.parms);

% Remove an observation from the dataset.
% --
% If acquisition or sample handling went awry for a sample or samples, removal
% is accomplished with the following function. It is desireable to run the
% command on every data structure that contains the observation. This
% example removes the first observation from the dataset.
cls.Y = rmobs(cls.Y, 1);
Ypls = rmobs(cls.Ypls, 1);
F.data = rmobs(F.data, 1);
S.data = rmobs(S.data, 1);

% Manually phase the dataset.
% --
% This function is useful once phase angle corrections are established, either
% through initial manual or automatic phasing routines not demonstrated here.
% In this example, a zero-order correction of -60 degrees and a first-order
% correction of 0 degrees per ppm are applied uniformly across all spectra.
S.phc0 = -60.0;
S.phc1 = 0.0;
S.data = phase(S.data, F.parms, S.phc0, S.phc1);

% Automatically phase the whole dataset.
% --
% This function is more useful in interactive processing, to get an initial
% estimate of the phase correction values, or just to bring everything into
% rough phase alignment.
[S.data, S.phc0, S.phc1] = autophase(S.data, F.parms);

% Automatically phase the mean of the dataset.
% --
% This function is more useful in interactive processing sessions to get an
% initial estimate of phase corrections that can then be used by the 'phase'
% function in a script like this. Phasing the mean is the best way to get a
% valid initial condition for the 'pscorr' function below.
[Smean, S.phc0, S.phc1] = autophase(mean(S.data)', F.parms);
S.data = phase(S.data, S.parms, S.phc0, S.phc1);

% Phase-scatter correct the dataset.
% --
% This is one of many available normalization routines. It is special because
% it corrects between-spectral variations in phase errors while normalizing.
% --
% OTHER AVAILABLE NORMALIZATION METHODS:
%   csnorm:    Constant Sum
%   pqnorm:    Probabilistic Quotient
%   histmatch: Histogram Matching
%   snv:       Standard Normal Variate
%   mscorr:    Multiplicative Scatter Correction
% --
S.data = pscorr(S.data);

% Reference the entire dataset's chemical shifts.
% --
% This function only modifies the chemical shift axis of the dataset, not the
% intensities values. In other words, this shifts every spectrum in the same
% way, rather than shifting each spectrum to match up.
% --
% In this example, a chemical shift value of -0.15 will be shifted to now be
% a value of 0.0.
S.ref.old = -0.15;
S.ref.new = 0.00;
S.ppm = refadj(S.ppm, S.ref.old, S.ref.new);

% Extract the real values from the spectral data for an initial data matrix.
% --
% The data up to this point has been complex in nature, having both real and
% imaginary components. Once the phasing routines are completed, the
% imaginary part of the data can safely be ignored.
X.data = realnmr(S.data, F.parms);
X.ppm = S.ppm;

% Manually find the indices of the portions of the dataset to remove.
% --
% The variable removal function relies on indices, not abscissa values, to
% perform its function. To get indices that correspond most closely with
% the chemical shift regions we wish to cut, this function is used.
% --
% This example removes three regions from the data: one from the most downfield
% end to 8.7 ppm, one from 4.9 to 4.6 ppm, and one from the most upfield end to
% 0.2 ppm.
i1 = findnearest(X.ppm, max(X.ppm));
i2 = findnearest(X.ppm, 8.7);
i3 = findnearest(X.ppm, 4.9);
i4 = findnearest(X.ppm, 4.6);
i5 = findnearest(X.ppm, 0.2);
i6 = findnearest(X.ppm, min(X.ppm));

% Remove the identified index ranges from the dataset.
% --
% Both the data matrix (X) and the abscissa vector (ppm) are modified in
% this process.
[X.data, X.ppm] = rmvar(X.data, X.ppm, [i1 : i2, i3 : i4, i5 : i6]);

% Optimally bin the reduced data matrix.
% --
% This function integrates differing-width chemical shift regions of the data
% into bins. The bin width is declared separately again for clarity. The
% function returns bin intensities, bin centers (optionally), and bin widths
% (optionally).
B.wbin = 0.025;
[B.data, B.ppm, B.widths] = binadapt(X.data, X.ppm, F.parms, B.wbin);

% Calculate a PCA model from the binned data.
% --
% This function creates a structure that contains information about the
% principal components of the binned dataset. The algorithm used is based
% on a singular value decomposition of B that allows for computationally
% expedient generation of cross-validated datasets, Rsq and Qsq.
pcaModel = pca(B.data);

% Add class information to the PCA model just generated.
% --
% This function adds the Y-matrix into the PCA model structure, allowing
% class-based functions (like scoresplot) to access the information and
% display more useful plots without requiring extra parameters.
pcaModel = addclasses(pcaModel, cls.Y);
pcaModel = addlabels(pcaModel, cls.labels);

