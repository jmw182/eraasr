%% Load RNEL test dataset

% loads variable data: ~40 trials x samples at 30kHz  x 128 channels
% stim times in stimTensor
% ~20 biphasic pulses at different frequencies per set (20, 50, 100, 300 Hz)

% input set number
%set = 3; % 3, 4, 5, 6
for set = 3 %[3 4 5 6]
iter = 1; % for multiple iterations of output files

% input data directory and filename
data_dir = 'D:\StimRecordTest\2022_02_15_data_for_eraasr';
fname = sprintf('Set%d_tensor_ERAASR.mat', set);

% output directory and filename
data_out_dir = 'D:\git\eraasr\data_out';
data_out_fname = sprintf('Set%d_Iter%d_eraasr_output.mat', set, iter);

% load input data
load(fullfile(data_dir, fname)); % loads set_tensor, stim_tensor

%% Setup ERAASR Parameters

opts = ERAASR.Parameters();
opts.Fs = 30000; % samples per second (Hz)
Fms = opts.Fs / 1000; % multiply to convert ms to samples

% thresholding to detect stim
opts.thresholdChannel = 16; % threshold which channel
opts.thresholdValue = -2000; % threshold channel at this value to find approximate start 
opts.thresholdHPCornerHz = 250; % High pass filter at this corner frequency before thresholding

% For alignment
opts.alignChannel = 1;
opts.alignUpsampleBy = 10; % supersample by this ratio before alignment
opts.alignWindowPre = Fms * 0.5;
opts.alignWindowDuration = Fms * 12;%double(range(stim_tensor,'all')) + (Fms*10); %Fms * 12;

% For Cleaning
opts.extractWindowPre = Fms * 20; %double(stim_tensor(1,1) - Fms*20); %Fms * 20; % number of samples to extract before the threshold crossing, should be sufficient to allow HP filtering 
opts.extractWindowDuration = double(range(stim_tensor,'all')) + (Fms*200); %Fms * 100; % total number of samples to extract, which includes the pre, during stim, and post stim windows
opts.cleanStartSamplesPreThreshold = Fms * 0.5; % number of samples before threshold crossing to include within the first pulse
        
opts.cleanHPCornerHz = 10; % light high pass filtering at the start of cleaning
opts.cleanHPOrder = 4; % high pass filter order 
opts.cleanUpsampleBy = 1; % upsample by this ratio during cleaning
opts.samplesPerPulse = Fms * 0.7; % duration of stim pulse
opts.nPulses = size(stim_tensor,2); % number of pulses per trial

% number of principal components for cleaning
opts.nPC_channels = 12;
opts.nPC_trials = 2;
opts.nPC_pulses = 6;

% when reconstructing each, omit this number of *TOTAL*
% channels/trials/pulses, including the one being reconstructed. 
% So 3 means omit me and my immediate neighbors
opts.omit_bandwidth_channels = 3;
opts.omit_bandwidth_trials = 1;
opts.omit_bandwidth_pulses = 1;

opts.alignPulsesOverTrain = true; %false; % do a secondary alignment within each train, in case you think there is pulse to pulse jitter. Works best with upsampling
opts.pcaOnlyOmitted = true; % if true, build PCs only from non-omitted channels/trials/pulses. if false, build PCs from all but set coefficients to zero post hoc

opts.cleanOverChannelsIndividualTrials = false;
opts.cleanOverPulsesIndividualChannels = false;
opts.cleanOverTrialsIndividualChannels = false;

opts.cleanPostStim = true; % clean the post stim window using a single PCR over channels

opts.showFigures = false; % useful for debugging and seeing well how the cleaning works
opts.plotTrials = 1; % which trials to plot in figures, can be vector
opts.plotPulses = 1; % which pulses to plot in figures, can be vector
opts.figurePath = pwd; % folder to save the figures
opts.saveFigures = false; % whether to save the figures
opts.saveFigureCommand = @(filepath) print('-dpng', '-r300', [filepath '.png']); % specify a custom command to save the figure

opts.quiet = false; % if true, don't print anything

%% Do alignment and cleaning procedure

[dataCleaned, extract] = ERAASR.cleanTrials(set_tensor, opts);

%% Plot the cleaned traces on a single trial

figure();
plot(squeeze(dataCleaned(1, :, :)));
box off;

%% Apply blanking before filtering
blank_opts.applyBlanking = true;
if blank_opts.applyBlanking
    blank_opts.method = 'linear_ramp';
    blank_opts.blank_samples = 50;
    dataBlanked = ERAASR.RNEL.apply_blanking(dataCleaned, stim_tensor, 'blank_samples', blank_opts.blank_samples, 'method', blank_opts.method);
else
    dataBlanked = dataCleaned;
end

%% Note before spike extraction
% It would presumably make sense to combine the stimulated trials with any 
% non-stimulated trials to ensure that the broadband signals are treated
% identically from this point forward

%% High pass filter the cleaned data
filt.cornerHP = 750;
filt.order = 1;
filt.filtfilt = false;
dataCleanedHP = ERAASR.highPassFilter(dataBlanked, opts.Fs, 'cornerHz', filt.cornerHP, 'order', filt.order, ...
    'subtractFirstSample', true, 'filtfilt', filt.filtfilt, 'showProgress', true);
        
%% Spike thresholding and waveform extraction

spikes.rmsThresh = -4.5 * ERAASR.computeRMS(dataCleanedHP, 'perTrial', false, 'clip', 60); % clip samples that sink outside +/- 60 uV

waveSamplesPrePost = [10 38];
[spikes.spikeTimes, spikes.waveforms] = ERAASR.extractSpikesCrossingThreshold(dataCleanedHP, spikes.rmsThresh, ...
    'mode', 'largestFirst', 'waveformSamplesPrePost', waveSamplesPrePost, 'lockoutPrePost', [9 30]);

%% Plot the mean spike waveforms
% Note that this will include some spikes outside of the stimulation period as is

nChannels = size(spikes.waveforms, 2);
nSamples = sum(waveSamplesPrePost);
spikes.waveformMeans = nan(nSamples, nChannels);
for iC = 1:nChannels
    spikes.waveformMeans(:, iC) = mean(cat(1, spikes.waveforms{:, iC}), 1);
end

figure();
plot(spikes.waveformMeans);
box off;

%% save data
out_path = fullfile(data_out_dir, data_out_fname);
out_vars = {'dataBlanked', 'dataCleaned', 'dataCleanedHP', 'spikes', 'opts', 'blank_opts', 'filt', 'extract', 'set', 'iter'};
save(out_path, out_vars{:}, '-v7.3');

end % set loop