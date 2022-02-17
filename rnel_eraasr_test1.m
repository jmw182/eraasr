%% Load RNEL test dataset

% loads variable data: ~40 trials x samples at 30kHz  x 128 channels
% stim times in stimTensor
% ~20 biphasic pulses at different frequencies per set (20, 50, 100, 300 Hz)

data_dir = 'C:\Data\StimRecordCableTest\ERAASR';
set = 3; % 3, 4, 5, 6
fname = sprintf('Set%d_tensor_ERAASR.mat', set);
load(fullfile(data_dir, fname)); % loads set_tensor, stim_tensor

%% Setup ERAASR Parameters

opts = ERAASR.Parameters();
opts.Fs = 30000; % samples per second
Fms = opts.Fs / 1000; % multiply to convert ms to samples

opts.thresholdHPCornerHz = 250;
opts.thresholdChannel = 16;
opts.thresholdValue = -2000;

opts.alignChannel = 1;
opts.alignUpsampleBy = 10;
opts.alignWindowPre = Fms * 0.5;
opts.alignWindowDuration = Fms * 12;%double(range(stim_tensor,'all')) + (Fms*10); %Fms * 12;

% 60 ms stim, align using 20 ms pre start to 110 post
opts.extractWindowPre = Fms * 20; %double(stim_tensor(1,1) - Fms*20); %Fms * 20;
opts.extractWindowDuration = double(range(stim_tensor,'all')) + (Fms*200); %Fms * 100;
opts.cleanStartSamplesPreThreshold = Fms * 0.5;
        
opts.cleanHPCornerHz = 10; % light high pass filtering at the start of cleaning
opts.cleanHPOrder = 4; % high pass filter order 
opts.cleanUpsampleBy = 1; % upsample by this ratio during cleaning
opts.samplesPerPulse = Fms * 0.7; %Fms * 3; % 3 ms pulses
opts.nPulses = size(stim_tensor,2); %20;

opts.nPC_channels = 12;
opts.nPC_trials = 2;
opts.nPC_pulses = 6;

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

%% Do alignment and cleaning procedure

[dataCleaned, extract] = ERAASR.cleanTrials(set_tensor, opts);

%% Plot the cleaned traces on a single trial

figure();
plot(squeeze(dataCleaned(1, :, :)));
box off;

%% Apply blanking before filtering
applyBlanking = true;
if applyBlanking
    method = 'linear_ramp';
    blank_samples = 50;
    dataBlanked = ERAASR.RNEL.apply_blanking(dataCleaned, stim_tensor, 'blank_samples', blank_samples, 'method', method);
else
    dataBlanked = dataCleaned;
end

%% Note before spike extraction
% It would presumably make sense to combine the stimulated trials with any 
% non-stimulated trials to ensure that the broadband signals are treated
% identically from this point forward

%% High pass filter the cleaned data

dataCleanedHP = ERAASR.highPassFilter(dataBlanked, opts.Fs, 'cornerHz', 750, 'order', 1, ...
    'subtractFirstSample', true, 'filtfilt', false, 'showProgress', true);
        
%% Spike thresholding and waveform extraction

rmsThresh = -4.5 * ERAASR.computeRMS(dataCleanedHP, 'perTrial', false, 'clip', 60); % clip samples that sink outside +/- 60 uV

waveSamplesPrePost = [10 38];
[spikeTimes, waveforms] = ERAASR.extractSpikesCrossingThreshold(dataCleanedHP, rmsThresh, ...
    'mode', 'largestFirst', 'waveformSamplesPrePost', waveSamplesPrePost, 'lockoutPrePost', [9 30]);

%% Plot the mean spike waveforms
% Note that this will include some spikes outside of the stimulation period as is

nChannels = size(waveforms, 2);
nSamples = sum(waveSamplesPrePost);
waveformMeans = nan(nSamples, nChannels);
for iC = 1:nChannels
    waveformMeans(:, iC) = mean(cat(1, waveforms{:, iC}), 1);
end

figure();
plot(waveformMeans);
box off;