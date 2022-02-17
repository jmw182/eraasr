function data_blanked = apply_blanking(data_tensor, stim_tensor, varargin)
% ERAASR.RNEL.apply_blanking(data_tensor, stim_tensor, ...)
% data_tensor is nTrials x nTime x nChannels tensor 
% (also support nTrials x 1 cell of nTime x nChannels matrices??)
% stim_tensor is nTrials x nPulses matrix
% (also support nTrials x 1 cell of nPulses x 1 vectors??)

% parse inputs
p = inputParser();
p.addParameter('blank_samples', 48, @isscalar) % how many samples to blank
p.addParameter('method', 'linear_ramp', @ischar) % method can be 'sample_and_hold' or 'linear_ramp'
p.addParameter('blank_offset', 2, @isscalar) % offset blank start from the stim times in stim_tensor. Can be positive or negative number of samples.
p.addParameter('showProgress', false, @islogical)
p.parse(varargin{:});
blank_samples = p.Results.blank_samples;
method = p.Results.method;
blank_offset = p.Results.blank_offset;
showProgress = p.Results.showProgress;

% intialize blanked tensor
data_blanked = data_tensor;
if blank_samples == 0
    return
end

% apply blank_offset to stim_tensor
stim_tensor = stim_tensor + blank_offset;

% determine number of trials and pulses
[nTrials, nPulses] = size(stim_tensor);
nChans = size(data_tensor, 3);

% loop over trials and pulses and apply blanking
if showProgress
    prog = ERAASR.Utils.ProgressBar(nTrials*nPulses, 'Applying blanking to data');
end

iCount = 0;
for iTrial = 1:nTrials
    for iPulse = 1:nPulses
        iCount = iCount+1;
        if showProgress
            prog.update(iCount);
        end
        
        % data_tensor is nTrials x nTime x nChannels tensor
        % stim_tensor is nTrials x nPulses matrix
        p0 = stim_tensor(iTrial, iPulse); % first sample to blank
        p_1 = p0 - 1; % last sample before blanking
        pEnd = p0 + blank_samples - 1;
        pIdx = p0:pEnd; % sample index to blank for this pulse
        switch method
            case 'sample_and_hold' 
                data_blanked(iTrial, pIdx, :) = repmat(data_tensor(iTrial, p_1, :), 1, blank_samples, 1);
            case 'linear_ramp'
                v0 = squeeze(data_tensor(iTrial, p_1, :));
                vEnd = squeeze(data_tensor(iTrial, pEnd + 1, :));
                for iChan = 1:nChans
                    data_blanked(iTrial, pIdx, iChan) = linspace(v0(iChan), vEnd(iChan), blank_samples);
                end
            otherwise
                error('blank method not implemented.');
        end

    end
end
if showProgress
    prog.finish();
end

end