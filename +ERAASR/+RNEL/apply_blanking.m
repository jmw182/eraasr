function data_blanked = apply_blanking(data_tensor, stim_tensor, varargin)
% ERAASR.RNEL.apply_blanking(data_tensor, stim_tensor, ...)
% data_tensor is either:
%  - nTrials x 1 cell of time x channels matrices
%  - nTrials x nTime x channels tensor
% stim_tensor is either
%  - nTrials x 1 cell of nPulses x 1 vectors
%  - nTrials x nPulses matrix
%
% Optional inputs include:
% blank_samples (how many samples to blank, default: 48)
% method ('linear_ramp' or 'sample_and_hold', default: 'linear_ramp')
%       blank_offset (how many samples to offset the stim_tensor times by to
%       initiate blanking. Can be positive (delays) or negative (precedes).
%       Default: 2)
% showProgress (default: false)

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


f_dataCell = iscell(data_tensor);
f_stimCell = iscell(stim_tensor);
%f_cell = f_dataCell || f_stimCell;
% if either input is a cell, make sure both are
% if f_cell
%     if f_dataCell && ~f_stimCell % convert stimTensor to cell
%         stim_tensor = num2cell(stim_tensor, 2); % split trials into cell of cells
%         stim_tensor = cellfun(@(x) x{:}', stim_tensor, 'UniformOutput', false); % extract innermost cells and transpose to match documentation
%         f_stimCell = true;
%     elseif f_stimCell && ~f_dataCell % convert dataTensor to cell
%         data_tensor = num2cell(data_tensor, [2 3]); % split trials into cell of cells
%         data_tensor = cellfun(@(x) x{:}', data_tensor, 'UniformOutput', false); % extract innermost cells and transpose to match documentation
%         data_blanked = data_tensor;
%         f_dataCell = true;
%     end
% end

% apply blank_offset to stim_tensor
if f_stimCell
    stim_tensor = cellfun(@(x) x + blank_offset, stim_tensor, 'UniformOutput', false);
else
    stim_tensor = stim_tensor + blank_offset;
end

% determine number of trials and pulses
if f_stimCell
    nTrials = length(stim_tensor);
    nPulses = cellfun('length', stim_tensor); % vector of number of pulses per trial
else
    [nTrials, nPulses] = size(stim_tensor);
    nPulses = repmat(nPulses, 1, nTrials); % expand to match format for f_stimCell
end

if f_dataCell
    nChans = cellfun(@(x) size(x,2), data_tensor); % vector of number of channels per trial
else
    nChans = size(data_tensor, 3);
    nChans = repmat(nChans, 1, nTrials); % expand to match format for f_dataCell
end

% loop over trials and pulses and apply blanking
if showProgress
    prog = ERAASR.Utils.ProgressBar(sum(nPulses), 'Applying blanking to data');
    iCount = 0;
end


for iTrial = 1:nTrials
    for iPulse = 1:nPulses(iTrial)

        if showProgress
            iCount = iCount+1;
            prog.update(iCount);
        end

        if f_stimCell % stim_tensor is nTrials x 1 cell of nPulses x 1 vectors
            p0 = stim_tensor{iTrial}(iPulse); % first sample to blank
        else % stim_tensor is nTrials x nPulses matrix
            p0 = stim_tensor(iTrial, iPulse); % first sample to blank
        end
        p_1 = p0 - 1; % last sample before blanking
        pEnd = p0 + blank_samples - 1;
        pIdx = p0:pEnd; % sample index to blank for this pulse

        if f_dataCell
            % data_tensor is nTrials x 1 cell of nTime x nChannels matrices

            switch method
                case 'sample_and_hold'
                    data_blanked{iTrial}(pIdx, :) = repmat(data_tensor{iTrial}(p_1, :), blank_samples, 1);
                case 'linear_ramp'
                    v0 = data_tensor{iTrial}(p_1, :);
                    vEnd = data_tensor{iTrial}(pEnd + 1, :);
                    for iChan = 1:nChans(iTrial)
                        data_blanked{iTrial}(pIdx, iChan) = linspace(v0(iChan), vEnd(iChan), blank_samples);
                    end
                otherwise
                    error('blank method not implemented.');
            end

        else % tensors (not cells)
            % data_tensor is nTrials x nTime x nChannels tensor

            switch method
                case 'sample_and_hold'
                    data_blanked(iTrial, pIdx, :) = repmat(data_tensor(iTrial, p_1, :), 1, blank_samples, 1);
                case 'linear_ramp'
                    v0 = squeeze(data_tensor(iTrial, p_1, :));
                    vEnd = squeeze(data_tensor(iTrial, pEnd + 1, :));
                    for iChan = 1:nChans(iTrial)
                        data_blanked(iTrial, pIdx, iChan) = linspace(v0(iChan), vEnd(iChan), blank_samples);
                    end
                otherwise
                    error('blank method not implemented.');
            end
        end

    end % end pulse loop
end % end trial loop

if showProgress
    prog.finish();
end

end