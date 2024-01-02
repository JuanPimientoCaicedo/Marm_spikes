


function wf = clean_waveforms(Directory, spikes_file, clusters_file, window, number_of_waveforms, sampling_frequency)
% --------------------------------------
% wf = clean_waveforms(Directory, spikes_file, clusters_file, window, number_of_waveforms)
% --------------------------------------
% Extracts individual spike waveforms from the raw datafile, for multiple
% clusters. Returns the waveforms and their means within clusters.
% Also corrects the DC offset for all the channels. 
% 
% For more information about Kilosort3, click: <https://github.com/MouseLand/Kilosort 
% MouseLand/Kilosort: Fast spike sorting with drift correction for up to a thousand 
% channels (github.com)>
% --------------------------------------
% INPUTS:
% Directory - KiloSort/Phy output folder
% 
% spikes_file - Vector of cluster spike times (in seconds) same length as .spikeClusters 
%               *(getWaveforms uses the spikes vector as samples)*
% 
% clusters_file -  Vector of cluster IDs (Phy nomenclature)   same length as 
%                  .spikeTimes
% 
% window - Number of samples before and after spiketime to include in
%          waveform.
%          It is necessary to get a large window before the event to get a 
%          good DC correction - recommended [-40, 41] 
%          *the DC correction is done with the first 10 samples 
%          of the window - That section must be spike free!*
% 
% number_of_waveforms - Number of waveforms per unit to pull out
% 
% sampling_frequency - sampling frequency of the spike_times vector,
%                      typically 30000 Hz.
% OUTPUTS:
%
% wf - the same ouptup than getWaveforms function
%      for more information: https://github.com/cortex-lab/spikes/blob/master/analysis/getWaveForms.m
%  --------------------------------------
% Required toolboxes
%
% * spikes Toolbox (<https://github.com/cortex-lab/spikes/ cortex-lab/spikes: 
% cortex lab code for electrophysiology (github.com)>)
% * npy-matlab Toolbox (<https://github.com/kwikteam/npy-matlab kwikteam/npy-matlab: 
% Experimental code to read/write NumPy .NPY files in MATLAB (github.com)>)
%  --------------------------------------
% getWaveForms function was corrected for indexing bug.
% 
% It is recommended to add these two lines to the function after line 49:
%   curSpikeTimes((curSpikeTimes + gwfparams.wfWin(1)) < 1) = [];
%   curSpikeTimes((curSpikeTimes + gwfparams.wfWin(end)) > size(mmf.Data.x,2)) = [];
%------------------------------------------
% Copyright (C) 2024 by Juan Pimiento
%------------------------------------------


% Get waveforms
good_clusters = unique(clusters_file);
gwfparams.dataDir = Directory;    % KiloSort/Phy output folder
apD = dir(fullfile(Directory, '*ap*.bin')); % AP band file from spikeGLX specifically
gwfparams.fileName = apD.name;         % .dat file containing the raw 
gwfparams.dataType = 'int16';            % Data type of .dat file (this should be BP filtered)
gwfparams.nCh = 385;                      % Number of channels that were streamed to disk in .dat file
gwfparams.wfWin = window;              % Number of samples before and after spiketime to include in waveform
gwfparams.nWf = number_of_waveforms;                    % Number of waveforms per unit to pull out
gwfparams.spikeTimes = ceil(spikes_file.*sampling_frequency); % Vector of cluster spike times (in samples) same length as .spikeClusters
gwfparams.spikeClusters = clusters_file;

wf = getWaveForms(gwfparams);

% Correct DC offset
for x = 1:length(good_clusters)
    mean_correction = mean(squeeze(wf.waveFormsMean(x,:,1:10)),2);
    wf.waveFormsMean(x,:,:) = squeeze(wf.waveFormsMean(x,:,:)) - mean_correction;
    for w = 1:gwfparams.nWf
        wf.waveForms(x,w,:,:) = squeeze(wf.waveForms(x,w,:,:)) - mean_correction;     
    end
    disp(['Completed DC correction from ' int2str(x) ' units of ' int2str(length(good_clusters)) '.']);
end

% Plot example
clu_of_int = randi(length(good_clusters)); % cluster of interest 
figure; 
imagesc(squeeze(wf.waveFormsMean(clu_of_int,:,:)))
set(gca, 'YDir', 'normal'); xlabel('time (samples)'); ylabel('channel number'); 
colormap(colormap_BlueWhiteRed); caxis([-1 1]*max(abs(caxis()))/2); box off;
title(['Voltage deflections of Cluster ',num2str(good_clusters(clu_of_int))])
[y, idx] = max(max(abs(squeeze(wf.waveFormsMean(clu_of_int,:,:))),[],2));
x_axis = gwfparams.wfWin(1):gwfparams.wfWin(2);
figure
plot(x_axis, squeeze(wf.waveForms(clu_of_int,:,idx,:)),"Color",[0.5 0.5 0.5])
hold on
plot(x_axis, squeeze(wf.waveFormsMean(clu_of_int,idx,:)),"LineWidth",1);
xlim([gwfparams.wfWin(1) gwfparams.wfWin(2)]);
title(['Waveform average of Cluster ',num2str(good_clusters(clu_of_int))]);
xlabel('time (samples)'); ylabel('voltage (\muV)');
hold off