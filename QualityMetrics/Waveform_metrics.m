


function [amp_median, Amp_Coeff_of_Disp, amp_mean, amp_std, amp_CV, cluster_amps] = Waveform_metrics(Directory, spikes, clusters,recording_duration,sampling_frequency,size)

% -----------------------------
% function [amp_median, Amp_Coeff_of_Disp, amp_mean, amp_std, amp_CV, cluster_amps] = Waveform_metrics(Directory, spikes, clusters,recording_duration,sampling_frequency)
% -----------------------------
% 
% Computes waveform amplitude metrics based on the absolute voltage of the
% main channel across time. This method considers drift because it uses a
% main channel update approach.
%
% This approach is computationally expensive, so a portion of the code
% relies in GPU computations.
%
% Part of these metrics were inspired in the SpikeInterface Documentation:
% https://spikeinterface.readthedocs.io/en/latest/modules/qualitymetrics.html
%
% -----------------------------
% REQUIRED TOOLBOXES
% - Parallel Computing Toolbox
% -----------------------------
% INPUTS:
%   
%   Directory - Kilosort/phy directory. Location of the raw *ap*.bin file.
%
%   spikes - spike times (in seconds) as a 1 dimensional vector
%            such as [0.1, 0.3, 0.7, 1, 1.6]. Same lenght as clusters.
%
%   clusters - spike clusters (in clusters ID) as a 1 dimensional vector
%              such as [12, 14, 12, 1, 2]. Same lenght as spikes
%              with the previous example, the spike times of cluster 12 
%              would be [0.1, 0.7].
%   
%   recording_duration - recording duration in seconds (scalar).
%
%   sampling_frequency - in Hz, typically 30000 (scalar).
%
%   size - number of waveforms considered for analysis:
%          "All" - All waveforms would be considered - this option is time 
%          consuming and computationally demanding.
%          "default" - Uses a maximum of 10000 spikes - (recommended for
%          quality metrics)
%
% OUTPUTS:
%   
%   amp_median - median of the spike amplitude distribution.
%   
%   Amp_Coeff_of_Disp - Quartile Coefficient of dispersion, computed as:
%                       (Q3 - Q1)/(Q3 + Q1).
%
%   amp_mean - mean of the spike amplitude distribution.
%   
%   amp_std - standard deviation of the spike amplitude distribution.
%
%   amp_CV - coefficient of variation of the spike amplitude distribution.
% 
%   cluster_amps - Cell array with the spike amplitudes of every cluster.
% 
%------------------------------------------
% Copyright (C) 2024 by Juan Pimiento
%------------------------------------------


%Initial parameters

good_clusters = unique(clusters);
bins = ceil(recording_duration/60);
window = [-40, 41];

% Initial parameters of getWaveforms function

gwfparams.dataDir = Directory;    % KiloSort/Phy output folder
apD = dir(fullfile(Directory, '*ap*.bin')); % AP band file from spikeGLX specifically
gwfparams.fileName = apD.name;         % .dat file containing the raw 
gwfparams.dataType = 'int16';            % Data type of .dat file (this should be BP filtered)
gwfparams.nCh = 385;                      % Number of channels that were streamed to disk in .dat file
gwfparams.wfWin = window;              % Number of samples before and after spiketime to include in waveform


%Creating empty variables

amp_median = nan(length(good_clusters),1);
Amp_Coeff_of_Disp = nan(length(good_clusters),1); % Quartile coefficient of dispersion
amp_mean = nan(length(good_clusters),1);
amp_std = nan(length(good_clusters),1);
amp_CV = nan(length(good_clusters),1);
cluster_amps = cell(length(good_clusters),2);

tic %start stopwatch (time consumption)

for z = 1:length(good_clusters)
   
    %Changing parameters of getWaveforms function
    gwfparams.spikeTimes = ceil(spikes(clusters == good_clusters(z)).*sampling_frequency); % Vector of cluster spike times (in samples) same length as .spikeClusters
    gwfparams.spikeClusters = clusters(clusters == good_clusters(z));
    
    if size == "default" & length(gwfparams.spikeTimes) > 10000
        gwfparams.nWf = 10000;
    else
    gwfparams.nWf = length(gwfparams.spikeTimes); % Number of waveforms per unit to pull out
    end

    wf = getWaveForms(gwfparams);
    
    mean_correction = mean(squeeze(wf.waveFormsMean(:,:,1:10)),2);
    wf.waveFormsMean = gpuArray(squeeze(wf.waveFormsMean) - mean_correction);
    array_info = whos("wf");
    if array_info.bytes < 8000000000
        wf.waveForms = gpuArray(squeeze(wf.waveForms));
    else
        wf.waveForms = squeeze(wf.waveForms);
    end


    for w = 1:gwfparams.nWf
        wf.waveForms(w,:,:) = squeeze(wf.waveForms(w,:,:)) - mean_correction;     
    end
    
    max_channel = nan(bins,1);
    binned_spikes = nan(bins,1);
    
    for x = 1:bins
        sl_window = [(x-1)*60,x*60]; % 60 seconds in a minute
        idx = ((wf.spikeTimeKeeps./sampling_frequency) > sl_window(1)) & ((wf.spikeTimeKeeps./sampling_frequency) <sl_window(2));
        window_mean = squeeze(mean((wf.waveForms(idx,:,:)),1));
        [~,max_channel(x)]= max(max(abs(window_mean),[],2));
        binned_spikes(x) = sum(idx);
    end
    
    max_channel(binned_spikes==0) = []; 
    binned_spikes(binned_spikes==0) = [];
    [binned_spikes,ia,~] = unique([1; cumsum(binned_spikes)]);
    ia = ia(2:end)-1;
    max_channel = max_channel(ia);
    
    central_point = ceil(length(window(1):window(2))./2); % To avoid contamination of surrounding spikes 
    central_window = central_point-20:central_point+20;
    amplitudes = nan(gwfparams.nWf,1);
    for q = 2:length(binned_spikes)
            blocks = abs(squeeze(wf.waveForms(binned_spikes(q-1):binned_spikes(q),max_channel(q-1),central_window)));
            amplitudes(binned_spikes(q-1):binned_spikes(q)) = max(blocks,[],2);
    end
    cluster_amps{z,1} = amplitudes;
    cluster_amps{z,2} = wf.spikeTimeKeeps./sampling_frequency;
    amp_median(z) = median(amplitudes);
    Amp_Coeff_of_Disp(z) = (prctile(amplitudes,75) - prctile(amplitudes,25))/(prctile(amplitudes,75) + prctile(amplitudes,25)); % Quartile coefficient of dispersion
    amp_mean(z) = mean(amplitudes);
    amp_std(z) = std(amplitudes);
    amp_CV(z) = amp_std(z)/amp_mean(z);

    disp(['Completed ' int2str(z) ' units of ' int2str(length(good_clusters)) '.']);
end

toc % stop stopwatch (time consumption)


end
