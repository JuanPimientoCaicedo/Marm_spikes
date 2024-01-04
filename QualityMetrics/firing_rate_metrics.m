function [firing_rate, firing_rate_mean, firing_rate_std, firing_rate_range, Total_spikes] = firing_rate_metrics(spikes,clusters,recording_duration,binsize)
% ---------------------------------
% [firing_rate, firing_rate_mean, firing_rate_std, firing_rate_range] = firing_rate_metrics(spikes,clusters,recording_duration,binsize)
% ---------------------------------
% computes firing rate metrics for individual clusters (putative neurons)
% for more information about quality metrics:
% https://spikeinterface.readthedocs.io/en/latest/modules/qualitymetrics.html
% ---------------------------------
% INPUTS:
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
%   binsize - binsize (in seconds) to compute firing rate metrics, including firing
%   rate mean, standard deviation and ranges. It is recommended to have a
%   small enough binsize, for example 10 seconds. 
% 
% OUTPUTS:
%   firing_rate - gross firing rate during recording (total
%   spikes/recording duration).
%
%   firing_rate_mean - mean firing rate of binned firing rate distribution.
%
%   firing_rate_std - standard deviation of the firing rate of binned 
%                     firing rate distribution.
%
%   firing_rate_range - 5th and 95th percentiles of binned firing rate 
%                       distribution.
%
%   Total_spikes - Cluster total spikes
%------------------------------------------
% Copyright (C) 2024 by Juan Pimiento
%------------------------------------------

%firing rate and total spikes across the whole recording

good_clusters = unique(clusters);
firing_rate = nan(length(good_clusters),1);
Total_spikes = nan(length(good_clusters),1);
for x = 1:length(good_clusters)
    fr = length(spikes(clusters==good_clusters(x)))./(recording_duration);
    Total_spikes(x) = length(spikes(clusters==good_clusters(x)));
    firing_rate(x) = fr;
end


%firing rate mean, range and std

Bins = round(recording_duration/binsize);
time_vector = linspace(0,recording_duration,Bins+1);
firing_rate_range = nan(length(good_clusters),2);
firing_rate_std = nan(length(good_clusters),1);
firing_rate_mean = nan(length(good_clusters),1);

for y = 1:length(good_clusters)
    fr_bins = (histcounts(spikes(clusters==good_clusters(y)),time_vector))./binsize;
    [firing_rate_std(y), firing_rate_mean(y)] = std(fr_bins);
    firing_rate_range(y,:) = prctile(fr_bins,[0, 95]);
end

end