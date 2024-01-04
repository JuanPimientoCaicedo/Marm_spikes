

function Pr = Presence_ratio(spikes, clusters, recording_duration, spacing)
% ---------------------------------
% Pr = Presence_ratio(spikes, clusters, recording_duration, spacing)
% ---------------------------------
% computes presence ratio for individual clusters (putative neurons)
% for more information about presence ratio:
% https://spikeinterface.readthedocs.io/en/latest/modules/qualitymetrics/presence_ratio.html
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
%   recording_duration - recording duration in seconds (scalar)
%  
%   spacing - time bins to compute Presence ratio in seconds (scalar)
%             60 seconds is recommended.
% 
% OUTPUTS:
%   Pr - 1 dimensional vector of lenght = unique(clusters), ordered in an
%        ascending fashion.
%        with the previous example, Pr would be computed in this order:
%        [1, 2, 12, 14]. each of these values is a cluster ID.
%------------------------------------------
% Copyright (C) 2024 by Juan Pimiento
%------------------------------------------


good_clusters = unique(clusters);
Pr = NaN(size(good_clusters));
Bins = round(recording_duration/spacing);
time_vector = linspace(0,recording_duration,Bins+1);

for x = 1:length(good_clusters)
    sp = spikes(clusters==good_clusters(x));
    binsCounts = histcounts(sp,time_vector);
    Pr(x) = sum(binsCounts>0)./Bins;
end