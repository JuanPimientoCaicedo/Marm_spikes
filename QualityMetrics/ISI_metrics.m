


function [ISI_fraction, ISI_falsePositiveRate] = ISI_metrics(spikes, clusters,ISI_interval)

% ---------------------------------
% [ISI_fraction, ISI_falsePositiveRate] = ISI_metrics(spikes, clusters,ISI_interval)
% ---------------------------------
% computes ISI violation metrics for individual clusters (putative neurons)
% for more information about ISI violations:
% https://github.com/cortex-lab/sortingQuality/tree/master
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
%   ISI_interval - Interval for an ISI to be considered an ISI violation,
%   there is no clear agreement about this value. Generally, the usual
%   range used is 1.5 - 2 ms.
% 
% OUTPUTS:
%   ISI_fraction - Typical metric of contamination: proportion of spikes
%   that are next to other spikes with an ISI of less than 
%   ISI_interval (function input).
% 
%   ISI_falsePositiveRate - go to ISIViolations function for more
%   information: https://github.com/cortex-lab/sortingQuality/blob/master/core/ISIViolations.m
%
% Both outputs are 1 dimensional vectors of lenght = unique(clusters), ordered in an
%        ascending fashion.
%        with the previous example, they would be computed in this order:
%        [1, 2, 12, 14]. each of these values is a cluster ID.
%
%------------------------------------------
% Copyright (C) 2024 by Juan Pimiento
%------------------------------------------


good_clusters = unique(clusters);
ISI_fraction = nan(length(good_clusters),1);
for x = 1:length(good_clusters)
    Unit_ISI = diff(spikes(clusters == good_clusters(x)));
    ISI_percentage = sum(Unit_ISI < ISI_interval)./length((Unit_ISI+1)); % added 1 to count all spikes
    ISI_fraction(x) = ISI_percentage;
end
ISI_falsePositiveRate = nan(length(good_clusters),1); % This value generally goes from 0 to 1. 0.5 is an acceptable value.

for y = 1:length(good_clusters)
    Unit_violations = ISIViolations(spikes(clusters == good_clusters(y)),0,ISI_interval);
    ISI_falsePositiveRate(y) = real(Unit_violations);
end
