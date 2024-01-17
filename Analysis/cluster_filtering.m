

function [Quality_Metrics_Table,spikes,clusters] = cluster_filtering(Quality_Metrics_Table,spikes,clusters)

%------------------------------------------
% [Quality_Metrics_Table,spikes,clusters] = cluster_filtering(Quality_Metrics_Table,spikes,clusters)
%------------------------------------------
% filters clusters that do not fulfill Quality metrics criteria
% *** not a universal function  ***
% *** it depends on the Quality_Metrics_Table configuration*** 
% *** look at the Run_QualityMetrics.mlx file in this toolbox ***
% 
%------------------------------------------
% INPUTS:
%   Quality_Metrics_Table - Quality metrics table organized in the way that
%                           Run_QualityMetrics.mlx file does it.
%
%   spikes - spike times (in seconds) as a 1 dimensional vector
%            such as [0.1, 0.3, 0.7, 1, 1.6]. Same lenght as clusters.
%
%   clusters - spike clusters (in clusters ID) as a 1 dimensional vector
%              such as [12, 14, 12, 1, 2]. Same lenght as spikes
%              with the previous example, the spike times of cluster 12 
%              would be [0.1, 0.7].
%
% OUTPUTS:
%
%   *** This outputs are exactly the same inputs but with "low quality units"
%       filtered ***
%
%   Quality_Metrics_Table - Quality metrics table organized in the way that
%                           Run_QualityMetrics.mlx file does it.
%
%   spikes - spike times (in seconds) as a 1 dimensional vector
%            such as [0.1, 0.3, 0.7, 1, 1.6]. Same lenght as clusters.
%
%   clusters - spike clusters (in clusters ID) as a 1 dimensional vector
%              such as [12, 14, 12, 1, 2]. Same lenght as spikes
%              with the previous example, the spike times of cluster 12 
%              would be [0.1, 0.7].
%
%------------------------------------------
% Written by Juan Pimiento (2024)
%------------------------------------------


%% Quality metrics criteria
% You can modify them depending on what you are interested in, this cutoffs
% are not universal
table_idx = Quality_Metrics_Table.firing_rate < 0.1 | ...
    Quality_Metrics_Table.Presence_ratio < 0.9 | ...
    Quality_Metrics_Table.ISI_fraction > 0.01 | ...
    Quality_Metrics_Table.ISI_fp > 0.5 | ... 
    Quality_Metrics_Table.amp_median < 20;

% Eliminate low quality units
Quality_Metrics_Table(table_idx,:) = [];
good_clusters = Quality_Metrics_Table.ClusterID;

% Clean clusters and spikes vectors from low quality units
spikes = spikes(ismember(clusters,good_clusters)); % first spikes vector!!
clusters = clusters(ismember(clusters,good_clusters));

end

