


function [results] = clean_spikes(Dir,method,samples)
% ------------------------------
% [results] = clean_spikes(Dir,method,samples)
% ------------------------------
% This function eliminates one of a pair of spikes with an artificial ISI violation 
% created by the template matching algorithm of kilosort3. This function is only 
% applied to manually curated spikes in phy, where clusters are classified as 
% follows: 1=MUA, 2=Good, 3=Unsorted.
% 
% The output of this function is automatically saved in a folder called clean_spikes 
% in the Kilosort data directory. With the same file convention as kilosort output.
% 
% For more information about Kilosort3, click: <https://github.com/MouseLand/Kilosort 
% MouseLand/Kilosort: Fast spike sorting with drift correction for up to a thousand 
% channels (github.com)>
% ------------------------------
% INPUTS:
% 
% Dir - The directory where the spike_times.npy and spike_clusters.npy are located.
% 
% method - you can use one of three options:
%          "keep_first" - eliminates the second spike
%          "keep_last" - eliminates the first spike
%          "random" - eliminates just one of the spikes in a random manner
% 
% samples - The number of samples that the function will look for ISI violations 
% - the recommended value is 3 (given that sampling rate is 30000Hz, 3 samples 
% would be 100 microseconds)
% 
% OUTPUTS:
% 
% Results - a table with the cluster identity in one column and the eliminated 
% spikes in the other.
% ------------------------------ 
% Required toolboxes
% This section requires:
% 
% * spikes Toolbox (<https://github.com/cortex-lab/spikes/ cortex-lab/spikes: 
% cortex lab code for electrophysiology (github.com)>)
% * npy-matlab Toolbox (<https://github.com/kwikteam/npy-matlab kwikteam/npy-matlab: 
% Experimental code to read/write NumPy .NPY files in MATLAB (github.com)>)
%  ------------------------------
% Eliminate duplicated spikes
% Kilosort often Double counts some spikes, as discussed in the following issues: 
% <https://github.com/MouseLand/Kilosort/issues/29 Double-counted spikes · Issue 
% #29 · MouseLand/Kilosort (github.com)> <https://github.com/MouseLand/Kilosort/issues/415 
% High rate of zero-lag autocorrelogram (ACG) peak · Issue #415 · MouseLand/Kilosort 
% (github.com)>.
% 
% One way of dealing with this problem is eliminating one of the two double 
% counted spikes, approach previously used by others: <https://spikeinterface.readthedocs.io/en/latest/api.html#spikeinterface.curation.remove_duplicated_spikes 
% API — SpikeInterface documentation>
%------------------------------------------
% Copyright (C) 2024 by Juan Pimiento
%------------------------------------------

% Upload data
sp = loadKSdir(Dir);
clean_st = [];
clean_clu = [];
good_clusters = sp.cids(sp.cgs == 2)'; % Phy cluster labels (1=MUA, 2=Good, 3=Unsorted)
eliminated_spikes = [];

% Eliminate double counted spikes

for x = 1:length(good_clusters)
    spike_times = sp.st(sp.clu == good_clusters(x));
    spike_cluster = sp.clu(sp.clu == good_clusters(x));
    ISI = diff(spike_times);
    ISI_violations = find(ISI < samples./sp.sample_rate);
    if isempty(ISI_violations) == 1
    elseif method == "keep_first"
    spike_times(ISI_violations+1) = [];
    spike_cluster(ISI_violations+1) = [];
    elseif method == "keep_last"
    spike_times(ISI_violations) = [];
    spike_cluster(ISI_violations) = [];
    elseif method == "random"
    rand_vector = randi([0,1],lenght(ISI_violations));
    spike_times(ISI_violations+rand_vector) = [];
    spike_cluster(ISI_violations+rand_vector) = [];
    end
clean_st = [clean_st;spike_times];
clean_clu = [clean_clu;spike_cluster];
eliminated_spikes = [eliminated_spikes;length(sp.st(sp.clu == good_clusters(x)))-length(spike_times)];
end

% order new vectors 
% Note that these new vectors are only 'good' clusters. MUA clusters are
% excluded to make vectors smaller

[clean_st,idx] = sort(clean_st,1,"ascend");
clean_clu = clean_clu(idx);
results = table(good_clusters,eliminated_spikes);
mkdir(fullfile(Dir,'Clean_Spikes'));
writeNPY(clean_st, fullfile(Dir,'Clean_Spikes','Clean_spike_times.npy'));
writeNPY(clean_clu, fullfile(Dir,'Clean_Spikes','Clean_spike_clusters.npy'));
disp(['output saved in' fullfile(Dir,'Clean_Spikes') '.']);