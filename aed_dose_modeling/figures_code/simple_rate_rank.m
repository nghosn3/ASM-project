function [all_ranks,all_aed_ranks,successes,successes_val,aed_chance,aed_load_data] = simple_rate_rank(aed_loads,sz_aed)

%% Initialize data things
npts = length(aed_loads);
all_ranks = cell(npts,1);
all_aed_ranks = cell(npts,1);

nchance = nan(npts,1);
all = nan(npts,1);
successes = nan(npts,1);
successes_val = nan(npts,1);

aed_chance = nan(npts,2);
aed_load_data=nan(npts,2);

for i = 1:npts
    curr_aed_sz = sz_aed{i};
    curr_aed = aed_loads{i};
    
    aed_load_data(i,:)=[nanmedian(curr_aed_sz) nanmedian(curr_aed)];
    aed_chance(i,:) = [nanmedian(curr_aed_sz),nanmedian(curr_aed)]; %median of seizure aed load, median of all aed loads
    
%     % remove nans
%     nan_things = isnan(curr_aed);
%     curr_aed_sz(nan_things) = [];
%     curr_aed_sz = logical(curr_aed_sz);    
%     curr_aed(nan_things) = [];
    
    % Get rank of SOZ
    % sort the rates in descending order
    [~,I] = sort(curr_aed,'descend');
    ranks = 1:length(curr_aed);
    ranks(I) = ranks;
    sz_ranks = sum((median(curr_aed_sz)) >curr_aed);
    successes_val(i) = nanmean(curr_aed_sz) < nanmean(curr_aed); %comparing mean value not rank
    successes(i) = sz_ranks< length(curr_aed)/2; 
    all_ranks{i} = ranks;
    all_aed_ranks{i} = sz_ranks;

end