function [preictal_aed_load,null_aed_loads] = get_avg_preictal_levels(ptIDs, all_meds, all_dose_curves, all_tHr, ieeg_offset,max_dur,emu_dur)

sz_cluster =0;

%seizures within a pt
all_aed_mean =zeros(1,length(ptIDs));
preictal_aed_load =cell(1,length(ptIDs));
null_aed_loads = cell(1,length(ptIDs));

for ipt =1:length(ptIDs)
    ptID = ['HUP' num2str(ptIDs(ipt))];
    [med_names,~,~] = parse_MAR(ptID,all_meds);
    
    % get the total AED dose over time
    drugs =zeros(length(med_names),max_dur*60); %450 hours of EMU stay in minutes
    for i =1:length(med_names)
        drug=all_dose_curves{ipt}{i};
        drug=drug./nanmax(drug); %normalize each drug curve
        if ~isempty(drug)
            dStart = round(all_tHr{ipt}{i}(1)*60)-1;
            drugs(i,dStart+1:dStart+length(drug))=drug;
        end
    end
    drug_sum=nansum(drugs,1);
    drug_sum=drug_sum./length(med_names); %normalize for number of drugs
    
    % cut off drug curve to only be length of emu stay
    % time = length(all_spike_rate{ipt})*10 + (sum(offsets)./60); %(estimated from spike data and eeg offsets from drug data
    
    drug_sum(drug_sum==0) =NaN; %not include all zeros in average, and in histogram
    
    %get seizure times
    offsets = ieeg_offset{2,ipt};
    ieeg_offset_datasets = ieeg_offset{1,ipt};
    [seizure_times] = convert_sz_to_emu_time(offsets,ieeg_offset_datasets,ptID);
    
    % convert seizure times to indices, so to minutes
    seizure_inds = round(seizure_times(:,1) *60);
    %remove seizure inds that are greater than emu recording
    time =emu_dur(ipt)*60;%number of minutes of emu stay
    seizure_inds(seizure_inds>time)=[];
    
    % here we want to prune the used seizure inds to pick the first seizure from a cluster
    if sz_cluster ==1 && ~isempty(seizure_inds)
        cutoff = 60*1;
        dSz_times =[cutoff+1; diff(seizure_inds)];
        leadSz_inds = dSz_times > cutoff; %seizures more than one hour apart
        seizure_inds =seizure_inds(leadSz_inds);
    end
    
    % med curves are sampled at one point per minute
    pt_preictal_1hr = zeros(1,length(seizure_inds));
    for i=1:length(seizure_inds)
        pt_preictal_1hr(i) = nanmean(drug_sum(seizure_inds(i)-60:seizure_inds(i)));        
    end
    
    % create random null distribution
    drug_sum(isnan(drug_sum))=[];
    drug_inds = 61:length(drug_sum);
    iters =10;
    pt_null = zeros(iters,length(seizure_inds));
    for i = 1:iters
        rand = randperm(length(drug_inds)); % start an hour in at least
        null_inds = drug_inds(rand(1:length(seizure_inds)));
        
        for j=1:length(null_inds)
            pt_null(i,j)= nanmean(drug_sum(null_inds(j)-60:null_inds(j)));
        end
        
    end
    preictal_aed_load{ipt} = pt_preictal_1hr;
    null_aed_loads{ipt} = pt_null(:);
end
end

%% simple statistical test:

% one-sample t-test

%p0 = [0, .1]; % what would the null standard deviation be?
%p1 = -.2;
%nout = sampsizepwr('t',p0,p1,.9);
%[~,ttest_p,~,stats] = ttest(preictal_AEDs_1hr_mean);


% biniomial distrubution test
% successes = raw_preictal_AEDs_1hr < all_aed_mean;
% binomial_p = binopdf(sum(successes),length(successes),0.5);

