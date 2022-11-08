%% run bootstrapped regression for all patients and all drugs
close all;clear;

% Get which patients
cohort_info = readtable('HUP_implant_dates.xlsx');
ptIDs = cohort_info.ptID;
weights = cohort_info.weight_kg;
baseline_sz_freqs = readtable('no_phi_baseline_sz_freq.xlsx');

% use only tapered patients? n=70
tapered = 0;

% get AED curves
load('MAR_032122.mat')

if tapered
    [tapered_pts,all_taper_info] = get_tapered_patients(ptIDs,all_meds);
    ptIDs = ptIDs(logical(tapered_pts));
    weights = weights(logical(tapered_pts));
end

% using the AED BPL model
[all_dose_curves,all_tHr,ptIDs,all_med_names,ieeg_offset,max_dur,emu_dur] = get_aed_curve_kg(ptIDs,weights,0,[]); % use min dose to normalize dose

% using only interpolation dosage
%[all_dose_curves,all_tHr,all_med_names,ieeg_offset,emu_dur] = get_dose_schedules(ptIDs);

%%
exclude_ativan=0;
rm_clusters =0;

medications = vertcat(all_med_names{:});
all_medications = unique(medications);
if exclude_ativan
    ativan_ind = contains(all_medications,'lorazepam');
    all_medications(ativan_ind)=[];
end

% make matrix with all variables (levels for each drug)
num_feats=length(all_medications)+4; % (baseline sz freq, time since last seizure, ptID, seizure bin)
sz=num_feats;
features = cell(num_feats,length(ptIDs));
winLen =60; %minutes

for i=1:length(ptIDs) % exclude HUP136
    ptID = ['HUP' num2str(ptIDs(i))];
    % for grabbing seizure times
    offsets = ieeg_offset{2,i};
    ieeg_offset_datasets = ieeg_offset{1,i};
    
    med_names = all_med_names{i};
    pt_curves = all_dose_curves{i};
    tHr=all_tHr{i};
    
    % exclude ativan in drug sum
    if exclude_ativan
        ativan_ind = contains(med_names,'lorazepam'); %#ok<*UNRCH>
        med_names(ativan_ind)=[];
        pt_curves(ativan_ind)=[];
        tHr(ativan_ind)=[];
    end
    
    % align drugs on same time scale
    end_times = cellfun(@(x)x(end),tHr);
    curve_lens = cellfun(@length,pt_curves);
    max_time = max(end_times);
    drugs =zeros(length(med_names),ceil(max_time*60))+1e-10;
    for n =1:length(med_names)
        drug=pt_curves{n};
        drug=drug./nanmax(drug); %normalize drug curve
        if ~isempty(drug)
            dStart = round(tHr{n}(1)*60)-1;
            drugs(n,dStart+1:dStart+length(drug))=drug;
        end
    end
    
    % get binned drug values
    nbins=ceil(length(drugs)./winLen);
    drug_wins = nan(1,nbins);
    
    for r = 1:length(med_names)
        ind =1;
        for c = 1:winLen:length(drugs)-winLen
            drug_wins(r,ind)=nanmean(drugs(r,c:c+winLen)); %creates leading zeros
            ind=ind+1;
        end
    end
    
    for n=1:length(med_names)
        ind = contains(all_medications,med_names{n});
        features(ind,i)={drug_wins(n,:)};
    end
    
    empty_inds =cellfun(@isempty,features(:,i));
    features(empty_inds,i)={zeros(1,nbins)+1e-10}; %adding small number to make columns full rank
    features(end-1,i)={i*ones(1,nbins)}; %ptID
    
    [seizure_times] = convert_sz_to_emu_time(offsets,ieeg_offset_datasets,ptID); % seizure times in hours of EMU stay
    
    % adjust seizure times to winLen minute bins
    seizure_times = seizure_times*60./winLen;
    seizure_times(seizure_times(:,1)>(max_time*60./winLen),:)=[]; %this is for you, HUP175!
    
    % shift seizure bins left so 1 = seizure in next time bin
    sz_bins = zeros(1,nbins+1);
    sz_inds = round(seizure_times(:,1));
    sz_bins(sz_inds) =1;
    sz_bins(1)=[];
    
    features(sz,i)={sz_bins};
    
    % add time since last seizure feature
    t_last_seizure = nan(1,nbins);
    seizure_inds = unique(round(seizure_times(:,1)));
    start = seizure_inds(1);
    ind=1;
    for j=start:seizure_inds(end)
        if ind<length(seizure_inds)
            if j>=seizure_inds(ind) && j < seizure_inds(ind+1)
                t_last_seizure(j)=j-seizure_inds(ind);
                
            else
                ind=ind+1;
                t_last_seizure(j)=j-seizure_inds(ind);
            end
        end
    end
    t_last_seizure(j+1:end) = (j+1:length(t_last_seizure))-seizure_inds(end);
    features(end-2,i)={t_last_seizure};
    
    % add baseline seizure frequency
    
    ind = baseline_sz_freqs.HUP_ID == ptIDs(i);
    this_sz_freq = baseline_sz_freqs.sz_per_month(ind);
    pt_baseline_sz_frequency = this_sz_freq*ones(1,nbins);
    features(end-3,i)={pt_baseline_sz_frequency};
    
    if rm_clusters
        % remove clustered seizures - more than Xhrs apart
        clustered_sz = [false; diff(seizure_times(:,1))< 3];
        cluster_inds = unique(sz_inds(clustered_sz));
        for f=1:height(features)
            features{f,i}(cluster_inds-1)=Inf; % nans are used in other features so need distinct identifier. inds -1 for bin shifting
        end
        inf_inds = cellfun(@isinf,features,'UniformOutput',false); inf_inds = inf_inds{1,i};
        
        for f = 1:height(features)
            temp_feats = features{f,i};
            temp_feats(inf_inds) =[];
            features{f,i}=temp_feats;
        end
        
    end
    
end


% get the medication counts, exclude the meds with less than 5 subjects - exclude subjects with none of the used meds
med_counts=cell(1,length(all_medications));
for n =1:length(all_medications)
    med_counts{n}=sum(contains(medications,all_medications{n}));
end
med_info = [all_medications med_counts'];
exclude_meds = [med_counts{:}]<5;
all_medications(contains(all_medications,'valproic acid'))={'valproic_acid'};

top_meds=all_medications(~exclude_meds);
med_feats=features;
med_feats(exclude_meds,:)=[];


% make table for regression
features = med_feats;
feat_names = [top_meds;{'baseline_sz_freq'};{'t_implant'}; {'ptID'}; {'seizure'}];
drug_table = table();
for i =1:length(feat_names)
    feat_name = feat_names{i};
    this_feat = features(i,:);
    this_feat = horzcat(this_feat{:})';
    drug_table.(feat_name)=this_feat;
end
drug_table.ptID = categorical(drug_table.ptID);

%% run with fake seizure times?

run_with_rand_sz = 0;
run_with_shifted_sz = 0;

if run_with_rand_sz
    for i=1:length(pts)
        pt_inds = double(drug_table.ptID) ==pts(i);
        seizures = drug_table.seizure(pt_inds);
        rand_inds = randperm(length(seizures));
        drug_table.seizure(pt_inds)= seizures(rand_inds);
    end
elseif run_with_shifted_sz
    for i=1:length(pts)
        pt_inds = double(drug_table.ptID) ==pts(i);
        seizures = drug_table.seizure(pt_inds);
        sz_inds = find(seizures==1);
        
        max_fshift = length(seizures)-sz_inds(end);
        forward_shift = round(rand(1)*max_fshift);
        
        max_bshift = sz_inds(1)-1;
        back_shift = round(rand(1)*max_bshift);
        
        shifted_inds = sz_inds + forward_shift - back_shift;
        new_sz = zeros(length(seizures),1); new_sz(shifted_inds)=1;
        drug_table.seizure(pt_inds)= new_sz;
    end
    
end

%% run bootstrapped logistic regression, get bootstrapped statistics - for null and AED model

% train model with bootstrapped samples, get bootstrapped estimates/odds ratio
iter = 1e3; % this takes awhile to run
% for ASM model
model = 1;
feat_names_new = [top_meds;{'baseline_sz_freq'};{'t_implant'}; {'ptID'}; {'seizure'}];

% save the coefficient estimates for each run
[coeff_names,coeff_stats,all_coefs]  = model_bootstrap_stats(drug_table,feat_names_new,model,iter,ptIDs); % coef stats: mean, 95-CI

disp('done with model coefficient estimates')
%for null model
null_names = [{'baseline_sz_freq'}; {'t_implant'}; {'ptID'}; {'seizure'}];
model=0;
[coeff_names_null,coeff_stats_null,~] = model_bootstrap_stats(drug_table,null_names,model,iter,ptIDs);
disp('done with null model coefficient estimates')

%% make odds ratio figure
%load('all_coeffs_and_stats_.mat');
ORs = exp(coeff_stats(:,1));
CIs = exp(coeff_stats(:,[2,3]));
inds = length(ORs):-1:1;

figure;
errorbar(ORs,inds,ORs-CIs(:,1),CIs(:,2)-ORs,'.k','horizontal','linewidth',2); hold on;
plot(ORs,inds,'.r','markersize',20);

xline(1,'--r','linewidth',1.5)

yticks(1:length(ORs))
yticklabels(coeff_names(inds))
xlabel('Odds Ratio')

sig_coeffs = coeff_stats(:,end) <0.05;
sig_coeffs_corrected = coeff_stats(:,end) <0.004;

%
plot(-1,inds(sig_coeffs),'*k');
%plot(-0.5,inds(sig_coeffs_corrected),'*b'); % bonferroni corrected


%% get AUC and operating point
% estimates

null_names = [{'baseline_sz_freq'};{'t_implant'}; {'ptID'}; {'seizure'}];

pts = 1:length(ptIDs);

pt_inds_master = cell(1,iter);
for i = 1:iter
    %random sample of 80% of patients
    pt_inds_master{i} =randperm(length(pts));
end

model = 1;
[OPs,AUCs,Xroc,Yroc,conf_matrices]  = model_bootstrap_ROC(drug_table,feat_names,model,iter,pt_inds_master);
disp('done with model AUCs')
model =0;
[OPs_null,AUCs_null,Xroc_null,Yroc_null,conf_matrices_null] =model_bootstrap_ROC(drug_table,null_names,model,iter,pt_inds_master);
disp('done with null model AUCs')

tail = 1; % one tailed p-value
AUCout = bootstrap_ci_and_p(AUCs,AUCs_null);
note = 'includes time since implant, and has all seizures (did not remove clustered seizures). null model has only baseline seizure frequency and time since implant';
% % save('all_coeffs_and_stats.mat');



