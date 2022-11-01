clear;
cd ..
curr_path = pwd;
addpath([pwd '/DATA'])
addpath([pwd '/figures_code'])

% load the aed metadata
aed_params = readtable('AED_metadata.xlsx');
aed_params.medication = arrayfun(@lower,aed_params.medication);

% Get which patients have AED data, load the data
cohort_info = readtable('HUP_implant_dates.xlsx');
ptIDs = cohort_info.ptID;
weights = cohort_info.weight_kg;

% load home medications
home_meds =  readtable('all_procedure_times_discharge_meds.xlsx','Sheet','Order Meds','VariableNamingRule','preserve');
load('MAR_032122.mat')

[all_dose_curves,all_tHr,ptIDs,all_med_names,ieeg_offset,max_dur,emu_dur] = get_aed_curve_kg(ptIDs,weights,0,0);

%% get which seizures are followed by ativan and the drug levels
all_seizures =table(); % add ptID, seizure ID, preictal AED load, and binary yes/no for ativan admin w/in 1hr

%% find the seizures followed by ativan administration
medications = vertcat(all_med_names{:});
all_medications = unique(medications);

aed_decrease = zeros(length(ptIDs),length(all_medications));

all_times = [];
seizure_offsets = [];
start_ind =1;

for ipt=1:length(ptIDs)
    
    ptID = ['HUP' num2str(ptIDs(ipt))];
    offsets = ieeg_offset{2,ipt};
    ieeg_offset_datasets = ieeg_offset{1,ipt};
    
    %get drug curves to grab daily levels
    [med_names,meds,~] = parse_MAR(ptID,all_meds);
    
    % find the earliest time a drug was administered, and align all drugs
    drugs =zeros(length(med_names),max_dur*60); %450 hours of EMU stay in minutes
    all_dstarts=[];
    for i =1:length(med_names)
        drug=all_dose_curves{ipt}{i};
        drug=drug./nanmax(drug); %normalize each drug curve
        if ~isempty(drug)
            dStart = round(all_tHr{ipt}{i}(1)*60)-1;
            all_dstarts=[all_dstarts dStart];
            drugs(i,dStart+1:dStart+length(drug))=drug;
        end
    end
    
    for n = 1:length(med_names)
        med_ind = find(strcmp(med_names(n),all_medications)); % idx of medication in table
        drug = drugs(n,:);
        drug=drug./nanmax(drug); %normalize each drug curve
        
        % med curves are sampled at one point per minute
        t0 = min(all_dstarts); %earliest time in min a drug was administered
        t1 = t0 + (24*60);
        t2 = t1 + (24*60);
        
        day1_decrease = (nanmean(drug(t0:t1))-nanmean(drug(t1:t2)))./ (nanmean(drug(t0:t1)));
        aed_decrease(ipt,med_ind) = day1_decrease;
    end
    
    % get seizure times
    [seizure_times] = convert_sz_to_emu_time(offsets,ieeg_offset_datasets,ptID); % seizure times in hours of EMU stay
    
    end_ind = start_ind + length(seizure_times);
    
    % convert seizure times to indices, so to minutes
    seizure_inds = round(seizure_times(:,1) *60);
    
    %Get ativan times
    ativan_inds = strcmp(meds.medication,'lorazepam');
    times = meds.admin_time(ativan_inds);
    
    time_to_closest_sz = nan(length(seizure_times(:,1)),2);
    for n=1:length(seizure_times(:,1))
        sz_diffs = seizure_times(n,1)-times;
        before_ativan = sz_diffs< 0;
        if ~isempty(times)
            if ~isempty(before_ativan) && ~(sum(before_ativan)==0)
                [mval,~]=min(abs(sz_diffs(before_ativan)));
                ind = find(abs(sz_diffs)==mval);
                time_to_closest_sz(n,:)=[times(ind(1)) seizure_times(n,1)];
            end
        end
    end
    all_seizures.ptID(start_ind:end_ind-1)=ptIDs(ipt)*ones(1,end_ind-start_ind);
    all_seizures.seizureEEC(start_ind:end_ind-1) = seizure_times(:,1);
    all_seizures.t_closest_ativan(start_ind:end_ind-1) =  time_to_closest_sz(:,1);
    
    start_ind = end_ind;
end


all_seizures.ativan_sz = double(all_seizures.t_closest_ativan - all_seizures.seizureEEC <= 1 &  all_seizures.t_closest_ativan - all_seizures.seizureEEC >=0);

% get binary of that patient had a convulsion or not
has_conv = false(length(ptIDs),1);
for i=1:length(ptIDs)
    pt_inds = all_seizures.ptID == ptIDs(i);
    has_conv(i) = any(all_seizures.ativan_sz(pt_inds));
end

% get baseline seizure frequencies
baseline_sz_freqs = readtable('no_phi_baseline_sz_freq.xlsx');
sz_freqs = zeros(length(ptIDs),1);
for i = 1:length(ptIDs)
    ind = baseline_sz_freqs.HUP_ID == ptIDs(i);
    sz_freqs(i) = baseline_sz_freqs.sz_per_month(ind);
end


%% run model for severity
tbl = table();
for n=1:length(all_medications)
    name = all_medications(n);
    tbl.(name{:})=aed_decrease(:,n);
end

% remove lorazepam 
ativan_ind = (strcmp(all_medications,'lorazepam'));
tbl(:,ativan_ind)=[];
all_medications(ativan_ind)=[];

med_counts=cell(1,length(all_medications));
for n =1:length(all_medications)
    med_counts{n}=sum(contains(medications,all_medications{n}));
end

med_info = [all_medications med_counts'];
exclude_meds = [med_counts{:}]<20;
all_medications(contains(all_medications,'valproic acid'))={'valproic_acid'};

top_meds=all_medications(~exclude_meds);
tbl(:,exclude_meds)=[];

% add baseline seizure frequency
tbl.baseline_sz_freq = sz_freqs./max(sz_freqs); % normalize like other feature
tbl1 = tbl;
tbl1.has_conv = double(has_conv);
%tbl.time_to_first_seizure = time_to_first_seizure;

% find -inf and negative values and make them zero for no decrease
% zero_inds = tbl1.asm_decrease == -inf | tbl1.asm_decrease <0;
% tbl1.asm_decrease(zero_inds)=0;


%mdl_severity = fitglm(tbl1,'distribution','binomial')

%% run model for length of stay
tbl2 = tbl;
tbl2.length_stay = hours(cohort_info.Explant_Date - cohort_info.Implant_Date);
[r,c] = find(table2array(tbl2(:,:) )== -Inf) ; % | tbl2.asm_decrease < 0
tbl2(r,c)={0};
%outliers = find(isoutlier(tbl2.baseline_sz_freq,'quartiles'));
%tbl2(39,:)=[];
mdl_LOS = fitlm(tbl2)

%% plot stuff
figure;

coefs = mdl_LOS.Coefficients.Estimate(2:end);
SEs = mdl_LOS.Coefficients.SE(2:end);
ORs = (coefs);
CIs = ([coefs+(1.96*SEs) coefs-(1.96*SEs)]);
coef_names=mdl_LOS.CoefficientNames(2:end);
inds = length(ORs):-1:1;

neg = ORs-CIs(:,2);
pos = CIs(:,1)-ORs;
errorbar(ORs,inds,neg,pos,'.k','horizontal','linewidth',2); hold on;
plot((ORs),inds,'.r','markersize',20); axis square;
xline(1,'--r','linewidth',1.5)
yticks(1:length(ORs))
yticklabels(coef_names(inds))
xlabel('model coefficient');
ylim([-0.5 6])
title('variable coefficients for predicting length of stay ')

model_data = [ORs CIs]