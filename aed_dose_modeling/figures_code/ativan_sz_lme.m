%% run a linear mixed model for seizures requiring ativan admin vs not - (results associated with figure 3)

cd ..
curr_path = pwd;
addpath([pwd '/DATA'])
addpath([pwd '/figures_code'])

% Get which patients have AED data, load the data
load('MAR_032122.mat')
cohort_info = readtable('HUP_implant_dates.xlsx');
ptIDs = cohort_info.ptID;
weights = cohort_info.weight_kg;

all_seizures =table(); % add ptID, seizure ID, preictal AED load, and binary yes/no for ativan admin w/in 1hr
[all_dose_curves,all_tHr,ptIDs,all_med_names,ieeg_offset,max_dur,emu_dur] = get_aed_curve_kg(ptIDs,weights);

%% find the seizures followed by ativan administration
all_times = [];
seizure_offsets = [];
start_ind =1;
for ipt=1:length(ptIDs)
    
    ptID = ['HUP' num2str(ptIDs(ipt))];
    offsets = ieeg_offset{2,ipt};
    ieeg_offset_datasets = ieeg_offset{1,ipt};
    
    %get drug curve to grab pre-ictal levels
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
    drug_sum=drug_sum./max(drug_sum); %normalize for number of drugs
    drug_sum(drug_sum==0) =NaN; %not include all zeros in average, and in histogram
    
    % get seizure times
    [seizure_times] = convert_sz_to_emu_time(offsets,ieeg_offset_datasets,ptID); % seizure times in hours of EMU stay
    
    % remove clustered seizures <2hrs apart
    cluster_diff = diff(seizure_times(:,1));
    noncluster_inds = [1; cluster_diff>1];
    seizure_times=seizure_times(logical(noncluster_inds),:);
    
    end_ind = start_ind + length(seizure_times);
    
    % convert seizure times to indices, so to minutes
    seizure_inds = round(seizure_times(:,1) *60);
    %remove seizure inds that are greater than emu recording
    %time =emu_dur(ipt)*60;%number of minutes of emu stay
    %seizure_inds(seizure_inds>time)=[];
    
    
    % med curves are sampled at one point per minute
    pt_preictal_1hr = zeros(1,length(seizure_inds));
    for i=1:length(seizure_inds)
        pt_preictal_1hr(i) = nanmean(drug_sum(seizure_inds(i)-60:seizure_inds(i)));        
    end
    
    %Get ativan times
    [~,meds,~] = parse_MAR(ptID,all_meds);
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
    all_seizures.preictal_aed_load(start_ind:end_ind-1) = pt_preictal_1hr;
    all_seizures.t_closest_ativan(start_ind:end_ind-1) =  time_to_closest_sz(:,1);
    
     start_ind = end_ind;
end

all_seizures.ativan_sz = double(all_seizures.t_closest_ativan - all_seizures.seizureEEC <= 1 &  all_seizures.t_closest_ativan - all_seizures.seizureEEC >=0);

%% run mixed effects model

% binary response of convulsive seizure
modelspec = 'ativan_sz ~ 1 + preictal_aed_load  + (1 | ptID)';
mdl = fitglme(all_seizures,modelspec,'distribution','binomial');

%%

ativan_sz_inds = all_seizures.ativan_sz == 1;
other_sz_inds = ~ativan_sz_inds;

data=[nanmean(all_seizures.preictal_aed_load(ativan_sz_inds)) nanmean(all_seizures.preictal_aed_load(other_sz_inds))];
stds = [std(all_seizures.preictal_aed_load(ativan_sz_inds)) std(all_seizures.preictal_aed_load(other_sz_inds))];

b=bar(data); hold on;
%b.FaceColor = 'flat';
% plot data points on top 
plot(.9,all_seizures.preictal_aed_load(ativan_sz_inds),'.k','markersize',10);
plot(1.9,all_seizures.preictal_aed_load(other_sz_inds),'.k','markersize',10);


er = errorbar([1 2],data,stds,'LineWidth',2);    
er.Color = [0 0 0];                            
er.LineStyle = 'none'; 

a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',14,'FontWeight','bold')
set(gca,'XTickLabelMode','auto');
xticklabels({'terminated by ativan','other seizures'});
ylabel('ASM load (normalized)')
title('ASM load','fontsize',18)

r = max([length(all_seizures.preictal_aed_load(ativan_sz_inds)) length(all_seizures.preictal_aed_load(other_sz_inds)) ]); c=2;
x = nan(r,c);
x(1:length(all_seizures.preictal_aed_load(ativan_sz_inds)),1) = all_seizures.preictal_aed_load(ativan_sz_inds);
x(1:length(all_seizures.preictal_aed_load(other_sz_inds)),2) = all_seizures.preictal_aed_load(other_sz_inds)';

[p,tbl,stats] = kruskalwallis(x,{'ativan','nonativan'},'on');
%multcompare(stats)





