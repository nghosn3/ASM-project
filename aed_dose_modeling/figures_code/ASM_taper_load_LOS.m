close all;clear;

% load the aed metadata
aed_params = readtable('AED_metadata.xlsx');
aed_params.medication = arrayfun(@lower,aed_params.medication);

% Get which patients have AED data, load the data
cohort_info = readtable('HUP_implant_dates.xlsx');
ptIDs = cohort_info.ptID;
weights = cohort_info.weight_kg;

% load medications
load('MAR_032122.mat')

[all_dose_curves,all_tHr,ptIDs,all_med_names,ieeg_offset,max_dur,emu_dur] = get_aed_curve_kg(ptIDs,weights,0,0); 

%% get which seizures are followed by ativan and the drug levels
all_seizures =table(); % add ptID, seizure ID, preictal AED load, and binary yes/no for ativan admin w/in 1hr

%% find the seizures followed by ativan administration
aed_decrease = zeros(3,length(ptIDs));

all_times = [];
seizure_offsets = [];
start_ind =1;
for ipt=1:length(ptIDs)
    
    ptID = ['HUP' num2str(ptIDs(ipt))];
    offsets = ieeg_offset{2,ipt};
    ieeg_offset_datasets = ieeg_offset{1,ipt};
    
    %get drug curve to grab daily levels
    [med_names,~,~] = parse_MAR(ptID,all_meds);
    
    % get the total AED dose over time
    drugs =zeros(length(med_names),round(emu_dur(ipt)*60)); %450 hours of EMU stay in minutes
    all_dstarts = [];
    for i =1:length(med_names)
        drug=all_dose_curves{ipt}{i};
        drug=drug./nanmax(drug); %normalize each drug curve
        if ~isempty(drug)
            dStart = round(all_tHr{ipt}{i}(1)*60)-1;
            all_dstarts = [all_dstarts dStart];
            drugs(i,dStart+1:dStart+length(drug))=drug;
        end
    end
    drug_sum=nansum(drugs,1);
    drug_sum=drug_sum./(length(med_names)); %normalize for number of drugs
    drug_sum(drug_sum==0) =NaN; %not include all zeros in average, and in histogram
    
    % med curves are sampled at one point per minute
    t0 = min(all_dstarts); %earliest time in min a drug was administered
    t1 = t0 + (24*60);
    t2 = t1 + (24*60);
    t3 = t2 + (24*60);
    t4 = t3 + (24*60);
    
    day1_decrease = (nanmean(drug_sum(t0:t1))-nanmean(drug_sum(t1:t2)))./ (nanmean(drug_sum(t0:t1)));
    
    % some patients have short EMU stays
    try
    day2_decrease = (nanmean(drug_sum(t0:t1))-nanmean(drug_sum(t2:t3)))./ (nanmean(drug_sum(t0:t1)));
    catch
        day2_decrease =NaN;
    end 
    
    try
    day3_decrease = (nanmean(drug_sum(t0:t1))-nanmean(drug_sum(t3:t4)))./ (nanmean(drug_sum(t0:t1)));
    catch
        day3_decrease =NaN;
    end 
    
    aed_decrease(1,ipt) = day1_decrease;
    
    % get seizure times
    [seizure_times] = convert_sz_to_emu_time(offsets,ieeg_offset_datasets,ptID); % seizure times in hours of EMU stay
    
    end_ind = start_ind + length(seizure_times);
    
    % convert seizure times to indices, so to minutes
    seizure_inds = round(seizure_times(:,1) *60);
    
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


%% run model for length of stay
tbl = table();
tbl.asm_decrease = aed_decrease(1,:)'; %decrease from day 1 to 2
tbl.baseline_sz_freq = (sz_freqs);%./max(sz_freqs); % normalize like other feature


tbl2 = tbl;
tbl2.length_stay = hours(cohort_info.Explant_Date - cohort_info.Implant_Date);
zero_inds = tbl2.asm_decrease == -Inf ; % | tbl2.asm_decrease < 0
tbl2(zero_inds,:)=[];
%outliers = find(isoutlier(tbl2.baseline_sz_freq,'quartiles'));
%tbl2(39,:)=[];
mdl_LOS = fitlm(tbl2)

%% plot stuff 
figure;
subplot(1,2,1);
x1=tbl2.asm_decrease;
x2 = tbl2.baseline_sz_freq;
y= tbl2.length_stay;
plot3(x1,x2,y,'.k','markersize',15); axis square; hold on;
 
ylim([0 300])
title('linear model for length of stay')
xlabel('ASM decrease');ylabel('baseline seizure frequency');zlabel('length of stay (hrs)');

% plot mesh on top
b = [mdl_LOS.Coefficients.Estimate];
x1fit = linspace(min(x1),max(x1),50);
x2fit = linspace(min(x2),max(x2),50);
[X1FIT,X2FIT] = meshgrid(x1fit,x2fit);
YFIT = b(1) + b(2)*X1FIT + b(3)*X2FIT;
mesh(X1FIT,X2FIT,YFIT)


subplot(1,2,2)
decrease_had_conv =x1(has_conv);
decrease_no_conv =x1(~has_conv);
len = max([length(decrease_had_conv) length(decrease_no_conv)]);
data = nan(2,len);
data(1,1:length( decrease_had_conv))= decrease_had_conv;
data(2,1:length( decrease_no_conv))= decrease_no_conv; 

boxplot(data',[{'had ativan sz'},{'no ativan sz'}]); axis square;
[p,h,stats]=ranksum(decrease_had_conv,decrease_no_conv)
save_path='/Users/ninaghosn/Documents/Litt_Lab/projects/Pioneer/AED-taper-networks/results-figures/';
%print([save_path 'fig04_los_conv.eps'],'-depsc2','-painters', '-tiff', '-r300', '-f')



function out = func2(x)
if length(x)>2
    out = x(3);
else
    out =NaN;
end
end