%% make figures!!

clear;

load('MAR_032122.mat')

% Get which patients
cohort_info = readtable('HUP_implant_dates.xlsx');
ptIDs = cohort_info.ptID;
weights = cohort_info.weight_kg;

% using the AED BPL model
[all_dose_curves,all_tHr,ptIDs,all_med_names,ieeg_offset,max_dur,emu_dur] = get_aed_curve_kg(ptIDs,weights);

% using only interpolation dosage
%[all_dose_curves,all_tHr,all_med_names,ieeg_offset,emu_dur] = get_dose_schedules(ptIDs);


% load AED parameters and
aed_params = readtable('AED_metadata.xlsx');
aed_params.medication = arrayfun(@lower,aed_params.medication);

%% 1A: single admin of different drugs: ativan, Keppra, Onfi or zonisamide, and CBZ
med_names = [{'lorazepam'},{'levetiracetam'},{'clobazam'},{'carbamazepine'}];
doses = [2 1500 10 400];
curves=cell(1,length(med_names));
for n=1:length(med_names)
    med_ind =(contains(aed_params.medication,med_names(n)));
    %get AED parameters
    F = aed_params(med_ind,:).F;
    vd=aed_params(med_ind,:).vd;
    ka=aed_params(med_ind,:).ka;
    tmax = aed_params(med_ind,:).t_max;
    % take mean of range for current model:
    tHalf = [aed_params(med_ind,:).t_half_e]; tHalf=strsplit(tHalf{1},',');tHalf = cellfun(@str2double,tHalf);
    tHalf = mean(tHalf);
    tInt = 24;
    dose = doses(n);
    c0=0;
    min_dose = aed_params(med_ind,:).min_dose_single_mg;
    %getAED curve
    [c_t,t] = get_single_dose_curve(c0,dose,tHalf,tInt,F,vd,tmax,ka);
    curves{n}=c_t;
end

axis square;
figure;
colors = lines(length(med_names));
for n=1:length(curves)
    plot(t,curves{n}./nanmax(curves{n}),'LineWidth',2,'Color',[colors(n,:) 0.6]); hold on;
end

legend_names= [{'lorazepam: 2mg'},{'levetiracetam: 1500mg'},{'clobazam: 10mg'},{'carbamazepine: 400mg'}];
legend(legend_names);
xlabel('time (hrs)','FontSize',14)
ylabel('blood plasma level relative to minimum dose','FontSize',14)
title('Single drug administration curve','FontSize',16);

%% for results associated with figure 3: rank of seizure preictal aed load relative to rest, binomial test
% Get the pre-ictal AED load times

[preictal_aed_load,null_aed_loads] = get_avg_preictal_levels(ptIDs, all_meds, all_dose_curves, all_tHr, ieeg_offset,max_dur,emu_dur);

% get the aed load for all time bins:
aed_loads=cell(1,length(ptIDs));
use_taper_patients = 0;
if use_taper_patients
    tapered_pts = get_tapered_patients(ptIDs,all_meds);
end

for i=1:length(ptIDs)
    pt_curves = all_dose_curves{i};
    med_names=all_med_names{i};
    tHr = all_tHr{i};
    
    drugs =zeros(length(med_names),ceil(emu_dur(i)*60)); %450 hours of EMU stay in minutes
    for n =1:length(med_names)
        drug=pt_curves{n};
        drug=drug./nanmax(drug); %normalize each drug curve
        if ~isempty(drug)
            dStart = round(tHr{n}(1)*60)-1;
            drugs(n,dStart+1:dStart+length(drug))=drug;
        end
    end
    drug_sum=nansum(drugs,1);
    drug_sum = drug_sum./length(med_names); %   RERUN FOR NORMALIZING 9/07/22
    % cut off drug curve to only be length of emu stay
    time =emu_dur(i)*60;%number of minutes of emu stay
    if time<length(drug_sum)
        drug_sum = drug_sum(1:time);
    end
    drug_sum(drug_sum==0)=[]; %delete rather than NaN it 
    
    % average in one hour bins
    nbins=ceil(length(drug_sum)./60);
    ind =1;
    drug_sum_wins = zeros(1,nbins);
    for j = 1:60:length(drug_sum)-60
        drug_sum_wins(ind)=mean(drug_sum(j:j+60));
        ind=ind+1;
    end
    aed_loads{i}=drug_sum_wins;
    
end
figure()
%[out,pval_binom_alt,successes] = plot_orders(aed_loads(logical(tapered)),preictal_aed_load(logical(tapered)));
subplot(1,2,1)
inds= ~cellfun(@isempty,preictal_aed_load); 
[out,pval_binom_alt,successes,successes_alt] = plot_orders(aed_loads(inds),preictal_aed_load(inds));
axis square;

subplot(1,2,2)
median_loads = cellfun(@median,aed_loads(inds));
sz_median_loads = cellfun(@median,preictal_aed_load(inds));
boxplot([median_loads' sz_median_loads'],[1 2]);hold on;
plot(0,0);
plot(1,median(median_loads),'.');
plot(2,mean(sz_median_loads),'.');
axis square;
xticklabels([{'median ASM load'},{'median pre-ictal ASM load'}]);
legend(['binomial test pval = ' num2str(pval_binom_alt)],['median ASM load = ' num2str(median(median_loads))],['median sz load = ' num2str(median(sz_median_loads))])


figure;axis square;
sz_means = horzcat(preictal_aed_load{:});
null_means = vertcat(null_aed_loads{:});
histogram(sz_means,'normalization','pdf');hold on;
histogram(null_means,'normalization','pdf')
legend({'pre-ictal AED load','null distribution'})
