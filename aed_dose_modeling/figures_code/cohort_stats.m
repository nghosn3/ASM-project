%% get cohort information - for figure 2 and supplement

close all;clear;

cd ..
curr_path = pwd;
addpath([pwd '/DATA'])
addpath([pwd '/figures_code'])
% loadthe aed metadata
aed_params = readtable('AED_metadata.xlsx');
aed_params.medication = arrayfun(@lower,aed_params.medication);

% Get which patients have AED data, load the data
cohort_info = readtable('HUP_implant_dates.xlsx');
ptIDs = cohort_info.ptID;
weights = cohort_info.weight_kg;

%  ptIDs = ptIDs(ptIDs==146);
%  weights = weights(ptIDs==146);

% get their medication data 
[all_dose_curves,all_Hr,ptIDs,all_med_names,ieeg_offset] = get_aed_curve_kg(ptIDs,weights); %EDIT!! TESTING SCRIPT

% get the seizure and SOZ localization information 
soz_info = readtable('Erin_szs_times.xlsx','VariableNamingRule','preserve','Sheet','SOZ');

% load home medications
load('MAR_032122.mat')

home_meds =  readtable('all_procedure_times_discharge_meds.xlsx','Sheet','Order Meds','VariableNamingRule','preserve');
emu_med_names = cellfun(@lower,home_meds.MED_NAME,'UniformOutput',false);
emu_med_names = cellfun(@strsplit,emu_med_names,'UniformOutput',false);
emu_med_names = cellfun(@(x)x{1},emu_med_names,'UniformOutput',false);

[gen_names] = get_AEDs_from_list(emu_med_names);
del_inds = cell2mat(cellfun(@isempty,gen_names,'UniformOutput',false));
home_meds(del_inds,:)=[];
home_meds.medication = gen_names(~del_inds);

%% medications in emu 
num_meds = zeros(1,length(ptIDs));
for i=1:length(ptIDs)
    med_names = all_med_names{i};
    med_names(contains(med_names,'lorazepam'))=[];
    num_meds(i)=length(med_names);
end

medications = vertcat(all_med_names{:});
all_medications = unique(medications);
% get the medication counts, exclude the meds with less than 5 subjects - exclude subjects with none of the used meds
med_counts=cell(1,length(all_medications));
for n =1:length(all_medications)
    med_counts{n}=sum(contains(medications,all_medications{n}));
end
med_info = [all_medications med_counts'];

%% seizures in the EMU 
num_sz = zeros(1,length(ptIDs));
for i=1:length(ptIDs)
    ptID = ['HUP' num2str(ptIDs(i))];
    [seizure_times,seizure_dataset] = get_seizure_times_from_sheet(ptID);
    num_sz(i)=height(seizure_times);
end

%% localization 
mtl = 0;
neocortical = 0;
other = 0;
for i=1:length(ptIDs)
    ptID = ['HUP' num2str(ptIDs(i))];
    pt_soz_info = soz_info(strcmp(soz_info.name,ptID),:);
    if contains(pt_soz_info.region,'mesial temporal')
        mtl=mtl+1;
    elseif contains(pt_soz_info.region,'neocort')
        neocortical=neocortical+1;
    else 
        other = other+1;
    end 
end 

%% home medications as proportion of minimum effective dose (proxy for overall load)
home_load = zeros(1,length(ptIDs));
for i = 1:length(ptIDs)
    ptID = ['HUP' num2str(ptIDs(i))];
    [med_names,meds,explant_date] = parse_MAR(ptID,all_meds);
    
    % grab home meds for pt
    pt_home_meds = home_meds(home_meds.HUP_ID==str2double(ptID(4:end)),:);
    % only for meds within 3days of explant date
    pt_home_meds= pt_home_meds(pt_home_meds.EVENT_TIME>= explant_date & pt_home_meds.EVENT_TIME<=explant_date+days(4),:);
    [~,ia] = unique(pt_home_meds.medication);
    pt_home_meds = pt_home_meds(ia,:);
    
    % find the dose of the home meds from the EMU meds if NaN
    %add meds given in first day ofstay
    missed_meds=cell(1,length(med_names));
    ind=1;
    if ~isempty(pt_home_meds)
        for j =1:length(med_names)
            if ~contains(pt_home_meds.medication,med_names(j))
                missed_meds{ind} =med_names(j);
            end
        end
    else
        pt_home_meds(1:length(med_names),'medication') = med_names;
    end
    
    % if dose for med not included in data
    home_med_names = unique(pt_home_meds.medication);
    for m=1:length(home_med_names)
        if isnan(pt_home_meds.DOSE(m)) || pt_home_meds.DOSE(m)==0
            temp_med_doses = [meds.dose(contains(meds.medication,home_med_names(m))); 0];
            pt_home_meds.dose(m)= temp_med_doses(1);
        else 
            pt_home_meds.dose(contains(pt_home_meds.medication,home_med_names(m)))= pt_home_meds.DOSE(m);
        end 
        % get medication frequency
        if ~isempty(pt_home_meds.SIG{m})
            [spd] = parse_med_frequencies(pt_home_meds.SIG(m));
            pt_home_meds.frequency(m) = spd;
        else 
            pt_home_meds.frequency(m)=NaN;
        end
    end 
    
    % calculate the load
    aed_load=0;
    for n=1:height(pt_home_meds)
        med_ind = contains(aed_params.medication,pt_home_meds.medication{n});
        if sum(med_ind)>0
            min_dose = aed_params.min_dose_single_mg(med_ind)*aed_params.min_dose_freq(med_ind);
            if ~isnan(pt_home_meds.frequency(n))
                min_freq=pt_home_meds.frequency(n);
            else 
                min_freq = aed_params.min_dose_freq(med_ind);
            end 
            dose = pt_home_meds.dose(n)*min_freq;
            if ~isnan(dose)
                aed_load=[aed_load  (dose./min_dose)];
            else
                disp([newline 'HUP' num2str(ptIDs(i)) ' ' pt_home_meds.medication{n} ' dose is nan: ' 'dose: ' num2str(pt_home_meds.dose(n)) ', frequency: ' num2str(pt_home_meds.frequency(n))])
            end
        end
        home_load(i)=nansum(aed_load);
    end
end 

figure;
histogram(home_load); hold on;
xline(mean(home_load)); legend([{['mean = ' num2str(mean(home_load))]},{['std = ' num2str(std(home_load))]}])
title('AED load at admission','fontsize',16)
xlabel('initial dose/min effective dose','fontsize',14);
ylabel('# patients','fontsize',14);
save_path='/Users/ninaghosn/Documents/Litt_Lab/projects/Pioneer/AED-taper-networks/results-figures/';
print([save_path 'home_aed_load_hist.eps'],'-depsc2','-painters', '-tiff', '-r300', '-f')

%% average number of drugs tapered 
num_tapered = zeros(1,length(ptIDs));
[tapered_pts,all_taper_info] = get_tapered_patients(ptIDs,all_meds);
tapered_ptIDs = ptIDs(logical(tapered_pts));


for i=1:length(ptIDs)
    taper_info=all_taper_info{i};
    if isfield(taper_info,'tapered')
        num_tapered(i)=sum([taper_info.tapered]);
    end
end

avg_num_tapered= mean(num_tapered)
std_num_tapered = std(num_tapered)

% for patients on less than or equal to 2 meds
aed_2 = num_meds <= 2;
aed_3_4 = num_meds > 2 & num_meds <=4;
aeds_5 = num_meds>4;

[mean(num_tapered(aed_2)) std(num_tapered(aed_2))]
[mean(num_tapered(aed_3_4)) std(num_tapered(aed_3_4))]
[mean(num_tapered(aeds_5)) std(num_tapered(aeds_5))]

% number of seizures for patients that were tapered vs not
taper_sz=num_sz(logical(tapered_pts))
nontaper_sz = num_sz(~logical(tapered_pts))

[p,h,stats] = ranksum(taper_sz,nontaper_sz)
figure;
subplot(1,2,1); boxplot(taper_sz); title('tapered seizures');ylim([0 100])
subplot(1,2,2);boxplot(nontaper_sz);title('nontapered seizures'); ylim([0 100])

%% baseline seizure frequency
baseline_sz_freqs = readtable('no_phi_baseline_sz_freq.xlsx');

tapered_ptIDs = ptIDs(logical(tapered_pts));

% get frequencies for ptIDs
baseline_sz_freqs.tapered = nan(height(baseline_sz_freqs),1);
for i=1:length(ptIDs)
    idx = baseline_sz_freqs.HUP_ID == ptIDs(i);
    baseline_sz_freqs.num_sz(idx) = num_sz(i);
    baseline_sz_freqs.tapered(idx) = tapered_pts(i);
    
end

tapered_freqs = baseline_sz_freqs.sz_per_month(baseline_sz_freqs.tapered==1);
nontapered_freqs = baseline_sz_freqs.sz_per_month(baseline_sz_freqs.tapered==0);


[p_baseline,h,stats] = ranksum(taper_sz,nontaper_sz);

figure;
x = log10(baseline_sz_freqs.sz_per_month);
y = log10(baseline_sz_freqs.num_sz);
taper_inds = (baseline_sz_freqs.tapered);

inf_inds =y==-Inf;
y(inf_inds)=[]; x(inf_inds)=[]; taper_inds(inf_inds)=[];
zero_inds = x==0;
x(zero_inds)=[]; y(zero_inds)=[]; taper_inds(zero_inds)=[];

taper_inds=logical(taper_inds);


%plot tapered patients in red
plot(x(taper_inds),y(taper_inds),'.r','markersize',25); hold on;
plot(x(~taper_inds),y(~taper_inds),'.b','markersize',25);
ylim([0 2.5]);xlim([0 4])
xlabel('log(baseline seizure frequency)');ylabel('log(#seizures in EMU)')
axis square;
coefs = polyfit(x',y',1);
pred = (coefs(1)*x) +coefs(2); hold on;
plot(x,pred,'k','Linewidth',2);
[r,p] = corr([x y]);
r = r(2,1);
p = p(2,1);
title(['p: ' num2str(p) ', R: ' num2str(r)]);
legend([{'taper'},{'no taper'}])



