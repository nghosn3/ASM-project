%generates figures for assessing percent decrease in ASMs
clear;

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

% get overall decrease and decrease in first 24hrs
[tapered_pts,all_taper_info] = get_tapered_patients(ptIDs,all_meds);

%%

% convert doses to proportion of minimum effective dose
taper_info = cell(1,length(all_taper_info)); %with adjusted for minimum effective dose
for i=1:length(all_taper_info)
    pt_info=all_taper_info{i};
    for n=1:width(pt_info)
        if ~isempty(pt_info(n).med_name)
            med_ind = contains(aed_params.medication,pt_info(n).med_name);
            min_dose = aed_params(med_ind,:).min_dose_single_mg * aed_params(med_ind,:).min_dose_freq ;
            doses = pt_info(n).dose_schedule;
            if iscell(doses)
                doses = doses{:}; %something weird with cell array
            end
            pt_info(n).min_dose_schedule =  {doses ./ min_dose};            
        end
    end
    taper_info{i} = pt_info;
end

%use_min_dose = 1;

aed_decrease = zeros(1,length(taper_info));
aed_day1_decrease = zeros(1,length(taper_info));
for i = 1:length(taper_info)
    
    func1 = @(x) x(2);
    min_func = @(x) min(x(2:end));
    curr_pt = taper_info{i};
    
    curr_pt_doses = [curr_pt(:).min_dose_schedule];
    day1_doses = cellfun(func1,curr_pt_doses,'UniformOutput',false);
    day2_doses = cellfun(@func2,curr_pt_doses,'UniformOutput',false);
    min_doses = cellfun(min_func,curr_pt_doses,'UniformOutput',false); %first half day excluded from analysis
    
    total_dose_day1 = nansum([day1_doses{:}]);
    total_min_dose = nansum([min_doses{:}]);
    % chnage will be in units of min effective dose
    day_1_change = ([day1_doses{:}] - [day2_doses{:}]);
    aed_decrease(i) = nansum([day1_doses{:}]- [min_doses{:}]);
    aed_day1_decrease(i) = nansum(day_1_change);
end
figure;
subplot(1,2,1)
histogram(aed_day1_decrease,8);title('ASM dose decrease on day 1'); axis square; xlabel('decrease #min effective dose')
xline(mean(aed_day1_decrease)); legend([{['mean = ' num2str(mean(aed_day1_decrease))]},{['std = ' num2str(std(aed_day1_decrease))]}])

subplot(1,2,2)
histogram(aed_decrease,8);title('overall ASM dose decrease'); axis square; xlabel('decrease in #min effective dose')
xline(mean(aed_decrease)); legend([{['mean = ' num2str(mean(aed_decrease))]},{['std = ' num2str(std(aed_decrease))]}])

save_path='/Users/ninaghosn/Documents/Litt_Lab/projects/Pioneer/AED-taper-networks/results-figures/';
print([save_path 'aed_mindecrease_ofminDose.eps'],'-depsc2','-painters', '-tiff', '-r300', '-f')

aed_decrease = zeros(1,length(taper_info));
aed_day1_decrease = zeros(1,length(taper_info));
for i =  1:length(taper_info)
    
    func1 = @(x) x(2);
    min_func = @(x) min(x(2:end));
    norm_func = @(x) x./max(x);
    curr_pt = taper_info{i};
    
    curr_pt_doses = [curr_pt(:).dose_schedule];
    curr_pt_doses_norm = cellfun(norm_func,curr_pt_doses,'UniformOutput',false);
    len = min(cellfun(@length,curr_pt_doses_norm));
    curr_pt_doses_norm = cellfun(@(x)x(1:len),curr_pt_doses_norm,'UniformOutput',false);%match lengths
    [min_daily_dose,day_ind] = min_func(sum(vertcat(curr_pt_doses_norm{:}),1)); day_ind = day_ind+1; % since excluding first day in indexing    
    
    day1_doses = cellfun(func1,curr_pt_doses,'UniformOutput',false);
    day2_doses = cellfun(@func2,curr_pt_doses,'UniformOutput',false);
    
    get_min_dose = @(x) x(day_ind);
    min_doses = cellfun(get_min_dose,curr_pt_doses,'UniformOutput',false); %first half day excluded from analysis
    
    all_aed_decreases = ([day1_doses{:}]- [min_doses{:}]) ./ [day1_doses{:}];
    inf_inds = all_aed_decreases ==-Inf; all_aed_decreases(inf_inds) = 0;
    aed_decrease(i) = nanmean(all_aed_decreases);
    day1_decrease = ([day1_doses{:}] - [day2_doses{:}]) ./ [day1_doses{:}];
    inf_inds = day1_decrease ==-Inf; day1_decrease(inf_inds) = 0;
    aed_day1_decrease(i) = nanmean(day1_decrease);
    
        
end
figure;
subplot(1,2,1)
histogram(aed_day1_decrease*100,8);title('ASM dose decrease on day 1'); axis square; xlabel('% decrease on day 1')
xline(mean(aed_day1_decrease*100)); legend([{['mean = ' num2str(mean(aed_day1_decrease*100))]},{['std = ' num2str(std(aed_day1_decrease*100))]}])

subplot(1,2,2)
histogram(aed_decrease*100,8);title('overall ASM dose decrease'); axis square; xlabel('overall % decrease in daily ASM dose')
xline(mean(aed_decrease*100)); legend([{['mean = ' num2str(mean(aed_decrease*100))]},{['std = ' num2str(std(aed_decrease*100))]}])

save_path='/Users/ninaghosn/Documents/Litt_Lab/projects/Pioneer/AED-taper-networks/results-figures/';
%print([save_path 'aed_%decrease_.eps'],'-depsc2','-painters', '-tiff', '-r300', '-f')


function out = func2(x)
if length(x)>2
    out = x(3);
else
    out =NaN;
end
end