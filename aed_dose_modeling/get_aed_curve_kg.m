function [all_dose_curves,all_tHr,ptIDs,all_med_names,all_ieeg_offset,max_dur,emu_dur] = get_aed_curve_kg(ptIDs,weights,use_min_dose,for_borel)

curr_path = pwd;
addpath([pwd '/DATA'])

if ~exist('use_min_dose','var')
    % third parameter does not exist, so default it to something
    use_min_dose = 0;
end

if ~exist('for_borel','var')
    % third parameter does not exist, so default it to something
    for_borel = 0;
end

%turn warning off
id = 'MATLAB:table:RowsAddedExistingVars';
warning('off',id)

% load all med data
load('MAR_032122.mat');

% load AED parameters and
aed_params = readtable('AED_metadata.xlsx');
aed_params.medication = arrayfun(@lower,aed_params.medication);

% load lab and home med data- data with med names to convert - update spreadsheets when full patient cohort arrives
home_meds =  readtable('all_procedure_times_discharge_meds.xlsx','Sheet','Order Meds','VariableNamingRule','preserve');
implant_info = readtable('all_procedure_times_discharge_meds.xlsx','VariableNamingRule','preserve');
implant_info.implant_time=hours(timeofday(implant_info.date_time));

% clean up medication name field in each data sheet:
temp_list = cellfun(@strsplit,home_meds.MED_NAME,'UniformOutput',false);
for i =1:length(temp_list)
    home_meds.MED_NAME(i) = temp_list{i}(1);
end

%get AEDs only
emu_med_names = cellfun(@lower,home_meds.MED_NAME,'UniformOutput',false);
emu_med_names = cellfun(@strsplit,emu_med_names,'UniformOutput',false);
emu_med_names = cellfun(@(x)x{1},emu_med_names,'UniformOutput',false);

[gen_names] = get_AEDs_from_list(emu_med_names);
del_inds = cell2mat(cellfun(@isempty,gen_names,'UniformOutput',false));
home_meds(del_inds,:)=[];
home_meds.medication = gen_names(~del_inds);

%% need to filter home_meds for only meds given at discharge, then calculate their frequency - use implant and explant dates 

%initialize dose curves for EMU stay
all_dose_curves =cell(1,length(ptIDs));
all_tHr =cell(1,length(ptIDs));
all_med_names=cell(1,length(ptIDs));

% get the ieeg_offset for plotting and the dataset the offset is for
all_ieeg_offset=cell(3,length(ptIDs));

emu_dur =zeros(1,length(ptIDs));
for ipt = 1:length(ptIDs)
    % patient medication administration
    ptID = ['HUP' num2str(ptIDs(ipt))];
    pt_weight = weights(ipt);
    disp(ptID);
    [med_names,meds,explant_date] = parse_MAR(ptID,all_meds);
    all_med_names(ipt) = {med_names};
    emu_dur(ipt) = hours(meds.date(end)-meds.date(1)) +24; %approx
    
    implant_date=implant_info(implant_info.HUP==str2double(ptID(4:end)),:);
    implant_date=implant_date(1,'date_time');
    
    %get offset in ieeg
    eeg_diff = (meds.admin_time*3600 - meds.OffsetSecondsInIeeg);
    eeg_round = unique(round(eeg_diff*10^5)/10^5);
    eeg_offsets = eeg_round(~isnan(eeg_round));
    
    eeg_datasets = unique(meds.dataset);
    eeg_datasets = eeg_datasets(contains(eeg_datasets,'D'));
    all_ieeg_offset(1,ipt) = {eeg_datasets};
    all_ieeg_offset(2,ipt) = {eeg_offsets};
    
    datasetStarts = unique(meds.DatasetStartIeeg);
    datasetStarts = datasetStarts(~isnat(unique(meds.DatasetStartIeeg)));
    all_ieeg_offset(3,ipt) = {datasetStarts};
    
    % initialize BPL curves
    c_all=cell(1,length(med_names));
    t_all=cell(1,length(med_names));
    
    % grab home meds for pt
    pt_home_meds = home_meds(home_meds.HUP_ID==str2double(ptID(4:end)),:);
    % only for meds within 3days of explant date
    pt_home_meds= pt_home_meds(pt_home_meds.EVENT_TIME>= explant_date & pt_home_meds.EVENT_TIME<=explant_date+days(4),:);
    [~,ia] = unique(pt_home_meds.medication);
    pt_home_meds = pt_home_meds(ia,:);
    
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
    
    %% if dose for med not included in data
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
        
    for n = 1:length(med_names)
        %disp(med_names{n});
        med_name_ind = contains(med_names,med_names{n}); % which med
        this_med_info = meds(contains(meds.medication,med_names(med_name_ind)),:); % get info of just that med
        
        % delete NaN or empty entries?
        if sum(isnan(this_med_info.dose)) >0
            if contains(this_med_info(isnan(this_med_info.dose),:).medication,'levetiracetam')
                this_med_info(isnan(this_med_info.dose),'dose')={1500};
            end 
            nan_med = this_med_info(isnan(this_med_info.dose),'medication');
            disp([ptID ' has medications with NaN dose: ' nan_med.medication{:}] )
            this_med_info(isnan(this_med_info.dose),:)=[];
        end
        
        %placeholder for adding the last timepoint in the emu stay, to get
        %the concentration to the end....using 450 can change later
        max_dur =450;
        dur =emu_dur(ipt);
        this_med_info.admin_time(end+1) = dur; %hrs
        this_med_schedule = diff(this_med_info.admin_time); %in hours
        
        % grab drug parameters for model:
        med_ind =contains(aed_params.medication,med_names(n));
        
        % if med_ind is zero, then drug not exist or its wine
        if sum(med_ind)==0
            disp(['warning: drug ' med_names{n} ' not found...']);
        else
            F = aed_params(med_ind,:).F;
            vd=aed_params(med_ind,:).vd * pt_weight; % EDIT 4/28 - include patient weight for vd
            ka=aed_params(med_ind,:).ka;
            tmax = aed_params(med_ind,:).t_max;
            % take mean of range for current model:
            tHalf = [aed_params(med_ind,:).t_half_e]; tHalf=strsplit(tHalf{1},',');tHalf = cellfun(@str2double,tHalf);
            tHalf = mean(tHalf);
            min_dose = aed_params(med_ind,:).min_dose_single_mg;
            
            %% get the initial concentration at EMU admission following surgery
            ind=contains(pt_home_meds.medication,med_names(n));
            if contains(pt_home_meds.medication(ind),'lorazepam')
                c0=0;
            elseif sum(ind)~=0 && ~isnan(pt_home_meds.dose(ind))
                dose = pt_home_meds.dose(ind); % later take into account nonuniform dose schedule
                frequency = pt_home_meds.frequency(ind);
                if isnan(pt_home_meds.frequency(ind))
                    frequency = 2; % ifthere is a dose but not a frequency, then say frequency is 2
                end
                % need to use implant date and time and date/time of first implant to pick c0 from c_ss
                [~,css]=get_c0(dose,frequency,tHalf,F,vd,tmax,ka); %this is curve of one day
                if this_med_schedule(1)==0 % if meds admisitered together
                    surg_tInt =this_med_schedule(2);
                else 
                    surg_tInt = this_med_schedule(1);
                end
                c0=css(round((mod(surg_tInt,24./frequency))*60)); %this is the c0 at the first admin,indexed from t=0 (midnight) so need to adjust based on sleep/wake estimate
            else
                c0=0;
            end
            
            %%
            c_drug=[];
            dose = this_med_info.dose(1);
            i=2;
            while i <= height(this_med_info)
                
                tInt = this_med_schedule(i-1);
                if tInt < 0.08 % 5 minutes apart
                    dose = dose + this_med_info.dose(i+1);
                    tInt = this_med_schedule(i);
                    i=i+1;
                end
                if use_min_dose
                    [c_t,~]= get_single_dose_curve(c0,dose./min_dose,tHalf,tInt,F,vd,tmax,ka);
                else
                    [c_t,~]= get_single_dose_curve(c0,dose,tHalf,tInt,F,vd,tmax,ka);
                end
                c_drug=[c_drug c_t];
                c0 = c_drug(end);
                dose = this_med_info.dose(i);
                i=i+1;
            end
            
            c_all(n) = {c_drug};
            all_dose_curves(ipt) ={c_all};
            
            tHr = linspace(this_med_info.admin_time(1),this_med_info.admin_time(end),length(c_drug));
            t_all(n) = {tHr};
            all_tHr(ipt) = {t_all};
            
            
        end
    end
end


end