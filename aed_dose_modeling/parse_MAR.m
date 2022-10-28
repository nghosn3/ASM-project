%% parse medication info from MAR (2/24/22) sheet 
function [med_names,meds,explant_date] = parse_MAR(ptID,all_meds)

pt_inds = all_meds.HUP_ID== str2double(ptID(4:end));
meds=all_meds(pt_inds,:);

% load implant dates for all patients and find info for ptID
implant_info=readtable('HUP_implant_dates.xlsx','VariableNamingRule','preserve');
implant_info.implant_time=implant_info.implant_time*24;

pt_implant = implant_info(implant_info.ptID == str2double(ptID(4:end)),:);
%if explant date is missing, assume they stayed for two weeks -- go back and collect explant date
if isnat(pt_implant.Explant_Date); pt_implant.Explant_Date = pt_implant.Implant_Date + days(21); end
%%
EMU_med_inds= (meds.TAKEN_TIME <= pt_implant.Explant_Date) & meds.TAKEN_TIME >= pt_implant.Implant_Date;
EMU_meds = meds(EMU_med_inds,:);

% filter AEDs only 
emu_med_names = cellfun(@lower,EMU_meds.MEDICATION_DISPLAY_NAME,'UniformOutput',false);
emu_med_names = cellfun(@strsplit,emu_med_names,'UniformOutput',false);
emu_med_names = cellfun(@(x)x{1},emu_med_names,'UniformOutput',false);

if strcmp(ptID,'HUP152')
    nonform = contains(emu_med_names,'non-formulary');
    emu_med_names(nonform)={'clorazepate'};
end 


[gen_names] = get_AEDs_from_list(emu_med_names);
del_inds = cell2mat(cellfun(@isempty,gen_names,'UniformOutput',false));
EMU_meds(del_inds,:)=[];
EMU_meds.generic_name = gen_names(~del_inds);

% filter medications that were not given 
given = contains(EMU_meds.MAR_ACT,'Given') | contains(EMU_meds.MAR_ACT,'New Bag') |  contains(EMU_meds.MAR_ACT,'Administered by Provider') ;
EMU_meds = EMU_meds(given,:);

% now get field values for meds table
meds=table();

meds.HUP = EMU_meds.HUP_ID;
meds.medication = EMU_meds.generic_name;
meds.route= EMU_meds.ADMINISTRATION_ROUTE;
meds.dose = EMU_meds.DOSE;
meds.dose_unit = EMU_meds.DOSE_UNIT;
meds.date = EMU_meds.TAKEN_TIME;
meds.ActionTime = timeofday(EMU_meds.TAKEN_TIME);

% get the hr of the total EMU stay for the admin
[~,sorted_inds] = sort(meds.date);
meds = meds(sorted_inds,:);
start_day = meds.date(1)-meds.ActionTime(1);
admin_times = hours(meds.date-start_day);
meds.admin_time = admin_times;

%% now get ieeg times for medication 
ieeg_info = readtable('HUP_ieeg_conversion.xlsx','VariableNamingRule' , 'preserve');
ieeg_info.dataset_start_EMU = datetime(ieeg_info.dataset_start_EMU);
pt_inds = ieeg_info.HUP == str2double(ptID(4:end));
pt_ieeg = ieeg_info(pt_inds,:); 

%temp dates of dataset start
t = pt_ieeg.dataset_start_ieeg;
dataset_starts_ieeg = t-timeofday(t); %ieeg time
day_offsets = cumsum([0; days(diff(dataset_starts_ieeg))]);
for i = 1:height(pt_ieeg)
    if isnat(pt_ieeg.dataset_start_EMU(1))
        pt_ieeg.dataset_start_EMU(i) = pt_ieeg.EMU_admission_date(i)+days(day_offsets(i)) + timeofday(pt_ieeg.dataset_start_ieeg(i));
    end
end

each_dataset_start =NaT(height(meds),1);
each_dataset_start_ieeg=NaT(height(meds),1);
each_dataset_name = cell(height(meds),1);

dataset_starts= pt_ieeg.dataset_start_EMU; %emu time
dataset_names = unique(pt_ieeg.dataset);
for n=1:length(dataset_starts)
    each_dataset_start_ieeg(meds.date > dataset_starts(n))=dataset_starts_ieeg(n);
    each_dataset_start(meds.date > dataset_starts(n))=dataset_starts(n);
    each_dataset_name(meds.date > dataset_starts(n)) = dataset_names(n);
end
% for the last entries, assume it is part of the last dataset; later if error occurs then admin is after end of ieeg recording
meds.OffsetSecondsInIeeg = seconds(meds.date-each_dataset_start);
meds.dataset=each_dataset_name; 
meds.dataset(cellfun(@isempty,meds.dataset))={'before start of ieeg recording'};
meds.DatasetStartIeeg = each_dataset_start_ieeg;

% other outputs
med_names = unique(meds.medication);
explant_date = pt_implant.Explant_Date;
end 
