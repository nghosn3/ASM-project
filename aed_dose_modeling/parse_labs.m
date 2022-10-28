%% parse lab data to find AED level labs from EMU stay
function [EMU_labs] = parse_labs(ptID,all_labs)

% load the lab data
addpath('/Volumes/borel.seas.upenn.edu/public/USERS/nghosn3/Pioneer/DATA');
addpath('/Volumes/borel.seas.upenn.edu/public/USERS/nghosn3/Pioneer/spikes-AED/aed_dose_modeling')

% get ALL patient medication information from master spreadsheet
%all_meds = readtable('all_inpatient_meds.xlsx','VariableNamingRule', 'preserve');
pt_inds = all_labs.HUP_ID== str2double(ptID(4:end));
labs=all_labs(pt_inds,:);

% load implant dates for all patients and find info for ptID
implant_info=readtable('HUP_implant_dates.xlsx','VariableNamingRule','preserve');
implant_info.implant_time=implant_info.implant_time*24;

pt_implant = implant_info(implant_info.ptID == str2double(ptID(4:end)),:);
%if explant date is missing, assume they stayed for two weeks -- go back and collect explant date
if isnat(pt_implant.Explant_Date); pt_implant.Explant_Date = pt_implant.Implant_Date + days(21); end
%%
EMU_lab_inds= (labs.SPECIMN_TAKEN_TIME <= pt_implant.Explant_Date) & labs.SPECIMN_TAKEN_TIME >= pt_implant.Implant_Date;
EMU_labs = labs(EMU_lab_inds,:);

% convert the result name to generic AED names if the test is for AED

% filter AEDs only 
emu_test_names = cellfun(@lower,EMU_labs.RESULT_NAME,'UniformOutput',false);
emu_test_names = cellfun(@strsplit,emu_test_names,'UniformOutput',false);
emu_test_names = cellfun(@(x)x{1},emu_test_names,'UniformOutput',false);

[gen_names] = get_AEDs_from_list(emu_test_names);
del_inds = cell2mat(cellfun(@isempty,gen_names,'UniformOutput',false));
EMU_labs(del_inds,:)=[];
EMU_labs.generic_name = gen_names(~del_inds);
EMU_labs.aed_name = gen_names(~del_inds);

end 