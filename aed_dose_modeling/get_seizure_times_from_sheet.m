% function that pulls all the seizure times for a patient, and returns the
% start and end times in HOURS
function [seizure_times,seizure_dataset] = get_seizure_times_from_sheet(ptID)

% load info from borel:
% addpath('/Volumes/borel.seas.upenn.edu/public/USERS/nghosn3/Pioneer/spikes-AED')
% addpath('/Volumes/borel.seas.upenn.edu/public/USERS/nghosn3/Pioneer/DATA')
all_seizures = readtable('Erin_szs_times.xlsx','VariableNamingRule','preserve','Sheet','Erin_szs_times');

% w = warning('query','last');
% id = w.identifier;
% warning('off',id)

dataset=cellfun(@(x) strsplit(x, '_'),all_seizures.IEEGname,'UniformOutput',false);
 
dataset_name = cell(length(dataset),1);
for i =1:length(dataset)
    dataset_name(i) = dataset{i}(end);
end

dataset_name(~contains(dataset_name,'D'))={'one file'};
all_seizures.Dataset=dataset_name;

%get the seizure times

pt_inds = contains(all_seizures.Patient,ptID);
pt_seizure_info = all_seizures(pt_inds,:);

seizure_times = [pt_seizure_info.start pt_seizure_info.end];
seizure_dataset = pt_seizure_info.Dataset;

end


