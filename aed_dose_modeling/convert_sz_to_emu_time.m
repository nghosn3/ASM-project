function [seizure_times] = convert_sz_to_emu_time(offsets,ieeg_offset_datasets,ptID)

%get seizure times
[seizure_times,seizure_dataset] = get_seizure_times_from_sheet(ptID);
if ~isempty(offsets)
    for j =1:height(seizure_times)
        % check which dataset the seizure is from, and add appropriate offset
        if isequal(seizure_dataset{j},'D01') || isequal(seizure_dataset{j},'one file')
            seizure_times(j,1)= (offsets(1)+(seizure_times(j,1)))./3600;
        else
            %ind = str2double(seizure_dataset{j}(end));
            ind = contains(ieeg_offset_datasets,['D0' seizure_dataset{j}(end)]);
            dataset_offset = offsets(ind);
            seizure_times(j,1)= (seizure_times(j,1) + dataset_offset)./3600; %convert to hours
        end
    end
end

end