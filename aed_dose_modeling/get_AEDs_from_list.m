function [gen_names] = get_AEDs_from_list(aed_label)

% load aed list to create map
aed_list=readtable('aed_list.csv');

% convert to lower case
aed_list.AED = cellfun(@lower,aed_list.AED,'UniformOutput',false);
aed_list.Generic = cellfun(@lower,aed_list.Generic,'UniformOutput',false);
%make map with keys
aed_map = containers.Map(aed_list.AED,aed_list.Generic);

%convert med names to all generic
gen_names =cell(length(aed_label),1);
for i = 1:length(aed_label)
    if contains(aed_label{i},aed_list.Generic)
        gen_names{i} =aed_label{i};
    elseif contains(aed_label{i},aed_list.AED)
        gen_names{i} =aed_map(aed_label{i});
        
    else
        gen_names{i} = '';
        
    end
end

end
