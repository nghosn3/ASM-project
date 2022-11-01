curr_path = pwd;
addpath([pwd '/DATA'])

kas=[];

aed_params = readtable('AED_metadata.xlsx');
aed_params.medication = arrayfun(@lower,aed_params.medication);

tHalfs = [aed_params.t_half_e];

for i=1:height(aed_params)
    tHalf=strsplit(tHalfs{i},',');
    tHalf = cellfun(@str2double,tHalf);
    tHalf = mean(tHalf);
    tmax = aed_params.t_max(i);
    ka=get_ka(tHalf,tmax);
    kas =[kas; double(ka)];
    
end
