function [tapered_pts,all_taper_info] = get_tapered_patients(ptIDs,all_meds)
tapered_pts = zeros(1,length(ptIDs));
all_taper_info=cell(1,length(ptIDs));


for i=1:length(ptIDs)
    ptID = ['HUP' num2str(ptIDs(i))];
    [med_names,meds,~] = parse_MAR(ptID,all_meds);
    try
    [taper_info,any_taper] = get_taper_info(meds, med_names);
    catch 
        disp(ptID) % what patient is erroring
    end 
    tapered_pts(i)=any_taper;
    all_taper_info{i}=taper_info;
end 

end