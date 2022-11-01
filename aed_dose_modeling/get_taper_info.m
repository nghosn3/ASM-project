%% Function to determine which drugs were tapered, when, and how much the drugs were reduced from original dose

%% inputs:
% all_drug_curve: the bpl curves for each drug a patient is on
% tHR: companion time vector for the drug curve.
% meds: structure with all the medication administration information for
% the patient
%
%% Outputs:
% taper_info: struct with information about drugs tapered
% any_taper: logical representing if the patients meds were tapered or not
% Fields:
% taper_start: when this taper period started
% taper_end: when this taper period ended
% drug: the name of the drug that was tapered
% total_reduction: total reduction in dose from start of taper to end of taper (%)

%%
function [taper_info,any_taper] = get_taper_info(meds, med_names)
%this function takes in the total AED load curve and determines the time
%points of taper (day taper starts)
taper_info = struct();
%ativan admins dont count for taper - rescue drugs
taper_ind = 1;

start_day=meds.date(1)-timeofday(meds.date(1));
for n = 1:length(med_names)
    if ~contains(med_names{n},'lorazepam')
        %find the specific medication
        med_name_ind = contains(med_names,med_names{n}); % which med
        this_med_info = meds(contains(meds.medication,med_names(med_name_ind)),:); % get info of just that med
        
        % remove half day at beginning (day of admission)
        rm_days = ((this_med_info.admin_time-24) ./24) < 0;
        this_med_info(rm_days,:)=[];
        
        if height(this_med_info.date) >0 && days(this_med_info.date(end)-this_med_info.date(1)) >0
            
            x = days(this_med_info.date-start_day);
            %x=ceil(x(end:-1:1)) +1;
            x=ceil(x); %x= x-x(1) +1;
            this_med_info.day = x;
            y = this_med_info.dose;
            dose_per_day=zeros(1,x(end));
            min_dose_per = zeros(1,x(end)); % get this by dividing dose_per_day by numberof admins on that day
            for i=1:height(this_med_info)
                day = this_med_info.day(i);
                dose_per_day(day)=nansum(y(x==day));
                min_dose_per(day) = nansum(y(x==day))./ nansum(x==day);
            end
            
            
            taper_info(n).dose_schedule = {dose_per_day};
            taper_info(n).min_dose_schedule = {min_dose_per};
            taper_info(n).days = {1:length(dose_per_day)};
            taper_info(n).med_name = med_names{n};
            
            % need length of time of taper, use the ind of the dose decrease to get
            % the date, subtract the start and end dates, or just use those to
            % index
            daily_doses=dose_per_day(2:end);
            dose_change = diff(daily_doses); % exclude first day (admission)
            no_dose = daily_doses ==0;
            dose_change(no_dose)=-1;
            if ~isempty(dose_change)
%                 sign_dose_change = sign(dose_change);
%                 diff_sign_dose = diff(sign(sign_dose_change));
%                 if ~isempty(diff_sign_dose)
%                     diff_sign_dose =[diff_sign_dose(1) diff_sign_dose];
%                 end
                
                if  sum(diff(find(sign(dose_change)==-1))==1)>=2 %if the drug was tapered, and for 2 days or more
                    
                    taper_info(taper_ind).tapered =1;
                    %taper_info(taper_ind).med_name = med_names{n};
                    taper_ind = taper_ind+1;
                end
            end
        end
    end
    
    
    %taper_start = when the taper began; should be based on a decrease of dose
    %taper_end = when the tapeer of the drug ends
    % meds_tapered
    
    
    % need to get with multiple taper events?
    % ddt = diff(total_dose)'; %find where curve is decreasing
    % decreasing_inds = ddt<0;
    % taper_days = days(decreasing_inds);
end

if isfield(taper_info,'tapered')
    
    any_taper = 1;
else
    any_taper=0;
    
end
end
