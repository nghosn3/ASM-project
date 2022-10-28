%% function to get the taper speed for each drug by fitting the decay to a logistic curve

%%outputs: 
% rate: the rate of exponential decay, representing how fast a drug was tapered based on the chnage in blood levels
% taper_days: which days the drug was tapered on - can get duration from this if needed
% tapered_med_names: drug names, without ativan
% start_end_dose: the dose they came in on, and min dose they reached in taper

function [rates,taper_period,tapered_med_names,start_end_doses,day1_change] = get_taper_speed(pt_curves,tHr,meds,med_names,implant_date)

% exclude ativan from taper - remove from pt_curves, tHr, and med_names
ativan_ind = contains(med_names,'lorazepam');
med_names(ativan_ind)=[];
pt_curves(ativan_ind)=[];
tHr(ativan_ind)=[];

% first find the days of the taper period by finding the daily dose of each
% drug - first day of teh first drug tapered to the last day of any drug (latest day)

% find the taper period (days)
taper_period = nan(length(med_names),2);
start_end_doses = nan(length(med_names),2);
day1_change=nan(1,length(med_names));
rates = nan(length(med_names)+1,2); % this is the rate constant from the exponential fit curve.
for n=1:length(med_names)
    med_name_ind = contains(med_names,med_names{n}); % which med
    this_med_info = meds(contains(meds.medication,med_names(med_name_ind)),:); % get info of just that med
    this_med_info(isnan(this_med_info.dose),:)=[];
    
    if  ~isempty(this_med_info) && days(this_med_info.date(end)-this_med_info.date(1)) >0 % the drug was taken for more than a day
        
        dates = this_med_info.date-this_med_info.ActionTime;
        x = (days(dates-implant_date))+1;
        this_med_info.day = x;
        y = this_med_info.dose;
        dose_per_day=zeros(1,x(end));
        for i=1:height(this_med_info)
            day = this_med_info.day(i);
            dose_per_day(day)=sum(y(x==day));
        end
        
        if dose_per_day(1)<=dose_per_day(2) && all(dose_per_day(1)<=dose_per_day) % they were not tapered
            initial_dose = dose_per_day(2); % their dose on their first full day of the EMU
            start_end_doses(n,:)=[initial_dose initial_dose]; % daily dose
            day1_change(n)=0;
            taper_period(n,:) = [0 0];
        else
             
            initial_dose =max([dose_per_day(1) dose_per_day(2)]); % if the second day is higher, indicates that the first day was a short EMU day
            day1_change(n) = 1- (dose_per_day(2)./initial_dose); % percent change within the first 24hrs
            [min_dose_achieved,min_ind] = min(dose_per_day);
            start_end_doses(n,:)=[initial_dose min_dose_achieved]; % daily dose
            end_taper = hours(hours(this_med_info.admin_time(1)+24)+days(min_ind));
            taper_period(n,:) = [this_med_info.admin_time(1)+24 end_taper];
        end
        % now fit the period of taper to an exponential curve, and get the rate constant to represent the agressiveness of taper
        if any(taper_period(n,:)==0)
            rates(n)=0;
        else 
        curve = pt_curves{n};
        time = tHr{n}; % time is in hrs
        
        [~,taper_start_ind] = min(abs(time-taper_period(n,1)));
        [~,taper_end_ind] = min(abs(time-taper_period(n,2)));
        
        % get the curve between that time 
        tapered_curve = curve(taper_start_ind:taper_end_ind);
        tapered_time = time(taper_start_ind:taper_end_ind);
        % fit curve to exponential decay 
        f = fit(tapered_time',tapered_curve','exp1');
        plot(f,tapered_time,tapered_curve);
        rates(n,1)=f.b;
        rates(n,2)=f.a;
        end 
        
    end
    
end

tapered_med_names=med_names;
end