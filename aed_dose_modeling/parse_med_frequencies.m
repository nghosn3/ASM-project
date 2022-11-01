function [spd] = parse_med_frequencies(med_string)

spd = NaN;

%maps text numbers like one to numbers like 1
num_str = containers.Map({'one' 'two' 'three' 'four' 'five' 'six' 'seven' ...
    'eight' 'nine' 'ten', 'once', 'twice'}, {'1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '1', '2'});

%split the string into words
med_split = lower(split(med_string));
%remove the second sentence
period = find(contains(med_split,'.'));
if ~isempty(period)
    if length(period)==1
        if isnan(str2double(med_split(period)))
            med_split = med_split(1:period);
        end 
    else 
    period = period(end); %if there is a decimal in a number
    med_split = med_split(1:period);
    end 
end

%keep track of where the times and numbers are
num_inds = zeros(1,length(med_split));

for i = 1:length(med_split)
    
    %check if the number for freq is represented in text form and replace as needed
    if isKey(num_str, lower(med_split{i}))
        med_split{i} = num_str(lower(med_split{i}));
    end
    
    % 1 represents a number (dose or frequency) and 2 represents a time phrase
    if ~isnan(str2double(med_split{i}))
        num_inds(i) = 1;
    elseif contains(lower(med_split{i}), 'hour') || contains(lower(med_split{i}), 'day') || contains(lower(med_split{i}), 'dai') || contains(lower(med_split{i}), 'morning') || contains(lower(med_split{i}), 'am') || contains(lower(med_split{i}), 'pm') || contains(lower(med_split{i}), 'mid') || contains(lower(med_split{i}), 'evening') || contains(lower(med_split{i}), 'bed') || contains(lower(med_split{i}), 'night')
        num_inds(i) = 2;
    end
end

%try to parse the frequency spd
%find where times can occur
time_idx = find(num_inds==2);


%if a time wasn't found, then assume tablets/capsules in a day or a dose
%(mg) are listed
if ~isempty(time_idx)
    %search for a XX tab phrase
    
    tab_idx = find(contains(med_split, 'tab')| contains(med_split, 'cap') | contains(med_split, 'mg') | contains(med_split, 'po'));
    time_idx = num_inds == 1; time_idx(tab_idx-1)=0; % numbers that are not tabs or doses
    %for each XX tab phrase, check that it is preceeded by a number
    total_freq = 0;
    
    if sum(time_idx)==1
        total_freq = str2double(med_split(time_idx));
        spd = total_freq;
    elseif ~isempty(tab_idx) && sum(time_idx)~=1
        for i = 1:length(tab_idx)
            if (num_inds(tab_idx(i)-1)==1)
                total_freq = total_freq + 1;
            end
        end
        if total_freq ~= 0
            spd = total_freq;
        end
        
    end
end 
    %else, if a number preceeds the first hour time phrase, then calculate the frequency
    for i=1:length(med_split)
        if contains(lower(med_split{i}), 'hour')
            %disp('contains hour')
            hrs = str2double(med_split{i-1});
            spd = 24 / hrs;
        end
    end
    
    

isneeded = 0;
if (contains(lower(med_string), 'need') || contains(lower(med_string), 'seiz') || contains(lower(med_string), 'sz'))
    spd = NaN;
    isneeded =1;
end

if isnan(spd) && ~isneeded
    disp([newline 'Unable to parse spd for SIG ' med_string{:}])
end



end


