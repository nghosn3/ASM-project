close all;clear;

curr_path = pwd;
addpath([curr_path '/DATA'])
addpath([curr_path '/figures_code'])
tic
% Get which patients

cohort_info = readtable('HUP_implant_dates.xlsx');
ptIDs = cohort_info.ptID;
weights = cohort_info.weight_kg;

% select patient to plot curves
ptID = 179;
weight = weights(ptIDs == ptID);



[all_dose_curves,all_Hr,ptIDs,all_med_names,ieeg_offset] = get_aed_curve_kg(ptID,weight); 
%%

for ipt = 1:length(ptIDs)
    % patient medication administration
    ptID = ['HUP' num2str(ptIDs(ipt))];
    offsets = ieeg_offset{2,ipt};
    ieeg_offset_datasets = ieeg_offset{1,ipt};
    
    % plot dose curve
    % subplot(m,n,ipt)
    figure('Position', [10 10 900 300])
    %plot(0,0);
    h=zeros(1,length(all_med_names{ipt}));
    med_colors = lines(length(all_med_names{ipt}));
    for i=1:length(all_med_names{ipt})
        curve=all_dose_curves{ipt}{i};
        if ~isempty(curve)
            h(i)=plot(all_Hr{ipt}{i},curve./max(curve),'LineWidth',2,'Color',[med_colors(i,:) .5]);hold on;
            %h(i)=plot(all_Hr{ipt}{i},curve./max(curve));hold on; %normalized to [0 1]
        else
            all_med_names{ipt}(i)=[];
            h(i)=NaN;
        end
    end
    h(isnan(h))=[];
    % plot the seizures -- in ieeg times
    
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
            xline(seizure_times(j,1),'--r','linewidth',2);hold on;
        end
    end
    % convert seizure times to indices, so to minutes
    seizure_inds = round(seizure_times(:,1) *60);
    
    
    ylabel('AED BPL (normalized)');xlabel('time (Hr)')
    %xlim([0 time(end)]);
    if ~isempty(seizure_times)
        legend(h(:),all_med_names{ipt}');
    else
        legend(all_med_names{ipt}');
    end
    title(ptID);
    
end
toc

