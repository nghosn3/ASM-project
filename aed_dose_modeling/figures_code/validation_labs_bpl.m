%% check correlations between measured lab values and predicted bpl values

close all;clear;

% Get which patients
cohort_info = readtable('HUP_implant_dates.xlsx');
ptIDs = cohort_info.ptID;
weights = cohort_info.weight_kg;

% get AED dose curves
[all_dose_curves,all_tHr,ptIDs,all_med_names,ieeg_offset,max_dur,emu_dur] = get_aed_curve_kg(ptIDs,weights); %#ok<*ASGLU>

% load AED parameters and
aed_params = readtable('AED_metadata.xlsx');
aed_params.medication = arrayfun(@lower,aed_params.medication);
aed_ref_ranges = readtable('AED_metadata.xlsx','Sheet','refRanges');

% LOAD DATA
load([pwd '/DATA/all_labs.mat'],'all_labs');
load('MAR_032122.mat');


%%
predictions = table();
for i = 1:length(ptIDs)
    
    ptID  = ['HUP' num2str(ptIDs(i))];
    curves = all_dose_curves{i};
    EMU_labs = parse_labs(ptID,all_labs);
    [med_names,meds,explant_date] = parse_MAR(ptID,all_meds);
    
    pred_levels = zeros(height(EMU_labs),1);
    % get the estimated level
    for n = 1:height(EMU_labs)
        med_ind = [];
        for r = 1:length(med_names)
            if strcmp(EMU_labs.generic_name(n),(med_names{r}))
                med_ind = r;
            end
        end
       
        if med_ind>0
            curve = curves{med_ind};
            % normalize curve?
            % curve = curve./max(curve);
            
            this_lab_time = EMU_labs.SPECIMN_TAKEN_TIME(n);
            if this_lab_time > meds.date(1)
                lab_emu_tHr = hours(this_lab_time-meds.date(1));
                lab_ind = round(lab_emu_tHr * 60);
                pred_levels(n) = curve(lab_ind);
            end
            
        end
        
    end
    % append to the table 
    EMU_labs.predicted_level = pred_levels;
    predictions =vertcat(predictions,EMU_labs);
    
end

%% remove NaN's - likely duplicates 
predictions(isnan(predictions.RESULTS),:)=[];
% remove values where predictions were 0 - level is likely due to
% pre-existing level and not from EMU ADMIN
predictions((predictions.predicted_level == 0),:)=[];


% unit conversion... gotta do it manually: convert everything to mg/L
all_units = unique(predictions.UNIT); 
%  {'mg/L' } and {'ug/mL'} are equal
%  {'ng/mL'} is 1e-3 of {'mg/L' }
ng_inds = strcmp(predictions.UNIT,'ng/mL');
predictions.RESULTS(ng_inds) = predictions.RESULTS(ng_inds) * 1e-3;
predictions.UNIT(ng_inds)= {'mg/L'};

% manual: messed up 
 predictions.generic_name(strcmp(predictions.generic_name,{'n-desmethylclobazam'})) = {'clobazam'};
% predictions(strcmp(predictions.UNIT,{'%'}),:) = [];


all_med_names = unique(predictions.generic_name);
figure(1);
for i = 1:length(all_med_names)
    subplot(3,5,i)
    inds = contains(predictions.generic_name,all_med_names(i));
    scatter(predictions.RESULTS(inds),predictions.predicted_level(inds)); axis square; hold on;
    maxval = max([predictions.RESULTS(inds); predictions.predicted_level(inds)]);
    try
        x = predictions.RESULTS(inds); y = (predictions.predicted_level(inds));
        [r,p] = corr([x y]);
        r = r(2,1);
        p = p(2,1);
        coefs = polyfit(x,y,1);
        pred = (coefs(1)*x) +coefs(2);
        %legend(['p: ' num2str(p) ', R: ' num2str(r)]);
        plot(x,pred,'k','Linewidth',2); 
    catch 
    end
    ylim([0 maxval]);
    xlim([0 maxval]);
    
    title([all_med_names{i} ', R = ' num2str(r)]);
    xlabel('measured lab result');ylabel('model predicted bpl')
end
%%
figure(2)

subplot(1,2,1)
scatter(log10(predictions.RESULTS+1),log10(predictions.predicted_level+1)); axis square; hold on;

% get linear correlation and line of best fit
x = log10(predictions.RESULTS+1); y = log10(predictions.predicted_level+1);
[r,p] = corr([x y]);
r = r(2,1);
p = p(2,1);
coefs = polyfit(x,y,1);
pred = (coefs(1)*(x) +coefs(2));
plot(x,pred,'k','Linewidth',2); hold off;

title(['all AED lab results, R = ' num2str(r) ', p =  ' num2str(p) ]); %xlim([0 60])
xlabel('log measured lab result (mg/L)');ylabel('log model predicted BPL (mg/L)')

subplot(1,2,2) 
scatter((predictions.RESULTS),(predictions.predicted_level)); axis square; hold on;
x = (predictions.RESULTS); y = (predictions.predicted_level);
[r,p] = corr([x y]);
r = r(2,1);
p = p(2,1);
coefs = polyfit(x,y,1);
pred = (coefs(1)*x) +coefs(2);
%legend(['p: ' num2str(p) ', R: ' num2str(r)]);
plot(x,pred,'k','Linewidth',2); hold off;

title(['all AED lab results, R = ' num2str(r) ', p =  ' num2str(p) ]); 
xlim([0 150]); ylim([0 150])
xlabel('measured lab result (mg/L)');ylabel('model predicted BPL (mg/L)')


% cohort 
ptIDs_labs = unique(predictions.HUP_ID);
labs_n = length(ptIDs_labs);
drugs_labs = unique(predictions.generic_name);
drugs_n = length(drugs_labs);


figure(3)
subplot(1,2,2) 
colors = turbo(length(all_med_names));
for i = 1:length(all_med_names)
    
    inds = contains(predictions.generic_name,all_med_names(i));
    plot(log10(predictions.RESULTS(inds)),log10(predictions.predicted_level(inds)),'.','markersize',20,'Color',colors(i,:));
    axis square; hold on;
end

ylim([0 2.5]);
xlim([0 2.5]);
xlabel('measured lab result');ylabel('model predicted bpl')
x = log10(predictions.RESULTS+1); y = log10(predictions.predicted_level+1);
[r,p] = corr([x y]);
r = r(2,1);
p = p(2,1);
coefs = polyfit(x,y,1);
pred = (coefs(1)*(x) +coefs(2));
plot(x,pred,'k','Linewidth',2); hold off;

title(['all AED lab results, R = ' num2str(r) ', p =  ' num2str(p) ]); %xlim([0 60])
xlabel('log measured lab result (mg/L)');ylabel('log model predicted BPL (mg/L)')
legend(all_med_names);


subplot(1,2,1) 
for i = 1:length(all_med_names)
    
    inds = contains(predictions.generic_name,all_med_names(i));
    plot((predictions.RESULTS(inds)),(predictions.predicted_level(inds)),'.','markersize',20,'Color',colors(i,:));
    axis square; hold on;
end

ylim([0 150]);
xlim([0 150]);
xlabel('measured lab result');ylabel('model predicted bpl')
x = (predictions.RESULTS+1); y = (predictions.predicted_level+1);
[r,p] = corr([x y]);
r = r(2,1);
p = p(2,1);
coefs = polyfit(x,y,1);
pred = (coefs(1)*(x) +coefs(2));
plot(x,pred,'k','Linewidth',2); hold off;

title(['all AED lab results, R = ' num2str(r) ', p =  ' num2str(p) ]); %xlim([0 60])
xlabel('measured lab result (mg/L)');ylabel('model predicted BPL (mg/L)')
legend(all_med_names);

%% bland aldman plot

aed_ref_ranges.Drug = arrayfun(@lower,aed_ref_ranges.Drug);

% convert units of ref ranges to mg/L 
%  {'mg/L' } and {'ug/mL'} are equal and {'ng/mL'} is 1e-3 of {'mg/L' }
ng_inds = strcmp(aed_ref_ranges.units,'ng/mL');
aed_ref_ranges.ref_range_min(ng_inds) = aed_ref_ranges.ref_range_min(ng_inds) * 1e-3;
aed_ref_ranges.ref_range_max(ng_inds) = aed_ref_ranges.ref_range_max(ng_inds) * 1e-3;
aed_ref_ranges.units(ng_inds)= {'mg/L'};


all_measured_z = zeros(height(predictions),1);
all_predicted_z = zeros(height(predictions),1);

for i =1:height(predictions)
    med_ind = find(strcmp(predictions.generic_name(i),aed_ref_ranges.Drug));
    range = aed_ref_ranges.ref_range_max(med_ind) -  aed_ref_ranges.ref_range_min(med_ind);
    ref_mean = mean([aed_ref_ranges.ref_range_max(med_ind) aed_ref_ranges.ref_range_min(med_ind)]);
    
    all_measured_z(i) = (predictions.RESULTS(i) - ref_mean) ./ (range);
    all_predicted_z(i) = (predictions.predicted_level(i) - ref_mean) ./ (range);
end

figure;
pred_names =unique(predictions.generic_name);
colors = turbo(length(pred_names));
for i = 1:length(pred_names)
    inds = find(contains(predictions.generic_name,pred_names(i)));
    jitter = rand(length(inds),1)./5;
    plot(zeros(length(inds),1)+i+jitter-mean(jitter),all_measured_z(inds)-all_predicted_z(inds),'.','Color',colors(i,:),'Markersize',25); hold on;
end

[B,I]=sort(predictions.generic_name);
boxplot(all_measured_z(I)-all_predicted_z(I),B)


legend(pred_names)
ylabel('measured - predicted (reference ranges)');
yline(0,'linewidth',3)
yline(-1,'linewidth',1)
yline(1,'linewidth',1)
xlim([0 14]);
xticks(1:length(pred_names)); xticklabels(pred_names);



ylim([-3 3])





