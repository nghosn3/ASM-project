%% find correlation between AED load and seizure severity
% look at seizure severity and tapering of individual drugs
% look at seizure severity and overall AED load

close all;clear;
curr_path = pwd;
addpath(curr_path)
addpath([curr_path '/DATA'])
addpath([curr_path '/figures_code'])

% Get which patients have AED data, load the data
load('MAR_032122.mat')
cohort_info = readtable('HUP_implant_dates.xlsx');
ptIDs = cohort_info.ptID;
weights = cohort_info.weight_kg;


% only choose leading seizure in a cluster
sz_cluster =0;
% only use tapered drugs to find drug load 
tapered=0;
exclude_ativan = 0;


% Get the patients that also have seizure severity scores
seizure_metadata = readtable('seizure_metadata_with_severity.xlsx', "ReadVariableNames" ,true);
sz_pts = unique(seizure_metadata.Patient);

% manual exclusion of patients that didnt have taper
% exclude_pts = [171 172 138 150 188 157];
% exclude_inds = zeros(1,length(sz_pts));
% for i = 1:length(exclude_pts)
%     pt = ['HUP' num2str(exclude_pts(i))];
%     pt_ind = contains(sz_pts,pt);
%     exclude_inds = exclude_inds + pt_ind';
% end
% sz_pts = sz_pts(~logical(exclude_inds));

% patients included in analysis
pt_inds=zeros(1,length(ptIDs));
for x = 1:length(ptIDs)
    ptID = ['HUP' num2str(ptIDs(x))];
    pt_inds(x) = sum(contains(sz_pts,ptID));
end

ptIDs = ptIDs(logical(pt_inds));
weights = weights(logical(pt_inds));



if tapered
    [tapered_pts,all_taper_info] = get_tapered_patients(ptIDs,all_meds);
    ptIDs = ptIDs(logical(tapered_pts));
     weights = weights(logical(tapered_pts));
end 

[all_dose_curves,all_tHr,ptIDs,all_med_names,ieeg_offset,max_dur,emu_dur] = get_aed_curve_kg(ptIDs,weights);
%[all_dose_curves,all_tHr,all_med_names,ieeg_offset,emu_dur] = get_dose_schedules(ptIDs);

%%
%all seizures within a pt
preictal_AEDs_1hr=[];
all_pt_inds=[];
sz_pt = [];
figure;
for x=1:length(ptIDs)
    ptID = ['HUP' num2str(ptIDs(x))];
    pt_inds = find(contains(seizure_metadata.Patient,ptID));
    pt_szData = seizure_metadata(pt_inds,:);
    seizure_times = pt_szData.SeizureEEC; %in seconds
    
   
    [med_names,meds,explant_date] = parse_MAR(ptID,all_meds);
   
    % get the datset name for the offset
    dataset=cellfun(@(x) strsplit(x, '_'),pt_szData.iEEGFilename,'UniformOutput',false);
    dataset_name = cell(length(dataset),1);
    for i =1:length(dataset)
        dataset_name(i) = dataset{i}(end);
    end
    dataset_name(~contains(dataset_name,'D'))={'one file'};
    seizure_dataset = dataset_name;
    
    % get the total AED dose over time
    pt_curves = all_dose_curves{x};
    tHr = all_tHr{x};
    
    % exclude ativan in drug sum
    if exclude_ativan
        ativan_ind = contains(med_names,'lorazepam');
        med_names(ativan_ind)=[];
        pt_curves(ativan_ind)=[];
        tHr(ativan_ind)=[];
    end
    drugs =zeros(length(med_names),ceil(emu_dur(x)*60)); 
    for n =1:length(med_names)
        drug=pt_curves{n};
        drug=drug./nanmax(drug); %normalize each drug curve
        if ~isempty(drug)
            dStart = round(tHr{n}(1)*60)-1;
            drugs(n,dStart+1:dStart+length(drug))=drug;
        end
    end
    drug_sum=nansum(drugs,1);
    drug_sum=drug_sum./length(med_names);
    % cut off drug curve to only be length of emu stay
    %drug_sum(find(drug_sum(1:ceil(emu_dur(x)*60))==0))=NaN; %delete rather than NaN it 
    time = (1:length(drug_sum))./60;

    offsets = ieeg_offset{2,x};
    ieeg_offset_datasets=ieeg_offset{1,x};
    for j =1:length(seizure_times)
        % check which dataset the seizure is from, and add appropriate offset
        if isequal(seizure_dataset{j},'D01') || isequal(seizure_dataset{j},'one file')
            seizure_times(j,1)= (offsets(1)+(seizure_times(j,1)))./3600; %hours
        else
            ind = contains(ieeg_offset_datasets,['D0' seizure_dataset{j}(end)]);
            dataset_offset = offsets(ind);
            seizure_times(j,1)= (seizure_times(j,1) + dataset_offset)./3600; %hours
        end
        
    end
    
    % convert seizure times to indices of drug curve, so to minutes
    seizure_inds = round(seizure_times(:,1) *60);
    
    % plot the drug sum and the seizures
    subplot(6,4,x) 
    plot(time,drug_sum); hold on;
    plot(time(seizure_inds),drug_sum(seizure_inds),'o'); title(ptID); hold off;
    
    % here we want pick first seizure in a cluster defined by cutoff
    if sz_cluster ==1
        cutoff = 60;
        dSz_times =[cutoff+1; diff(seizure_inds)];
        leadSz_inds = dSz_times > cutoff; %seizures more than one hour apart
        seizure_inds =seizure_inds(leadSz_inds);
        %pt_SS = pt_SS(leadSz_inds);
        pt_inds = pt_inds(leadSz_inds);
    end
    all_pt_inds=[all_pt_inds pt_inds'];
     
    % med curves are sampled at one point per minute
    pt_preictal_1hr = zeros(1,length(seizure_inds));
    for i=1:length(seizure_inds)
        pt_preictal_1hr(i) = nanmean(drug_sum(seizure_inds(i)-60:seizure_inds(i)));
    end
    
    preictal_AEDs_1hr=[preictal_AEDs_1hr pt_preictal_1hr];
    sz_pt = [sz_pt ptIDs(x)];
end

%% correlation between AED load and seizure severity 
all_data = seizure_metadata(all_pt_inds,:);
all_data.preictal_AEDs_1hr = preictal_AEDs_1hr';

pt_szTypes = seizure_metadata.SeizureType(all_pt_inds);
severity_scores = seizure_metadata.SeizureSeverity(all_pt_inds);

% remove nan inds
naninds = isnan(preictal_AEDs_1hr);

%save('severity_scores_aeds.mat','preictal_AEDs_1hr','pt_szTypes','severity_scores','sz_pt')
save('severity_scores_aeds.mat','all_data')

[r,p] = corr([preictal_AEDs_1hr(~naninds)' severity_scores(~naninds)]);
r = r(2,1);
p = p(2,1);

figure();
%linear correlation for severity scores (duration*convexHull)

% seperate colors by seizure type
szTypes = unique(pt_szTypes);
%szTypes= {'FBTCS','FAS','FIAS','Subclinical'}; % corrected order according to plot data

type_colors = lines(length(szTypes)); % [red, green, blue, purple]
colors = cell(1,length(all_pt_inds));
for i=1:length(szTypes)
    inds = contains(pt_szTypes,szTypes{i});
    colors(inds)={type_colors(i,:)};    
end

title('Total AED load, all drugs','fontsize',16)
for i = 1:length(all_pt_inds)
    plot(severity_scores(i),preictal_AEDs_1hr(i),'.','Linewidth',2,'Color',colors{i},'markersize',20);hold on;
end
ylabel('normalized AED load','fontsize',14)
xlabel('seizure severity score','fontsize',14)
%linear fit 
coefs = polyfit(severity_scores(~naninds)',preictal_AEDs_1hr(~naninds)',1);
aed_pred = (coefs(1)*severity_scores(~naninds)) +coefs(2);
plot(severity_scores(~naninds),aed_pred,'k','Linewidth',2);
legend(['p: ' num2str(p) ', R: ' num2str(r)]);

% convex hull - t-score them for correlation
hull_vol = seizure_metadata.TotalVolume_cm_3_(all_pt_inds);
mean_hull = mean(hull_vol);
std_hull=std(hull_vol);
hull_tscores = (hull_vol-mean_hull)./(std_hull./sqrt(length(hull_vol)));

[r,p] = corr([preictal_AEDs_1hr' hull_vol]);
r = r(2,1);
p = p(2,1);

figure()
subplot(1,2,1)
title('Total AED load, all drugs','fontsize',16);
for i = 1:length(all_pt_inds)
    plot(hull_vol(i),preictal_AEDs_1hr(i),'.','Linewidth',2,'Color',colors{i},'MarkerSize',10);hold on;
end
ylabel('normalized AED load','fontsize',14)
xlabel('Convex hull (vol)','fontsize',14)
axis square;
%linear fit
coefs = polyfit(hull_vol',preictal_AEDs_1hr',1);
aed_pred = (coefs(1)*hull_vol) +coefs(2);
legend(['p: ' num2str(p) ', R: ' num2str(r)]);
plot(hull_vol,aed_pred,'k','Linewidth',2);
% duration
durations = seizure_metadata.SeizureDuration(all_pt_inds);
mean_dur = mean(durations);
std_dur=std(durations);
dur_tscores = (durations-mean_dur)./(std_dur./sqrt(length(durations)));

[r,p] = corr([preictal_AEDs_1hr' durations]);
r = r(2,1);
p = p(2,1);

subplot(1,2,2)
title('Total AED load, all drugs','fontsize',16)
for i = 1:length(all_pt_inds)
    plot(durations(i),preictal_AEDs_1hr(i),'.','Linewidth',2,'Color',colors{i},'MarkerSize',10);hold on;
end
ylabel('normalized AED load','fontsize',14)
xlabel('duration (s)','fontsize',14);
axis square;
%linear fit
coefs = polyfit(durations',preictal_AEDs_1hr',1);
aed_pred = (coefs(1)*durations) +coefs(2);
legend(['p: ' num2str(p) ', R: ' num2str(r)]);
plot(durations,aed_pred,'k','Linewidth',2);

%% make boxplots for the seizure types and do distribution tests

fas_inds = strcmp(pt_szTypes,'FAS');
fbtc_inds = strcmp(pt_szTypes,'FBTCS');
fias_inds = strcmp(pt_szTypes,'FIAS');
sub_inds = strcmp(pt_szTypes,'Subclinical');


fas_aeds = preictal_AEDs_1hr(fas_inds);
fbtc_aeds = preictal_AEDs_1hr(fbtc_inds);
fias_aeds = preictal_AEDs_1hr(fias_inds);
sub_aeds = preictal_AEDs_1hr(sub_inds);

% get durations
fas_dur = all_data.SeizureDuration(fas_inds);
fias_dur = all_data.SeizureDuration(fias_inds);
fbtc_dur = all_data.SeizureDuration(fbtc_inds);
sub_dur = all_data.SeizureDuration(sub_inds);


r = max([length(fbtc_dur) length(fas_dur) length(fias_dur) length(sub_dur) ]); c=length(szTypes);
x = nan(r,c);
x(1:length(fbtc_dur),1) = fbtc_dur;
x(1:length(fas_dur),2) = fas_dur;
x(1:length(fias_dur),3) = fias_dur;
x(1:length(sub_dur),4) = sub_dur;

% do the stats for duration
[p,tbl,stats] = kruskalwallis(x,szTypes,'on');

multcompare(stats)

data=[mean(fbtc_dur) mean(fas_dur) mean(fias_dur) mean(sub_dur)];
stds = [std(fbtc_dur) std(fas_dur) std(fias_dur) std(sub_dur)];
% 
% b=bar(data); hold on;

% %b.FaceColor = 'flat';
% b.CData =type_colors;
% % plot data points on top 
% plot(.9,fbtc_dur,'.k','markersize',10);
% plot(1.9,fas_dur,'.k','markersize',10);
% plot(2.9,fias_dur,'.k','markersize',10);
% plot(3.9,sub_dur,'.k','markersize',10);
% 
% er = errorbar(1:length(szTypes),data,stds,'LineWidth',2);    
% er.Color = [0 0 0];                            
% er.LineStyle = 'none'; 
% 
% a = get(gca,'XTickLabel');  
% set(gca,'XTickLabel',a,'fontsize',14,'FontWeight','bold')
% set(gca,'XTickLabelMode','auto');
% xticklabels(szTypes);
% ylabel('duration (s))')
% title('sz duration of clinical seizure types','fontsize',18)

% for fbtcs vs all other, duration
figure;
data=[mean(fbtc_dur) mean([fas_dur; fias_dur; sub_dur])];
stds = [std(fbtc_dur) std([fas_dur; fias_dur; sub_dur])];

b=bar(data); hold on;
%b.FaceColor = 'flat';
b.CData =type_colors;
% plot data points on top 
plot(.9,fbtc_dur,'.k','markersize',10);
plot(1.9,[fas_dur; fias_dur; sub_dur],'.k','markersize',10);


er = errorbar([1 2],data,stds,'LineWidth',2);    
er.Color = [0 0 0];                            
er.LineStyle = 'none'; 

a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',14,'FontWeight','bold')
set(gca,'XTickLabelMode','auto');
xticklabels({'convulsive seizure','other'});
ylabel('duration (sec)')
title('duration of clinical seizure types','fontsize',18)

r = max([length(fbtc_dur) length([fas_dur; fias_dur; sub_dur]) ]); c=2;
x = nan(r,c);
x(1:length(fbtc_dur),1) = fbtc_dur;
x(1:length([fas_dur; fias_dur; sub_dur]),2) = [fas_dur; fias_dur; sub_dur];

[p,tbl,stats] = kruskalwallis(x,{'fbtcs','other'},'on');



%% aed stats
r = max([length(fbtc_aeds) length(fas_aeds) length(fias_aeds) length(sub_aeds) ]); c=length(szTypes);
x = nan(r,c);
x(1:length(fbtc_aeds),1) = fbtc_aeds;
x(1:length(fas_aeds),2) = fas_aeds;
x(1:length(fias_aeds),3) = fias_aeds;
x(1:length(sub_aeds),4) = sub_aeds;


% do the stats 
[p,tbl,stats] = kruskalwallis(x,szTypes,'on');
multcompare(stats)

% plot stuff
figure;

data=[mean(fbtc_aeds) mean(fas_aeds) mean(fias_aeds) mean(sub_aeds)];
stds = [std(fbtc_aeds) std(fas_aeds) std(fias_aeds) std(sub_aeds)];

b=bar(data); hold on;
%b.FaceColor = 'flat';
b.CData =type_colors;
% plot data points on top 
plot(.9,fbtc_aeds,'.k','markersize',10);
plot(1.9,fas_aeds,'.k','markersize',10);
plot(2.9,fias_aeds,'.k','markersize',10);
plot(3.9,sub_aeds,'.k','markersize',10);


% 
% er = errorbar(1:length(szTypes),data,stds,'LineWidth',2);    
% er.Color = [0 0 0];                            
% er.LineStyle = 'none'; 

a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',14,'FontWeight','bold')
set(gca,'XTickLabelMode','auto');
xticklabels(szTypes);
ylabel('AED load (normalized)')
title('AED load in clinical seizure types','fontsize',18)

% figure for convulsive sz vs. not 
figure;

data=[mean(fbtc_aeds) mean([fas_aeds fias_aeds sub_aeds])];
stds = [std(fbtc_aeds) std([fas_aeds fias_aeds sub_aeds])];

b=bar(data); hold on;
%b.FaceColor = 'flat';
b.CData =type_colors;
% plot data points on top 
plot(.9,fbtc_aeds,'.k','markersize',10);
plot(1.9,[fas_aeds fias_aeds sub_aeds],'.k','markersize',10);


er = errorbar([1 2],data,stds,'LineWidth',2);    
er.Color = [0 0 0];                            
er.LineStyle = 'none'; 

a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',14,'FontWeight','bold')
set(gca,'XTickLabelMode','auto');
xticklabels({'convulsive seizure','other'});
ylabel('AED load (normalized)')
title('AED load in clinical seizure types','fontsize',18)

r = max([length(fbtc_aeds) length([fas_aeds fias_aeds sub_aeds]) ]); c=2;
x = nan(r,c);
x(1:length(fbtc_aeds),1) = fbtc_aeds;
x(1:length([fas_aeds fias_aeds sub_aeds]),2) = [fas_aeds fias_aeds sub_aeds]';

[p,tbl,stats] = kruskalwallis(x,{'fbtcs','other'},'on');
multcompare(stats)

%% run a mixed effects model on the seizures where the predictor variable is the AED load and the response is convulsion or not. random effect of patient number

% feats_train = table();
% 
% sz_labels = all_data.SeizureCategory;
% conv = contains(sz_labels,'FBTCS');
% 
% patients = all_data.Patient;
% grab_ptID = @(x) str2double(x(4:end));
% pts = cellfun(grab_ptID,patients);
% 
% feats_train.convulsion = double(conv); % 1 for convulsion, 0 for nonconvulsion
% feats_train.ptID = pts;
% feats_train.aed_load = all_data.preictal_AEDs_1hr;
% feats_train.seizure_severity = all_data.SeizureSeverity;
% 
% % binary response of convulsive seizure
% modelspec = 'convulsion ~ 1 + aed_load  + (1 | ptID)';
% mdl = fitglme(feats_train,modelspec,'Distribution','binomial');

% continuous response of seizure severity
% modelspec = 'seizure_severity ~ 1 + aed_load  + (1 | ptID)';
% mdl = fitlme(feats_train,modelspec);
% 
% % which patients
% pts = unique(Patients)
% for i = 1:length()



