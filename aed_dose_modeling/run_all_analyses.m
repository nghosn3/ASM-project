% RUN ALL RESULTS FOR ASM DOSE MODELING AND AND ANALYSES
close all;clear;
%% run cohort statistics and get results for taper strategy used in the epilpesy monitoring unit

%starts in /aed_dose_modeling
curr_path = pwd;
addpath(curr_path)
addpath([curr_path '/DATA'])
addpath([curr_path '/figures_code'])

run cohort_stats.m
run fig_01A_and_sz_rank.m
run aed_percent_decrease_fig.m

%% run validation of blood plasma estimates using lab measured values 
 run validation_labs_bpl.m
 
 %% run mixed model analyses
 
 % run modeling testing the relationship between severe seizures (terminated by rescue therapy) and ASM blood levels
 run ativan_sz_lme.m
 
 % run model testing if baseline seizure frequency and ASM taper affect length of stay in the epilepsy monitoring unit:
 
 % 1) using decrease in overall ASM load as predictor:
 run ASM_taper_load_LOS.m
 
 % 2) using the taper of the top 4 ASMs as predictors in model: 
 run all_ASMs_LOSmodel_load.m