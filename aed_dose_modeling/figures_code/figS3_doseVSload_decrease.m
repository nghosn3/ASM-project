load('asm_dose_decrease_tbl.mat')
load('asm_load_decrease_tbl.mat')

asm_load = tbl_load.asm_decrease;
asm_dose = tbl_decrease.asm_decrease;

% remove inf values from asm_dose
inf_inds = asm_dose == -Inf;
asm_load(inf_inds)=[];
asm_dose(inf_inds)=[];

[RHO,PVAL] = corr([asm_load asm_dose])
scatter(asm_load,asm_dose); hold on;
yvals = asm_dose;
xvals = asm_load;

coefs = polyfit(xvals,yvals,1);
aed_pred = (coefs(1)*xvals) +coefs(2);
plot(xvals,aed_pred,'k','linewidth',2)
ylabel('decrease in ASM dose');
xlabel('decrease in modeled ASM load')
legend(['R = ' num2str(RHO(2,1))],['p = ' num2str(PVAL(2,1))])
axis square;