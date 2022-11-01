function [out,pval_binom_alt,successes,successes_val] = plot_orders(rates,sozs)

myColours = [0, 0.4470, 0.7410;...
    0.8500, 0.3250, 0.0980;...
    0.4660, 0.6740, 0.1880;...
    0.4940, 0.1840, 0.5560;...
    0.6350, 0.0780, 0.1840];
grayColor = 0.75*[1 1 1];
markersize = 3;
npts = length(rates);


%[all_ranks,all_soz_ranks,nchance,all,successes] = get_ranks(things,sozs,rates,which,min_rate);
[all_ranks,all_soz_ranks,successes,successes_val,soz_chance,aed_load_data] = simple_rate_rank(rates,sozs);
% Alternate statistical test
%{
In this test, I see whether the median rate is higher for SOZ than chance
for each patient (whereas in the prior I see if the median rate RANKING is
higher for SOZ than chance). I think this is better because when I take
ranks matlab will force equal things to have different ranks, and so this
seems more accurate. The p values are slightly different but both <0.001
%}
if 0
    soz_chance(any(isnan(soz_chance),2),:) = []; % remove nans
    npts_alt = size(soz_chance,1);
    alt_successes = soz_chance(:,1) > soz_chance(:,2);
    pval_binom_alt = 2*binocdf(sum(alt_successes==0),npts_alt,0.5);
    
end

%% Re-order by time bin ...# of electrodes
num_bins = cellfun(@length,all_ranks); % the ranks should be the time bin

% remove those with no elecs
% no_elecs = num_bins == 0;
% all_ranks(no_elecs) = [];
% all_soz_ranks(no_elecs) = [];
% successes(no_elecs) = [];

npts = length(num_bins);

[num_bins,I] = sort(num_bins);
all_ranks = all_ranks(I);
all_soz_ranks = all_soz_ranks(I);


%% Plot stuff
for i = 1:npts
    plot(i,all_ranks{i},'o','color',grayColor,'markersize',markersize);
    hold on
    sp = plot(i,nanmedian(all_soz_ranks{i}),'p',...
        'markerfacecolor',[0.2660 0.6040 0.2880],'markeredgecolor',[0.2660 0.6040 0.2880],...
        'linewidth',2,'markersize',markersize+5);
    cp = plot(i,median(1:num_bins(i)),'o','markeredgecolor',[0.3 0.3 0.3],...
        'markerfacecolor',[0.25 0.25 0.25],...
        'linewidth',2,'markersize',markersize+2);
    
end


pval_binom_alt = 2*binocdf(sum(successes_val==0),length(successes),0.5);
p_signed = signrank(aed_load_data(:,1),aed_load_data(:,2));

% Get median across all patients
all_all_soz_ranks = [];
for i = 1:npts
    all_all_soz_ranks = [all_all_soz_ranks;(all_soz_ranks{i})'];
    
end
median_rank = median(all_all_soz_ranks);

out.n = length(successes);
out.nsuc = sum(successes==1);
out.pval = pval_binom_alt;
out.p_signed = p_signed;
out.median_rank = median_rank;

xlim([0 npts])
xl = xlim;
xbar = xl(1) + 1.02*(xl(2)-xl(1));
xtext = xl(1) + 1.04*(xl(2)-xl(1));
newxl = [xl(1) xl(1) + 1.05*(xl(2)-xl(1))];
plot([xbar xbar],[1 num_bins(end)],'k-','linewidth',2)
if pval_binom_alt >= 0.05
text(xtext-1,(1+num_bins(end))/2,get_asterisks(pval_binom_alt,1),'rotation',90,...
        'horizontalalignment','center','fontsize',16)

else
    text(xtext,(1+num_bins(end))/2,get_asterisks(pval_binom_alt,1),'rotation',90,...
        'horizontalalignment','center','fontsize',20)
end
xlim(newxl)
legend([sp,cp],{'Median Sz AED load', 'Chance'},'fontsize',15,'location','northwest')
ylabel('time bin (ranked)','fontsize',14); xlabel('Patient','fontsize',14)

end