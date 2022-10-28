function out = bootstrap_ci_and_p(A,B)
%{
Boos, Dennis D. "Introduction to the bootstrap world." Statistical science 18.2 (2003): 168-174.
https://www.stat.umn.edu/geyer/old/5601/examp/tests.html
Sample A and Sample B
Sample A-B = C
Can calculate CI of C
I want two-tailed p-value of is A > B (or B>A). I will calculate
1: (sum(C>=0)+1)/(length(C)+1)
2: (sum(C<=0)+1)/(length(C)+1)
Take the mininum of [1,2] and multiply this by 2 to get a 2-sided pvalue
%}

if ~exist('B','var') % if only one thing
    %% Get mean and 95% CI for each thing
    Amean = nanmean(A);
    A_CI_95 = prctile(A,[2.5,97.5]);
    p_A_bigger = (sum(A>=0)+1)/(length(A)+1);
    p_B_bigger = (sum(A<=0)+1)/(length(A)+1);
    p_value  =2*min([p_A_bigger,p_B_bigger]);
    out.mean = Amean;
    out.CI_95 = A_CI_95;
    out.p = p_value;
    
else
    
    %% Get mean and 95% CI for each thing
    Amean = mean(A);
    Bmean = mean(B);
    A_CI_95 = prctile(A,[2.5,97.5]);
    B_CI_95 = prctile(B,[2.5,97.5]);
    
    %% Get p-value for difference
    C = A-B;
    %assert(~any(isnan(C))) % throw an error if any nans
    
    p_A_bigger = (sum(C>=0)+1)/(length(C)+1);
    p_B_bigger = (sum(C<=0)+1)/(length(C)+1);
    %p_value  =2*min([p_A_bigger,p_B_bigger]); % two tailed p-value
    
    p_value = (sum(B>=A)+1)/(length(A)+1); %one-tailed p-value
    
    
    
    out.mean = [Amean Bmean];
    out.CI_95 = [A_CI_95;B_CI_95];
    out.p = p_value;
    
end


end

