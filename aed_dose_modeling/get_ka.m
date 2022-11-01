
function [ka] = get_ka(tHalf,tmax)
%% determine parameters  
%input should be half life, then figure out what k_el is 
kel = -log(0.5)./tHalf;

% solve for ka
syms ka
tmax_eq = (log(ka)-log(kel))./(ka-kel)==tmax;
ka=vpasolve(tmax_eq,ka);


%calculate the drug concentration over the time interval (in hours, sampled every minute between dose and next dose)

end 

