%% function to generate the aed blood levels during a set interval of time for a single specified dose 
%% single compartment, exp decay model with rapid absorption
% something to consider later: for drugs that interact, we can multiply the
% values of one drug by some fraction when another drug is present that
% increases the metabolism of that drug 

%inputs:
% dose in mg - should be administered dose + current blood level
% tHalf - elimintaion half life in hours
% tInt - non-relative interval in hours (determined from time between medication admins of specific drug)
% F = bioavailability fraction
% vd - volume of distribution in liters
%% we want to use first order kinetics for plasma concentration over time, assuming rapid absorption 
% later, we want to include rising phase for slow absorption determined by
% time to maximal plasma concentration 


function [c_t,t] = get_single_dose_curve(c0,dose,tHalf,tInt,F,vd,tmax,ka)
%% determine parameters  
%input should be half life, then figure out what k_el is 
tInt = [0 tInt];
t = linspace(tInt(1),tInt(2),tInt(2)*60); % one point per minute for each hour in the interval
kel = -log(0.5)./tHalf;

%calculate the drug concentration over the time interval (in hours, sampled every minute between dose and next dose)
c_t = c0*exp(-kel.*t)+((ka*F*dose)./(vd*(ka-kel)))*(exp(-kel.*t)-exp(-ka.*t));

end 

% 
% % solve for ka
% syms ka
% tmax_eq = (log(ka)-log(kel))./(ka-kel)==tmax;
% kab=vpasolve(tmax_eq,ka);



%calculate the drug concentration over the time interval (in hours, sampled every minute between dose and next dose)
% **THIS IS THE SINGLE COMPARTMENT, ELIMINATION ONLY (RAPID ABSORPTION) MODEL
% c_0 = dose*F; % here we want to probably multiply by the bioavailability 
% c_t = c_0*exp(-kel.*t);
% 
% % solve for cmax, then solve for ka
% cmax = ((F*dose)/vd) * (1-exp(-kel*tmax));
% 
% syms ka %F D vd kel t
% eqn = ((ka*F*dose)/(vd*(ka-kel)))*(exp(-kel*tmax)-exp(-ka*tmax)) == cmax;
% kab=vpasolve(eqn,ka);
