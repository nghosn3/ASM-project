% for getting the coefficient estimates and the correspondinf g p-values
function [coeff_names,coeff_stats,all_coefs] = model_bootstrap_stats(drug_table,feat_names,model,nb,ptIDs)

% get bootstrapped samples - sample with replacement
pt_inds_master = cell(1,nb);
for i = 1:nb
    pt_inds_master{i} = datasample(1:length(ptIDs),length(ptIDs)); %sameple with replacement
end

ncoeffs = length(feat_names)-2; %exclude ptID and seizure
pt_stats = nan(nb,ncoeffs);
coeff_names_expected = feat_names(1:end-2);

for ib = 1:nb
    train_pts = pt_inds_master{ib};
    
    % get features for the specific sample
    
    feats_train = [];
    for j = 1:length(train_pts)
        feats_train = [feats_train; drug_table(contains(cellstr(drug_table.ptID),num2str(train_pts(j))),:)];
    end
    
    % run null model or ASM model?
    if model
        try
            [mout] = run_logistic_on_AED(feats_train,feat_names);
            coeff_names = mout.CoefficientNames(2:end);
            assert(isequal(coeff_names',coeff_names_expected))
            %if ~isfield(mout,'labels'), continue; end
            for ic = 1:ncoeffs
                pt_stats(ib,ic) = (mout.Coefficients{ic+1,2}); % *RE-RUN try estimate - odds ratio - exp(estimate)
            end
            
            coeff_stats = nan(ncoeffs,4); % mean, lower CI, higher CI, p (meaningless here)
            %OP_stats = nan(ncoeffs,3); % mean, lower CI, higher CI, p (meaningless here)
            for ic = 1:ncoeffs
                tout = bootstrap_ci_and_p(squeeze(pt_stats(:,ic)));
                coeff_stats(ic,:) = [tout.mean,tout.CI_95,tout.p];
                
                
            end
         catch
             disp(['asm model: error on iter ' num2str(ib)]);
         end
    else
         try
            [mout] = run_NULL_logistic_on_AED(feats_train,feat_names);
            
            coeff_names = mout.CoefficientNames(2:end);
            assert(isequal(coeff_names',coeff_names_expected))
            %if ~isfield(mout,'labels'), continue; end
            for ic = 1:ncoeffs
                pt_stats(ib,ic) = (mout.Coefficients{ic+1,2}); % *RE-RUN try estimate - odds ratio - exp(estimate)
            end
            
            
            
         catch
             disp(['null model: error on iter ' num2str(ib)]);
         end
    end
      
end

coeff_stats = nan(ncoeffs,4); % mean, lower CI, higher CI, p (meaningless here)
%OP_stats = nan(ncoeffs,3); % mean, lower CI, higher CI, p (meaningless here)
for ic = 1:ncoeffs
    tout = bootstrap_ci_and_p(squeeze(pt_stats(:,ic)));
    coeff_stats(ic,:) = [tout.mean,tout.CI_95,tout.p]; 
end

all_coefs.names = coeff_names;
all_coefs.coefs = pt_stats;
