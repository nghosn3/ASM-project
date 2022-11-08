function [OPs,AUCs,Xroc,Yroc,conf_matrices] = model_bootstrap_ROC(drug_table,feat_names,model,nb,pt_inds_master)

pt_op = nan(nb,1);
pt_AUCs = nan(nb,1);
conf_matrices = cell(nb,1);
all_xroc=zeros(1000,nb);
all_yroc=zeros(1000,nb);
for ib = 1:nb
    
    curr_pts = pt_inds_master{ib};
    split = round(length(curr_pts)*.8);
    
    train_pts = curr_pts(1:split);
    test_pts = curr_pts(split+1:end);
    
    feats_train = [];
    for j = 1:length(train_pts)
        feats_train = [feats_train; drug_table(contains(cellstr(drug_table.ptID),num2str(train_pts(j))),:)];
    end
    
    % run null model or ASM model?
    if model
        [mout] = run_logistic_on_AED(feats_train,feat_names); %drug table should be the training data, test outside
    else
        
        [mout] = run_NULL_logistic_on_AED(feats_train,feat_names);
    end
    
    if ~isempty(mout)
    % get testing set features - same testing set for all iterations
    feats_test = [];
    for j = 1:length(test_pts)
        feats_test = [feats_test; drug_table(contains(cellstr(drug_table.ptID),num2str(test_pts(j))),:)];
    end
    labels = feats_test.seizure; %seizure labels
    
    % get model prediction on test set
    sz_pred = predict(mout,feats_test);
    
    % from sz_pred and labels get predictions for each class
    class_sz = sz_pred(logical(labels));
    class_no_sz = sz_pred(~labels);
    
%     disc = 0.01; % operating point at 0.01 - determined by looking atthe distrubution of prediction probabilities
%     conf_matrix = make_conf_matrix(class_sz,class_no_sz,disc);
%     
    % get ROC curve to get the AUC and operating point
    intervals= linspace(0, 1, 1000);
    [Xrocs, Yrocs, ~, pt_AUCs(ib)]= perfcurve(labels,sz_pred,1,'Xvals',intervals);
    try
    all_xroc(1:length(Xrocs),ib)=Xrocs;
    all_xroc(1:length(Yrocs),ib)=Yrocs;
    catch
        disp('error')
    end 
    
    % Get operating point at 0.1 FPR
    [~,ind]=min(abs(Xrocs-.1));
    pt_op(ib) = Yrocs(ind);
    
    end 
%  conf_matrices{ib} = conf_matrix;
  
end

AUCs = (pt_AUCs);
OPs = (pt_op);
Xroc = all_xroc;
Yroc = all_yroc;



end