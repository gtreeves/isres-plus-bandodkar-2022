function sens_coeff = sensitivity_analysis(fcn,X,F,modelname,addParams)

nParams = size(X,2);
ntotal  = size(X,1);

sens_coeff = zeros(ntotal,nParams+1);

for i = 1:ntotal
    x  = 10.^X(i,:);
    f  = F(i);
    disp(['Set: ',num2str(i),'/',num2str(ntotal)])
    
    x2                  = repmat(x,[nParams,1]);
    d                   = x + 0.1*x;
    x2(1:nParams+1:end) = d;
    extraParams         = repmat(addParams,[nParams,1]);
    
    f1                  = feval(fcn, [log10(x2), extraParams],modelname);
    coeff               = (f1 - f)/f/0.1;
    sens_coeff(i,:)     = [coeff',f];
end

%{
labels = {'\lambda_u', '\lambda_w', ...
    '\lambda_v', '\sigma_u', '\sigma_w', '\sigma_v', '\mu_u', '\mu_w', '\mu_v',...
	'\gamma', '\psi', '\alpha', '\beta', '\phi', '\beta_0'};
if nParams == 16
    labels(end+1) = {'C_0'}; 
end
labels(end+1) = {'Error'}; 

sens_coeff2 = array2table(sens_coeff,'VariableNames',labels);
%}

end

