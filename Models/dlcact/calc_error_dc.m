function [F,penalty,K,Soln] = calc_error_dc(Params,dl0,cact0,reltol,abstol)
% Calculate the error (and penalty) of the dsat model of dlcact. Params 
% are in log10 format.

warning('off','all')

if ~exist('dl0','var')
	dl0 = 1;
end
if ~exist('cact0','var')
	cact0 = 1;
end
if ~exist('reltol','var')
	reltol = 1e-7;
end
if ~exist('abstol','var')
	abstol = 1e-7;
end




% Experimental data - dl-Venus
C    = dlVenusData('nuclear');   
C    = C/1000;
npts = length(C);


% 
% Initialize variables
%
n_sets  = size(Params,1);
F       = nan(n_sets,1);
K       = F;
Soln    = cell(n_sets,1);
penalty = nan(n_sets,2);
%opts    = optimset('Display','off');


parfor(i = 1:n_sets)
    
    % 
    % Run model
    %
    params = 10.^(Params(i,:)); 
    
    try
        solnwt    = run_dsat(params,dl0,cact0,reltol,abstol);
        Soln{i}   = solnwt;
        
        % Simulation data as a column vector
        D         = [solnwt.nuclearDorsal.NC11(:);
                     solnwt.nuclearDorsal.NC12(:);
                     solnwt.nuclearDorsal.NC13(:);
                     solnwt.nuclearDorsal.NC14(:)];
        Dlnuc     = [solnwt.dlNuc.NC11(:);
                     solnwt.dlNuc.NC12(:);
                     solnwt.dlNuc.NC13(:);
                     solnwt.dlNuc.NC14(:)];
        DlCactnuc = [solnwt.dlCactNuc.NC11(:);
                     solnwt.dlCactNuc.NC12(:);
                     solnwt.dlCactNuc.NC13(:);
                     solnwt.dlCactNuc.NC14(:)]; 
        
        %
        % Calculate Error
        %
        Beta  = mean(D.*C)/mean(D.^2);
        epsln = (C-Beta*D);
        gamma = sum(epsln.^2);
        F(i)  = gamma;
        
                
               
        % 
        % Calculating Penalties
        %
        % Calculate Penalties 
        % 1. dlNuc should be closer to gregsData('nuclear') than dlCactNuc else, penalize.
        % 2. The dorsal-most nucleus must have less free dl than dl/Cact
        % complex else, penalize
        Beta1        = mean(Dlnuc.*C)/mean(Dlnuc.^2);
        Beta2        = mean(DlCactnuc.*C)/mean(DlCactnuc.^2);
        epsln1       = (C - Beta1*Dlnuc);
        epsln2       = (C - Beta2*DlCactnuc);
        phi          = [(sum(epsln1.^2) - sum(epsln2.^2))/npts, ...
                        solnwt.dlNuc.NC14(end,end) - solnwt.dlCactNuc.NC14(end,end)];
          
 
        penalty(i,:) = phi;
                  
    catch ME
        disp(ME.message)
        F(i)            = nan;
        penalty(i,:)    = nan;
        Soln{i}         = nan;   
    end           
end


