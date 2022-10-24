function [protein,XT] = run_dsat(params, dl0, cact0, reltol, abstol)
% This function runs the deconvolution model with Toll saturation (dsat)
% from NC 10 -- NC 14.
% params:     a vector of parameters in linear space.
% dl0,cact0:  the initial concentrations (normalized) of dl,cact (resp)
% T:          a vector that delineates time points between mitoses and interphases.

warning off MATLAB:ode15s:IntegrationTolNotMet


if ~exist('reltol','var')
	reltol = 1e-7;
end
if ~exist('abstol','var')
	abstol = 1e-7;
end
if ~exist('dl0','var')
	dl0 = 1;
end
if ~exist('cact0','var')
	cact0 = 1;
end

nParams     = length(params);
addParams   = [];
if nParams ~= 18
    if nParams == 15
       kappa       = 1e6;                  % toll param...make large if you want unsaturable Toll
       delta       = 1e-100;               % neg fbk params 1
       K           = 0.05;                 % neg fbk params 2  
       addParams   = [kappa, delta, K];
    end
    if nParams == 16
       delta       = 1e-100;               % neg fbk params 1
       K           = 0.05;                 % neg fbk params 2  
       addParams   = [delta, K];
    end
end


%
% Load dl-Venus data
%
load('Models/dlcact/Mats/hybrid_embryo.mat','data')


%
% Initialize variables - Proteins, time & space
%
names = {'dlNuc','dlCyt','dlCactNuc','dlCactCyt','cactNuc','cactCyt'};
ts    = {'T10','T10m','T11','T11m','T12','T12m','T13','T13m','T14'};
xs    = {'X10','X10m','X11','X11m','X12','X12m','X13','X13m','X14'};


% time
% Tspan = struct( ...
%     'NC10' ,data.t(1:16),  ...
%     'NC10m',data.t(16:33), ...
%     'NC11' ,data.t(1:16)    + data.t(33), ...
%     'NC11m',data.t(17:33)   + data.t(33), ...
%     'NC12' ,data.t(34:59)   + data.t(33), ...
%     'NC12m',data.t(60:76)   + data.t(33), ...
%     'NC13' ,data.t(77:126)  + data.t(33), ...
%     'NC13m',data.t(127:148) + data.t(33), ...
%     'NC14' ,linspace(35.6,90,184)' + data.t(33));



Tspan = struct( ...
    'NC10' ,data.t(1:16),  ...
    'NC10m',data.t(16:33), ...
    'NC11' ,data.t(1:16)    + data.t(33), ...
    'NC11m',data.t(16:34)   + data.t(33), ...
    'NC12' ,data.t(34:59)   + data.t(33), ...
    'NC12m',data.t(59:77)   + data.t(33), ...
    'NC13' ,data.t(77:126)  + data.t(33), ...
    'NC13m',data.t(126:149) + data.t(33), ...
    'NC14' ,linspace(data.t(149),data.t(149+183),184)' + data.t(33));
ncs  = {'NC10','NC10m','NC11','NC11m', ...
        'NC12','NC12m','NC13','NC13m','NC14'};
nncs = length(ncs);


% space/nuclei
m = 0;
M = struct( ...
    'NC10',13, 'NC10m',13, ...
    'NC11',19, 'NC11m',19, ...
    'NC12',26, 'NC12m',26, ...
    'NC13',36, 'NC13m',36, 'NC14',51);
sub_X = struct( ...
    'NC10',[], 'NC10m',[], ...
    'NC11',[], 'NC11m',[], ...
    'NC12',[], 'NC12m',[], ...
    'NC13',[], 'NC13m',[], 'NC14',[]);


% Protein
sub_pro = struct( ...
    'NC10',zeros(13,16), 'NC10m',zeros(13,18), ...
    'NC11',zeros(19,16), 'NC11m',zeros(19,17), ...
    'NC12',zeros(26),    'NC12m',zeros(26,17), ...
    'NC13',zeros(36,50), 'NC13m',zeros(36,22), 'NC14',zeros(51,184));
protein = struct( ...
    'dlNuc',    sub_pro, 'dlCyt',    sub_pro, ...
    'dlCactNuc',sub_pro, 'dlCactCyt',sub_pro, ...
    'cactNuc',  sub_pro, 'cactCyt',  sub_pro, ...
    'nuclearDorsal',sub_pro,                  ...
	'params',params,'addParams',addParams,'T',Tspan, 'X',sub_X, 'M',M, 'dl0',dl0, 'cact0',cact0);

params      = [params, addParams];

%
% nc 14 grids
%
mgrid = M.NC14;
xgrid = linspace(0,1,mgrid)';
tlen  = structfun(@length,Tspan);
%tlen1 = [0;cumsum(tlen)];
% vtlen = [tlen1(3)+1:tlen1(4), tlen1(5)+1:tlen1(6), ...
%     tlen1(7)+1:tlen1(8), tlen1(9)+1:tlen1(10)];     % nc 11 to nc 14 without mitosis
% Tmat  = struct2cell(Tspan);
% Tmat  = cell2mat(Tmat);


% Initialize Ymat
% Ymat = zeros(6*mgrid,tlen1(end));


%
% loop over all nuclear cycles
%
fhandle = @ftns;
options = odeset('RelTol',reltol,'AbsTol',abstol,'Jacobian',@jac);
X       = cell(nncs,1); 
T       = X;
for i = 1:nncs
    nc = ncs{i};
    
    %
	% Stage-dependent factors
    %
	m_old           = m;
	m               = M.(nc);
	tspan           = Tspan.(nc);
	x               = linspace(0,1,m)';
	protein.X.(nc)  = x;

    
    %
	% Transport matrix
    %{
	e        = ones(m,1);
	P        = spdiags([e -2*e e],[-1 0 1],m,m);
	P(1,2)   = 2; 
    P(m,m-1) = 2; 
    %} 
    
    % this is optimized for parallel calculation
    ii = [1:m, 1:m-1, 2:m];     %rows:    main diagonal, right diagonal, left diagonal
    jj = [1:m, 2:m, 1:m-1];     %columns: main diagonal, right diagonal, left diagonal
    vv = [-2*ones(1,m), 2, ones(1,2*m-4), 2]; 
    P  = sparse(ii,jj,vv);
    
    %
	% Initial conditions
    %
    if i == 1                               % the ultimate initial conditions	
        un0 = zeros(m,1);
        uc0 = zeros(m,1);
        wn0 = dl0 * ones(m,1);
        wc0 = dl0 * ones(m,1);
        vn0 = cact0 * ones(m,1);
        vc0 = cact0 * ones(m,1);
        
        stage = 'interphase';
		[An,Am,~,Vn,Vc] = nuclearSize(1,'static', m, stage);
	elseif mod(i,2) == 0                    % transition from interphase to mitosis
		un0 = zeros(m,1); 
        uc0 = (Vn*un(:,end) + Vc*uc(:,end))/(Vn + Vc);
        wn0 = un0; 
        wc0 = (Vn*wn(:,end) + Vc*wc(:,end))/(Vn + Vc);
        vn0 = un0;
		vc0 = (Vn*vn(:,end) + Vc*vc(:,end))/(Vn + Vc);
		stage = 'mitosis';
		[An,Am,~,Vn,Vc] = nuclearSize(1,'static',m,stage);
	elseif mod(i,2) ~= 0                    % transition from mitosis to interphase
		x_old = linspace(0,1,m_old)';       % Interpolating onto the new x-mesh.  Also, Nuc species = Cyt  
		uc0   = interp1(x_old,uc(:,end),x); 
        un0   = uc0;
		wc0   = interp1(x_old,wc(:,end),x); 
        wn0   = wc0;
		vc0   = interp1(x_old,vc(:,end),x); 
        vn0   = vc0;
        stage = 'interphase';
		[An,Am,~,Vn,Vc] = nuclearSize(1, 'static', m, stage);
    end
    
    
    
    %
    % Solve ode
    %
	Y0   = [un0; uc0; wn0; wc0; vn0; vc0;];
	soln = ode15s(fhandle, tspan, Y0, options, params,x,P,Vn,Vc,An,Am,stage);
    if soln.stats.nsteps == 0
		error('soln.stats.nsteps was zero')
    end
    
    
    %
    % Evaluate solution
    %
	Y  = deval(soln, tspan);
	un = Y(1:m,:);
	uc = Y(1*m+1:2*m,:);
	wn = Y(2*m+1:3*m,:);
	wc = Y(3*m+1:4*m,:);
	vn = Y(4*m+1:5*m,:);
	vc = Y(5*m+1:6*m,:);
    
    
    %
    % Storing data
    %
    T{i}        = tspan;
	X{i}        = x;
	XT.(xs{i})  = x;
	XT.(ts{i})  = tspan;
%     ymat = zeros(6*mgrid,tlen(i));
    for j = 1:length(names)
        y    = Y((j-1)*m+1:j*m,:);
        %ymat((j-1)*mgrid+1:j*mgrid,:) = interp1(x,y,xgrid,'cubic');
		protein.(names{j}).(nc) = y;
    end
    %Ymat(:,tlen1(i)+1:tlen1(i+1)) = ymat;
    %protein.Ymat.(nc) = ymat;
	
    protein.soln.(nc) = soln;
    

    %
    % Store nuclear and total Dl in separate variables
    %
	protein.('nuclearDorsal').(nc)  = ...
		protein.dlNuc.(nc) + protein.dlCactNuc.(nc);
	protein.('totalDorsal').(nc)    = ...
		(Vn*(protein.dlNuc.(nc) + protein.dlCactNuc.(nc)) + ...
		Vc*(protein.dlCyt.(nc) + protein.dlCactCyt.(nc))) / (Vc + Vn);
    %protein.('nuclearDorsalMat').(nc) = ymat(1:mgrid,:) + ymat(2*mgrid+1:3*mgrid,:);
   
end

% protein.Ymat = Ymat;
% protein.Tmat = Tmat;

% nuclearDorsal = Ymat(1:mgrid,:) + Ymat(2*mgrid+1:3*mgrid,:);
% 
% protein.nuclearDlmat = nuclearDorsal(:,vtlen);
% protein.Tmat_nuclear = Tmat(vtlen);








 


