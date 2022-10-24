function [C,dC,Cstruct,dCstruct,Cmat] = dlVenusData(options,varargin)

% dlVenusData(options) - a function that creates a column vector of the
% dl-Venus data from Reeves et al., 2012, useful for calculating error
%
% Output:
% C  : the column vector of data
% dC : the weights of each data point
% NC : a structure with nuclear cycle fieldnames and corresponding dlVenus
% data as a nuclei x time matrix
%
% Input:
% options: 'nuclear' - nuclear dl fluorescence
%          'totalDl' - total dl fluorescence
% varargin: (name-value pairs)
% IncludeMitosis: boolean variable that allows you to include mitosis data.
% Default is false.
% SingleMatrix: boolean variable that gives you a single matrix with data
% from all nuclear cycles interpolated on nuclear cycle 14 spatial grid.
% Default is false.

load('Mats/hybrid_embryo.mat')

% options
switch options
    case 'nuclear'
        A  = data.A;
        B  = data.B;
        dA = data.dA;
        dB = data.dB;
    case {'totalDl','totalDorsal','total'}
        A  = data.A2;
        B  = data.B2;
        dA = data.dA;
        dB = data.dB;
end


% additional arguments sent to function
args.IncludeMitosis = false;            %default values
args.SingleMatrix   = false;
if nargin > 1
    if mod(length(varargin),2)==0
        for i = 1:2:length(varargin)
            args.(varargin{i})=varargin{i+1};
        end
    else
        error('Varargin requires name-value pairs. Some options not specified')
    end
end

% number of nuclei structure
nuclei = struct('M',[19 26 36 51]);

%% NC 11
% N11T = t(1:16);
N11A = A(1:16);
N11B = B(1:16);
dA11 = dA(1:16);
dB11 = dB(1:16);

sig = 0.15;
m = -0.07;
M = nuclei.M(1);
x = linspace(0,1,M);
N11D = repmat(N11A,1,M).*exp(-repmat(x,length(N11A),1).^2/2/sig^2) + ...
    repmat(N11B,1,M) + m*N11A*x;
N11D = fliplr(rot90(N11D,3));

dC11 = exp(-repmat(x,length(N11A),1).^2/2/sig^2).*repmat(dA11,1,length(x)) + ...
	m*dA11*x + m*repmat(x,length(dA11),1) + repmat(dB11,1,length(x));

%% NC 11 - Mitosis
N11Am = A(16:34);
N11Bm = B(16:34);
dA11m = dA(16:34);
dB11m = dB(16:34);

sig = 0.15;
m = -0.07;
M = nuclei.M(1);
x = linspace(0,1,M);
N11Dm = repmat(N11Bm,1,M);
N11Dm = fliplr(rot90(N11Dm,3));

dC11m = exp(-repmat(x,length(N11Am),1).^2/2/sig^2).*repmat(dA11m,1,length(x)) + ...
	m*dA11m*x + m*repmat(x,length(dA11m),1) + repmat(dB11m,1,length(x));


%% NC 12
% N12T = t(34:59);
N12A = A(34:59);
N12B = B(34:59);
dA12 = dA(34:59);
dB12 = dB(34:59);
M = nuclei.M(2);
x = linspace(0,1,M);
N12D = repmat(N12A,1,M).*exp(-repmat(x,length(N12A),1).^2/2/sig^2) + ...
    repmat(N12B,1,M) + m*N12A*x;
N12D = fliplr(rot90(N12D,3));
dC12 = exp(-repmat(x,length(N12A),1).^2/2/sig^2).*repmat(dA12,1,length(x)) + ...
	m*dA12*x + m*repmat(x,length(dA12),1) + repmat(dB12,1,length(x));

%% NC 12 - Mitosis
% N12T = t(34:59);
N12Am = A(59:77);
N12Bm = B(59:77);
dA12m = dA(59:77);
dB12m = dB(59:77);
M = nuclei.M(2);
x = linspace(0,1,M);
N12Dm = repmat(N12Bm,1,M);
N12Dm = fliplr(rot90(N12Dm,3));
dC12m = exp(-repmat(x,length(N12Am),1).^2/2/sig^2).*repmat(dA12m,1,length(x)) + ...
	m*dA12m*x + m*repmat(x,length(dA12m),1) + repmat(dB12m,1,length(x));


%% NC 13
% N13T = t(77:126);
N13A = A(77:126);
N13B = B(77:126);
dA13 = dA(77:126);
dB13 = dB(77:126);
M = nuclei.M(3);
x = linspace(0,1,M);
N13D = repmat(N13A,1,M).*exp(-repmat(x,length(N13A),1).^2/2/sig^2) + ...
    repmat(N13B,1,M) + m*N13A*x;
N13D = fliplr(rot90(N13D,3));
dC13 = exp(-repmat(x,length(N13A),1).^2/2/sig^2).*repmat(dA13,1,length(x)) + ...
	m*dA13*x + m*repmat(x,length(dA13),1) + repmat(dB13,1,length(x));

%% NC 13 - Mitosis
% N13T = t(77:126);
N13Am = A(126:149);
N13Bm = B(126:149);
dA13m = dA(126:149);
dB13m = dB(126:149);
M = nuclei.M(3);
x = linspace(0,1,M);
N13Dm = repmat(N13Bm,1,M);
N13Dm = fliplr(rot90(N13Dm,3));
dC13m = exp(-repmat(x,length(N13Am),1).^2/2/sig^2).*repmat(dA13m,1,length(x)) + ...
	m*dA13m*x + m*repmat(x,length(dA13m),1) + repmat(dB13m,1,length(x));

%% NC 14
% N14T = t(149:332);
N14A = A(149:332);
N14B = B(149:332);
dA14 = dA(149:332);
dB14 = dB(149:332);

M = nuclei.M(4);
x = linspace(0,1,M);
N14D = repmat(N14A,1,M).*exp(-repmat(x,length(N14A),1).^2/2/sig^2) + ...
    repmat(N14B,1,M) + m*N14A*x;
N14D = fliplr(rot90(N14D,3));
dC14 = exp(-repmat(x,length(N14A),1).^2/2/sig^2).*repmat(dA14,1,length(x)) + ...
	m*dA14*x + m*repmat(x,length(dA14),1) + repmat(dB14,1,length(x));


%% Fluorescence data
% check whether to include mitosis data
if args.IncludeMitosis
    C        = [N11D(:);N11Dm(:);N12D(:);N12Dm(:);N13D(:);N13Dm(:);N14D(:)];
    dC       = [dC11(:);dC11m(:);dC12(:);dC12m(:);dC13(:);dC13m(:);dC14(:)];
    Cstruct  = struct('NC11',N11D,'NC11m',N11Dm,'NC12',N12D,'NC12m',N12Dm, ...
                      'NC13',N13D,'NC13m',N13Dm,'NC14',N14D);
    dCstruct = struct('NC11',dC11','NC11m',dC11m','NC12',dC12','NC12m',dC12m', ...
                      'NC13',dC13','NC13m',dC13m','NC14',dC14');
else
    C        = [N11D(:);N12D(:);N13D(:);N14D(:)];
    dC       = [dC11(:);dC12(:);dC13(:);dC14(:)];
    Cstruct  = struct('NC11',N11D,'NC12',N12D,'NC13',N13D,'NC14',N14D);
    dCstruct = struct('NC11',dC11','NC12',dC12','NC13',dC13','NC14',dC14');
end


% check if single matrix is requested
if args.SingleMatrix
   fields   = fieldnames(Cstruct);
   xspan    = linspace(0,1,size(Cstruct.(fields{end}),1));
   Cmat     = [];   
   for i = 1:length(fields)
      nc        = fields{i}; 
      ns        = size(Cstruct.(nc),1); 
      Cmat.(nc) = interp1(linspace(0,1,ns),Cstruct.(nc),xspan,'spline');
   end
end




