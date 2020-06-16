%% Wrap function for integration of 7 box model
%
% This function calls the ODE-solver and returns:
% - t: a vector containing time points (days)
% - y: a matrix containing the 14 state variables (first the DOC in the 
%       seven boxes, then the bacterial biomass, according to indices 
%       "Jdom" and "Jbac" from setup_model.m) for each time point
% - PE, PO, PD: structures containing parameters defined in setup_model.m,
%       setup_ocean.m, and setup_run.m, respectively

function [t, y, PE, PO, PD] = wrap_boxmodel(varargin)

% Print model progress (or suppress by 'silent')
silent = 0;
if any(strcmp(varargin, 'silent'))
    silent = 1;
end
if silent == 0
    fprintf('\nPreparing...')
end

% Default: no one-way flow (circulation instead)
oneway = 0;
if any(strcmp(varargin, 'oneway'))
    oneway = 1;
end

% Optional input argument: Number of boxes
nb  = 7; % Default value
if any(strcmp(varargin, 'nb'))
    ind  = find(strcmp(varargin, 'nb'));
    nb = varargin{ind+1};
end

% Optional input argument: Circulation strength factor
Psi_factor  = 1; % Default value
if any(strcmp(varargin, 'Psi_factor'))
    ind  = find(strcmp(varargin, 'Psi_factor'));
    Psi_factor = varargin{ind+1};
end

% Optional input argument: Length of spin-up period
yspin  = 990; % Default value
if any(strcmp(varargin, 'yspin'))
    ind  = find(strcmp(varargin, 'yspin'));
    yspin = varargin{ind+1};
end
PD.yspin = yspin;

% Optional input argument: Length of run period
yrun  = 10; % Default value
if any(strcmp(varargin, 'yrun'))
    ind  = find(strcmp(varargin, 'yrun'));
    yrun = varargin{ind+1};
end

% Optional input argument: total simulation length
if any(strcmp(varargin, 'yend'))
    ind  = find(strcmp(varargin, 'yend'));
    yend = varargin{ind+1};
    yrun  = 0.1*yend;
    yspin = 0.9*yend;
end
yend = yspin+yrun;

% Model parameters 
PE = setup_model('nb', nb);
PO = setup_ocean(PE, 'Psi_factor', Psi_factor, 'oneway', oneway);

% Box colors
if nb == 7
    PD.cols = parula(7);
    PD.size_scatter = [50 38 24 14 7 3 1];
    PD.BoxAbbr = {'AA', 'SA', 'LL', 'NA', 'TC', 'NADW', 'AABW'};
elseif nb == 2
    PD.cols = [mycolors('red'); mycolors('blue')];
    PD.size_scatter = [38 3];
    PD.BoxAbbr = {'Surface', 'Deep'};
else
    PD.cols = mycolors('green');
    PD.size_scatter = 10;
    PD.BoxAbbr = {'Box'};
end

%% DOM model parameters

%%%%%%%%%%
% Optional: variation of paramters by specified factors
%%%%%%%%%%

numM   = repmat(100 , nb, 1); % number of different DOM compounds
if any(strcmp(varargin, 'numM_var'))
    ind  = find(strcmp(varargin, 'numM_var'));
    numM = numM.*varargin{ind+1}';
end
if any(strcmp(varargin, 'numD_var'))
    ind  = find(strcmp(varargin, 'numD_var'));
    numM = numM.*varargin{ind+1}';
end

nsubs  = repmat(3 , nb, 1);  % number of substrates each bacteria can take up
if any(strcmp(varargin, 'nsubs_var'))
    ind  = find(strcmp(varargin, 'nsubs_var'));
    nsubs = nsubs.*varargin{ind+1}';
end

% Bacterial parameters, uptake according to Michaelis Menten kinetics
r_max   = ones(nb, 1); % maximum uptake rate (vmax, 1/d)
if any(strcmp(varargin, 'rmax_var'))
    ind  = find(strcmp(varargin, 'rmax_var'));
    r_max = r_max.*varargin{ind+1}';
end
if any(strcmp(varargin, 'r_max_var'))
    ind  = find(strcmp(varargin, 'r_max_var'));
    r_max = r_max.*varargin{ind+1}';
end

K   = repmat(10 , nb, 1); % Half-saturation constant [mmol/m³]
if any(strcmp(varargin, 'K_var'))
    ind  = find(strcmp(varargin, 'K_var'));
    K = K.*varargin{ind+1}';
end

r_mort1 = repmat(0.02 , nb, 1); % linear bacterial mortality rate [1/d]
if any(strcmp(varargin, 'r_mort1_var'))
    ind  = find(strcmp(varargin, 'r_mort1_var'));
    r_mort1 = r_mort1.*varargin{ind+1}';
end
if any(strcmp(varargin, 'rmort_var'))
    ind  = find(strcmp(varargin, 'rmort_var'));
    r_mort1 = r_mort1.*varargin{ind+1}';
end

% Assignment fractions
eta     = repmat(0.2 , nb, 1); % bacterial growth efficiency, fraction of uptake converted to biomass
if any(strcmp(varargin, 'eta_var'))
    ind  = find(strcmp(varargin, 'eta_var'));
    eta = eta.*varargin{ind+1}';
end

beta    = 0.14./(1-eta);  % fraction of non-biomass uptake that is excreted 
%(i.e. beta*(1-eta) = fraction of uptake excreted as transformed DOM)
if any(strcmp(varargin, 'beta_var'))
    ind  = find(strcmp(varargin, 'beta_var'));
    beta = beta.*varargin{ind+1}';
end

% If Q10-factor is provided, implement temperature-dependency of uptake
if any(strcmp(varargin, 'Q10'))
    ind  = find(strcmp(varargin, 'Q10'));
    
    % Temperature dependence after Blackford et al 2004 (fast decrease at
    % low temperatures, slower at high temps due to enzyme inhibition)
    Q10 = varargin{ind+1}';
    fT = Q10.^((PO.T'-10)/10) - Q10.^((PO.T'-32)/3);
    r_max = r_max.*fT;
    PD.fT = fT;
    
end
if any(strcmp(varargin, 'Q10m'))
    ind  = find(strcmp(varargin, 'Q10m'));
    
    % Temperature dependence
    Q10m = varargin{ind+1}'; 
    fTm = Q10m.^((PO.T'-10)/10);
    assert(all(fTm>0) && all(~isinf(fTm)), 'fTm is negative or inf')
    warning('Enzyme inhibition left out for mortality')
    r_mort1 = r_mort1.*fTm;
    
%     % Temperature dependence after Blackford et al 2004 (fast decrease at
%     % low temperatures, slower at high temps due to enzyme inhibition)
%     Q10m = varargin{ind+1}';
%     fTm = Q10m.^((PO.T'-10)/10) - Q10m.^((PO.T'-32)/3);
%     r_mort1 = r_mort1.*fTm;
    
end

% Uptake factor (mashing all relevant parameters together)
uf   = (nsubs./numM).*(r_max./K);

%%%%%%%%%%
% SUPPLY
%%%%%%%%%%

% Production vector of DOC [mmol/m³/d]
Svec = 1e03*(24*60*60)*[2.95e-10,5.57e-10,1.58e-09,1.21e-09,1.242e-10,7.035e-11,2.9e-11]';
PD.Svec_def = Svec;

% Seasonal production vector of DOC [mmol/m³/d]
Svec_monthly = [0.0420955839269280,0.0758099959394937,0.123038355378408,0.0253153831889454,0.0107308800000000,0.00607824000000000,0.00250560000000000;0.0318971510903888,0.0675945769716677,0.126112059155261,0.0305228499296810,0.0107308800000000,0.00607824000000000,0.00250560000000000;0.0182012505161531,0.0543277710463580,0.128139985950087,0.0417541124314494,0.0107308800000000,0.00607824000000000,0.00250560000000000;0.0122305674578453,0.0392470600792018,0.132615521174355,0.0802144716735983,0.0107308800000000,0.00607824000000000,0.00250560000000000;0.0121355324364422,0.0276786726141654,0.135913375992018,0.156142981811362,0.0107308800000000,0.00607824000000000,0.00250560000000000;0.0121355324364422,0.0311231559191571,0.139209509325054,0.198306959266959,0.0107308800000000,0.00607824000000000,0.00250560000000000;0.0121355324364422,0.0237646444811208,0.144042318684890,0.208190616872910,0.0107308800000000,0.00607824000000000,0.00250560000000000;0.0120404974150393,0.0256615667688887,0.148100486371975,0.169077872605272,0.0107308800000000,0.00607824000000000,0.00250560000000000;0.0151317848129810,0.0358764809377451,0.146660085412137,0.121734967070900,0.0107308800000000,0.00607824000000000,0.00250560000000000;0.0226210118960864,0.0516748776573133,0.141598520132259,0.0789287503715146,0.0107308800000000,0.00607824000000000,0.00250560000000000;0.0323071448362051,0.0682285035733169,0.135650273284659,0.0442710474934833,0.0107308800000000,0.00607824000000000,0.00250560000000000;0.0431716627435129,0.0766627578091411,0.133291305711099,0.0347932153412143,0.0107308800000000,0.00607824000000000,0.00250560000000000];

% Optional input: supply (else using default from PO)
if any(strcmp(varargin, 'Svec'))
    ind  = find(strcmp(varargin, 'Svec'));
    Svec = varargin{ind+1};
elseif any(strcmp(varargin, 'Svec')) && any(strcmp(varargin, 'seasonal'))
    error('seasonal supply cannot be combined with specified supply!')
end

% Modify supply if 2- or 1-box model is chosen
if nb == 2
    if length(Svec) == 7
        Svec = [sum(Svec(1:4)'.*(PO.volboxweight(1:4)/sum(PO.volboxweight(1:4)))); ...
            sum(Svec(5:7)'.*(PO.volboxweight(5:7)/sum(PO.volboxweight(5:7))))];
    else
        error('Wrong dimension of supply vector')
    end
elseif nb == 1
    Svec = mean(Svec);
end

% Calculate equilibrium concentrations (equation 3 and 4)
Bstar   = eta.*Svec./(r_mort1.*(1-beta).*(1-eta));
Dstar   = r_mort1.*K.*numM./(eta.*nsubs.*r_max);
x = sprintf('%1.2f\t', Dstar);
y = sprintf('%1.2f\t', Bstar);
if silent == 0
    fprintf('\nD* [mmolC/m³]: %s', x)
    fprintf('\nB* [mmolC/m³]: %s', y)
end

% save all variables into structure
PD.nsubs   = nsubs;
PD.r_max   = r_max;
PD.K       = K;
PD.r_mort1 = r_mort1;
PD.eta     = eta;
PD.beta    = beta;
PD.uf      = uf;
PD.Svec    = Svec;
PD.Dstar   = Dstar;
PD.Bstar   = Bstar;
PD.numM    = numM;

% if supply should be seasonal, take monthly supply rates
if any(strcmp(varargin, 'seasonal'))
    if nb==7
        PD.Svec = Svec_monthly'; 
    elseif nb==2
        PD.Svec = [mean(Svec_monthly(:,1:4),2)'; mean(Svec_monthly(:,5:7),2)']; 
    elseif nb==1
        PD.Svec = mean(Svec_monthly,2)'; 
    end
else
    PD.Svec = repmat(Svec, 1, 12);
end

% construct x and v for interpolation of monthly Svec values at every
% second
x = [0 (1:(365/12):365) 365]; % first day of every month
v = NaN(nb,12+2);
for b=1:nb % go through the boxes
    v(b,:) = [mean(PD.Svec(b,[1 12])) PD.Svec(b,:) mean(PD.Svec(b,[1 12]))]; % values from Svec corresponding to x
end
PD.x = x;
PD.v = v;

%% Initial state

% DOC reference data from CLIVAR cruises (via "load_clivar_DOC.m")
DOC_clivar = [44.7, 47.4, 65.7, 54.2, 47.2, 42.7, 40.7]; % [µmolC/kg] / [mmolC/m³]
PD.DOC_reference  = DOC_clivar; % save to structure

% Biomass reference data
biomass_data = [0.3776, 0.6454, 0.4897, 0.4932, 0.1626, 0.0921, 0.0380]; % [mmolC/m³]
PD.biomass_reference = biomass_data;

% Initial state
if nb == 7
    DOM_o = DOC_clivar; %repmat(10, 1, 7) *1e-3; % initial ocean DOM (mmol/m3)
    Bac_o = biomass_data;%repmat(0.1, 1, 7)*1e-3;% initial ocean bacterial biomass carbon (mmol/m3)
elseif nb == 2
    DOM_o = [mean(DOC_clivar(1:4)*1e-3) mean(DOC_clivar(5:7)*1e-3)]; 
    Bac_o = [mean(biomass_data(1:4)*1e-3) mean(biomass_data(5:7)*1e-3)];
elseif nb == 1
    DOM_o = mean(DOC_clivar*1e-3); 
    Bac_o = mean(biomass_data*1e-3);
end

% Optional input argument: Initial B and D
if any(strcmp(varargin, 'DOM_o'))
    ind  = find(strcmp(varargin, 'DOM_o'));
    DOM_o = varargin{ind+1};
end
if any(strcmp(varargin, 'Bac_o'))
    ind  = find(strcmp(varargin, 'Bac_o'));
    Bac_o = varargin{ind+1};
end

PD.DOM_o = DOM_o;
PD.Bac_o = Bac_o;

% package initial values into matrix  
y0mat          = PE.m0; 
y0mat(PE.Idom) = DOM_o;
y0mat(PE.Ibac) = Bac_o;

%% Integrate differential equations

if silent == 0
    fprintf('\nIntegrating model forward...')
end
y0    = y0mat(PE.Ires); 
trun  = 0 : 1 : yend*365; %[d] 
options = odeset('abstol', 1e-18, 'reltol', 1e-13,...
    'NonNegative', ones(size(y0)));
[t,y] = ode15s(@carbon_climate_derivs, trun, y0', options, PE, PO, PD); 


if silent == 0
    fprintf('\nDone.\n')
end


end