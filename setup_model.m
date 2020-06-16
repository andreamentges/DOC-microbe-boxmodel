function Var = setup_model(varargin)

nd = 2; % number state variable 'domains' (Ocean DOM, Ocean Bacteria)

% Optional input argument: Number of boxes
nb  = 7; % Default value
if any(strcmp(varargin, 'nb'))
    ind  = find(strcmp(varargin, 'nb'));
    nb = varargin{ind+1};
end

%% Earth parameters 

re       = 6371e3; % Earth radius [m]
Ae       = 4*pi*re^2; % Earth area [m²]
fla      = 0.3; % fraction land area
foc      = 1-fla; % fraction ocean area
ps       = 1013.5 * (1e2); % mean surface pressure [Pa]

%% Structure and Indices of state variable arrays

% generic matrices
m0 = zeros(nb,nd)*nan; 
m1 = m0+1;

% indices of pools in generic matrices
Idom = sub2ind(size(m0),1:nb,zeros(1,nb)+1);
Ibac = sub2ind(size(m0),1:nb,zeros(1,nb)+2);
Ires = cat(2,Idom,Ibac);

% indices of pools in vector of m0(Ires)
Jdom = 0*nb+1:1*nb; 
Jbac = 1*nb+1:2*nb;

%% save all variable into structure

Var=v2struct;

end




