function [dydt] = carbon_climate_derivs(t,y, PE, PO, PD)
%% Create local copies of state variables

Mloc = y(PE.Jdom)'; % DOC [mmol/m³]
Bloc = y(PE.Jbac)'; % Bacterial carbon biomass [mmol/m³]


%% update basic quantities

Psi   = PO.Psi_o; % [1/d]
Svec_m = NaN(PE.nb,1);

% Interpolate for current second between monthly values (smooth)
for b = 1:PE.nb
    d_in_year = mod(t, 365);
    Svec_m(b) = interp1(PD.x, PD.v(b,:), d_in_year);
end


%% DOM ecosystem model [basic unit: mmolC/m³/d]

Beco = PD.eta.*PD.uf.*Bloc'.*Mloc' - PD.r_mort1.*Bloc';
Meco = PD.beta.*(1-PD.eta).*PD.uf.*Mloc'.*Bloc' + PD.r_mort1.*Bloc' - PD.uf.*Mloc'.*Bloc';


%% Compute Tendencies

dMdt = (Psi*Mloc' + Meco) + Svec_m; % DOC [mmol/m³]
dBdt = Psi*Bloc' + Beco; % Bacteria [mmol/m³]

%% matrix of derivatives

dydtmat = PE.m0;
dydtmat(PE.Idom) = dMdt; % Ocean Dissolved Organic Carbon [mmolC/m³]
dydtmat(PE.Ibac) = dBdt; % Ocean Bacterial Biomass Carbon [mmolC/m³]
dydt = dydtmat(PE.Ires)';

end
