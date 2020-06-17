%% Main file to prepare figures of DOC-microbe boxmodel
%
% This is the computer code associated to the following manuscript:
%
% Title:﻿Microbial physiology governs the oceanic distribution of dissolved 
% organic carbon in a scenario of equal degradability
%
% Authors: A. Mentges, C. Deutsch, C. Feenders, S. T. Lennartz, B. Blasius,
% T. Dittmar
%
% Correspondence: andrea.mentges@uol.de

clear variables
close all

%% Default boxmodel (section 3.1)

% Run default model simulation (default simulation time is 1000 years, 
% including a 990 year spin-up phase and a 10 year model run)
[t, y, PE, PO, PD] = wrap_boxmodel();

% Plot resulting concentrations over time
plot_boxmodel(t, y, PE, PO, PD, 'time')

% Plot resulting concentrations in the boxes (spatially)
plot_boxmodel(t, y, PE, PO, PD, 'boxes')

% Save DOC and biomass concentration at the end of the simulation
D.def = y(end, PE.Jdom);
B.def = y(end, PE.Jbac);

% Query D* and B*
D.star = PD.Dstar';
B.star = PD.Bstar';

% Query reference concentrations for DOC and biomass (see section 2.4)
D.ref = PD.DOC_reference;
B.ref = PD.biomass_reference;

%% Modification: Seasonal supply (section 3.2)

% Run model with seasonally varying DOC supply
[t, y, PE, PO, PD] = wrap_boxmodel('seasonal');

% Plot concentrations over time
plot_boxmodel(t, y, PE, PO, PD, 'time')

% Save time-averaged mean concentrations of the model run (last 90 years)
D.s_seasonal = mean(y(t>90, PE.Jdom));
B.s_seasonal = mean(y(t>90, PE.Jbac));

% Plot over time (only the last 5 years)
plot_boxmodel(t, y, PE, PO, PD, 'shorttime', 5)

% Plot resulting concentrations in the boxes (spatially)
plot_boxmodel(t, y, PE, PO, PD, 'boxes')
p = findobj(gcf,'Type','patch');
p(2).CData = D.s_seasonal;
p(1).CData = B.s_seasonal;

%% Modification: Temperature-dependence (section 3.3)

% Run model with temperature-dependent uptake
uptake_Q10 = 2;
[t, y, PE, PO, PD] = wrap_boxmodel('Q10', uptake_Q10);

% save final concentrations
D.temp = y(end, PE.Jdom);
B.temp = y(end, PE.Jbac);

% Plot resulting concentrations over time
plot_boxmodel(t, y, PE, PO, PD, 'time')

% Plot resulting concentrations in the boxes (spatially)
plot_boxmodel(t, y, PE, PO, PD, 'boxes')


%% Modification: Depth-dependent physiology (section 3.4)

% Increased mortality
[t, y, PE, PO, PD] = wrap_boxmodel('rmort_var', [2 2 2 2 1 1 1], 'yspin', 1000); 
fprintf('%1.2f\t',y(end, PE.Jdom))
D.variant_rmort = y(end, PE.Jdom);
B.variant_rmort = y(end, PE.Jbac);
plot_boxmodel(t, y, PE, PO, PD, 'time')
plot_boxmodel(t, y, PE, PO, PD, 'boxes')

% Reduced growth efficiency
[t, y, PE, PO, PD] = wrap_boxmodel('eta_var', [1/2 1/2 1/2 1/2 1 1 1], 'yspin', 1000); 
fprintf('%1.2f\t',y(end, PE.Jdom))
D.variant_eta = y(end, PE.Jdom);
B.variant_eta = y(end, PE.Jbac);
plot_boxmodel(t, y, PE, PO, PD, 'time')
plot_boxmodel(t, y, PE, PO, PD, 'boxes')

% Reduced number of substrates
[t, y, PE, PO, PD] = wrap_boxmodel('nsubs_var', [1/2 1/2 1/2 1/2 1 1 1], 'yspin', 1000);
D.variant_nsubs = y(end, PE.Jdom);
B.variant_nsubs = y(end, PE.Jbac);
plot_boxmodel(t, y, PE, PO, PD, 'time')
plot_boxmodel(t, y, PE, PO, PD, 'boxes')

% Increased DOC diversity
[t, y, PE, PO, PD] = wrap_boxmodel('numD_var', [2 2 2 2 1 1 1], 'yspin', 1000); 
D.variant_n = y(end, PE.Jdom);
B.variant_n = y(end, PE.Jbac);
plot_boxmodel(t, y, PE, PO, PD, 'time')
plot_boxmodel(t, y, PE, PO, PD, 'boxes')


%% Table 3: Overview table of concentrations

[D.ref' D.star' D.def' D.s_seasonal' D.temp' D.variant_rmort']
[B.ref' B.star' B.def' B.s_seasonal' B.temp' B.variant_rmort']

%% Figure 3: Simulations comparison plot

DOM_conc = [D.ref; D.def; D.temp; D.variant_rmort];
Bac_conc = [B.ref; B.def; B.temp; B.variant_rmort];
[cx, cy] = plot_boxmodel(t, y,PE, PO, PD, 'noplot');
n = size(DOM_conc,1);

figure('color', 'white', 'position', [545,194,686,748])
for i = 1:n

    subplot(n,2,i*2-1, 'YDir', 'reverse')
    patch(cx, cy, DOM_conc(i,:))
    ylabel('Depth [m]'), axis tight
    c = colorbar; ylabel(c, sprintf('DOC\n[mmolC/m³]'))
    c.Position([1 3 4]) = [0.38 0.02 0.158];
    colormap(parula(100))
    caxis([27 65])
    c.Ticks = 30:10:70;
    ax = gca;
    ax.XTick = [];
    ax.YTick = [0, 1000, 2000, 3000];
    if i == n
        ax.XTick = ax.XLim;
        ax.XTickLabel = {'South', 'North'};
    end
    ca = gca;
    ca.Position(3:4) = [0.23 0.1577];
    mycmap = parula(100); colormap(gca, mycmap(51:end,:))

    subplot(n,2,i*2, 'YDir', 'reverse')
    patch(cx, cy, Bac_conc(i,:))
    axis tight
    c = colorbar; ylabel(c, sprintf('Bacterial biomass\n[mmolC/m³]'))
    c.Ticks = (0:0.5:2);
    c.Position([1 3 4]) = [0.826 0.02 0.158];
    mycmap = parula(100); colormap(gca,mycmap(10:end-10,:))
    caxis([0 2])
    ax = gca;
    ax.XTick = [];
    ax.YTick = [];
    if i == n
        ax.XTick = ax.XLim;
        ax.XTickLabel = {'South', 'North'};
    end
    ca = gca;
    ca.Position(3:4) = [0.23 0.1577];

end

x = 0.005;
ax = axes('Units', 'Normal', 'Position', [.075 .075 .85 .86+x], 'Visible','off');
set(get(ax,'Title'),'Visible','on')
title('Reference DOC and biomass concentration',...
    'fontweight', 'bold', 'FontAngle', 'italic');
ax = axes('Units', 'Normal', 'Position', [.075 .075 .85 .64+x], 'Visible','off');
set(get(ax,'Title'),'Visible','on')
title('Default run',...
    'fontweight', 'bold', 'FontAngle', 'italic');
ax = axes('Units', 'Normal', 'Position', [.075 .075 .85 .42+x], 'Visible','off');
set(get(ax,'Title'),'Visible','on')
title('Modification: Temperature-dependent DOC uptake',...
    'fontweight', 'bold', 'FontAngle', 'italic');
ax = axes('Units', 'Normal', 'Position', [.075 .075 .85 .20+x], 'Visible','off');
set(get(ax,'Title'),'Visible','on')
title('Modification: Sea surface parameterization',...
    'fontweight', 'bold', 'FontAngle', 'italic');

set(findall(gcf, '-property', 'fontsize'), 'fontsize', 10)

annotation('textbox', [0.038, 0.95, 0, 0], 'string', 'A', 'Fontsize', 10, 'Fontweight', 'bold')
annotation('textbox', [0.53, 0.95, 0, 0], 'string', 'B', 'Fontsize', 10, 'Fontweight', 'bold')
annotation('textbox', [0.038, 0.73, 0, 0], 'string', 'C', 'Fontsize', 10, 'Fontweight', 'bold')
annotation('textbox', [0.53, 0.73, 0, 0], 'string', 'D', 'Fontsize', 10, 'Fontweight', 'bold')
annotation('textbox', [0.038, 0.505, 0, 0], 'string', 'E', 'Fontsize', 10, 'Fontweight', 'bold')
annotation('textbox', [0.53, 0.505, 0, 0], 'string', 'F', 'Fontsize', 10, 'Fontweight', 'bold')
annotation('textbox', [0.038, 0.285, 0, 0], 'string', 'G', 'Fontsize', 10, 'Fontweight', 'bold')
annotation('textbox', [0.53, 0.285, 0, 0], 'string', 'H', 'Fontsize', 10, 'Fontweight', 'bold')

%% Fig. S5: All three modifications

% Run model with all three modifications at once (seasonal supply,
% temperature-dependent uptake (Q10=2), and doubled mortality in the four
% surface boxes
[t, y, PE, PO, PD] = wrap_boxmodel('seasonal', 'Q10', 2, 'rmort_var', [2 2 2 2 1 1 1]);

% Plot concentrations over the last 5 years
plot_boxmodel(t, y, PE, PO, PD, 'shorttime', 5)

% Save time-averaged mean
D.allmods = mean(y(t>90, PE.Jdom));
B.allmods = mean(y(t>90, PE.Jbac));

% Plot concentrations across the boxes (spatially)
plot_boxmodel(t, y, PE, PO, PD, 'boxes')
p = findobj(gcf,'Type','patch');
p(2).CData = D.allmods;
p(1).CData = B.allmods;
a = findobj(gcf, 'type', 'axes');
caxis(a(2), [22 100])
cbar = findall(gcf, 'Tag', 'Colorbar');
cbar(2).Ticks = 25:25:100;

%% Table S3: Turnover time of DOC and water masses

% Turnover time of DOC in equilibrium
DOC_loss = (PD.nsubs./PD.numM).*(PD.r_max./PD.K).*PD.Bstar.*PD.Dstar; % [mmolC m-3 d-1] This is the uptake in equilibrium
turnover_DOC = PD.Dstar./DOC_loss; % [d]
turnover_DOC_y = turnover_DOC/365; % [y]

% Turnover time of water in the boxes
flux = PO.Psi_Sv;
water_loss = -flux(logical(eye(size(flux)))); % water lost from box [m3 s-1]
turnover_water = PO.V'./water_loss; % [s];
turnover_water_y = turnover_water/(60*60*24*365); % [y]


