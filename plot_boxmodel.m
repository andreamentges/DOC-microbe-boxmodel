%% Plotting function for boxmodel output
%
% Optional input argument (varargin) can be one of two strings:
% - 'time': returns a plot of DOC and biomass over time
% - 'boxes': returns a plot of DOC and biomass color-coded in the seven
%           ocean boxes

function [cx, cy] = plot_boxmodel(t, y, PE, PO, PD, varargin)

% calendar year
ty = t/(365); 

% get box coordinates
if PE.nb == 7
    % smallest box is square, respectively set width of all boxes
    W = repmat(sqrt(min(PO.A)), 1, 7); % width of all boxes [m]

    % calculate length of all boxes from their width
    L = PO.A./W; % length of all boxes [m]
    H = PO.H; % height of all boxes [m]
    assert(all(W.*L - PO.A==0), 'wrong width or length')

    % Try 2: NADW
    V_NADW = PO.V(6);
    h = (H(3)+H(5))-H(2);
    hx = (V_NADW-L(2)*W(2)*h-L(4)*W(4)*h)/((L(2)+L(3)+L(4))*W(2));
    NADW_bottom = 1000+hx; % NADW bottom depth

    % AABW
    V_AABW = PO.V(7);
    hy = (V_AABW-L(1)*W(1)*(750+hx))/(L(1)*W(1)+(L(2)+L(3)+L(4))*W(2));
    AABW_bottom = 1000+hx+hy;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Set up coordinates of boxes 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % starting at lower left corner, in clockwise direction.
    %(replicate last point if < 8 coordinates so all vectors are same size)
    AAx = [0 0 L(1) L(1) L(1) L(1) L(1) L(1)];
    AAy = [H(1) 0 0 H(1) H(1) H(1) H(1) H(1)];
    SAx = [L(1) L(1) L(1)+L(2) L(1)+L(2) L(1)+L(2) L(1)+L(2) L(1)+L(2)...
        L(1)+L(2)];
    SAy = [H(2) 0 0 H(2) H(2) H(2) H(2) H(2)];
    LLx = [L(1)+L(2) L(1)+L(2) sum(L(1:3)) sum(L(1:3)) sum(L(1:3)) ...
        sum(L(1:3)) sum(L(1:3)) sum(L(1:3))];
    LLy = [H(3) 0 0 H(3) H(3) H(3) H(3) H(3)];
    NAx = [sum(L(1:3)) sum(L(1:3)) sum(L(1:4)) sum(L(1:4)) sum(L(1:4)) ...
        sum(L(1:4)) sum(L(1:4)) sum(L(1:4))];
    NAy = [H(4) 0 0 H(4) H(4) H(4) H(4) H(4)];
    TCx = [L(1)+L(2) L(1)+L(2) sum(L(1:3)) sum(L(1:3)) sum(L(1:3)) ...
        sum(L(1:3)) sum(L(1:3)) sum(L(1:3))];
    TCy = [H(3)+H(5) H(3) H(3) H(3)+H(5) H(3)+H(5) H(3)+H(5) H(3)+H(5)...
        H(3)+H(5)];
    NADWx = [L(1) L(1) L(1)+L(2) L(1)+L(2) sum(L(1:3)) sum(L(1:3))...
        sum(L(1:4)) sum(L(1:4))];
    NADWy = [NADW_bottom H(2) H(2) H(3)+H(5) H(3)+H(5) H(4) H(4) NADW_bottom];
    AABWx = [0 0 L(1) L(1) sum(L(1:4)) sum(L(1:4)) sum(L(1:4)) sum(L(1:4))];
    AABWy = [AABW_bottom H(1) H(1) NADW_bottom NADW_bottom AABW_bottom...
         AABW_bottom AABW_bottom];
    cx = [AAx; SAx; LLx; NAx; TCx; NADWx; AABWx]';
    cy = [AAy; SAy; LLy; NAy; TCy; NADWy; AABWy]';
        
    elseif PE.nb == 2
        
        nb = PE.nb;
        W = sqrt(min(PO.A)); % width of all boxes [m]
        L = PO.A./W; % length of all boxes [m]
        H = PO.H; % height of all boxes [m]
        AAx = [0 0 L(1) L(1)];
        AAy = [H(1) 0 0 H(1)];
        TCx = [0 0 L(1) L(1)];
        TCy = [H(2) H(1) H(1) H(2)];
        cx = [AAx; TCx]';
        cy = [AAy; TCy]';
        
    elseif PE.nb == 1
        
        nb = PE.nb;
        W = sqrt(min(PO.A)); % width of all boxes [m]
        L = PO.A./W; % length of all boxes [m]
        H = PO.H; % height of all boxes [m]
        cx = [0 0 L(1) L(1)];
        cy = [H(1) 0 0 H(1)];
        
end

%% Plot development of all variables over time
if any(strcmp(varargin, 'time'))
    
    lw = 1;
    
    % change default colors
    co = PD.cols;
    set(groot, 'defaultAxesColorOrder', co, 'defaultAxesLineStyleOrder','-|--|:')

    figure('color', 'white', 'position', [444,425,515,412])
    subplot(3,2,1)
    plot(ty(ty<PD.yspin), y(ty<PD.yspin, PE.Jdom), 'linewidth', lw)
    ylabel(sprintf('DOC\n[mmol/m³]')), axis([0 inf 0 inf])
    title('spin-up phase')
    subplot(3,2,2)
    plot(ty(ty>=PD.yspin), y(ty>=PD.yspin, PE.Jdom), 'linewidth', lw)
    axis tight
    title('model run')
    
    subplot(3,2,3)
    plot(ty(ty<PD.yspin), y(ty<PD.yspin, PE.Jbac), 'linewidth', lw)
    ylabel(sprintf('Bacteria\n[mmol/m³]')), axis([0 inf 0 inf])
    subplot(3,2,4)
    plot(ty(ty>=PD.yspin), y(ty>=PD.yspin, PE.Jbac), 'linewidth', lw)
    axis tight

    l = legend(PD.BoxAbbr);
    l.Position = [0.77896,0.09,0.1257,0.1175];
    
      
    set(findall(gcf, '-property', 'fontsize'), 'fontsize', 10)
end

%% Plot over short time interval (last x years)
if any(strcmp(varargin, 'shorttime'))
    ind  = find(strcmp(varargin, 'shorttime'));
    
    if isnumeric(varargin{ind+1})
        tsp = varargin{ind+1};
    else
       tsp = 10; % plot last 10 years
    end
    
    t = t/365;
    
    lw = [1.3 1.3 1.3 1.3 1.5 1.7 2.1 2.5];
    lstyles = {'-', '-', '-', '-', '-', '-', '-'};
    
    idx = max(t)-tsp; % last tsp years of simulation
    figure('color', 'white', 'position', [368,497,1140,208])
    for i = 1:PE.nb
        subplot(1,3,1)
        hold on
        plot(t(t>=idx), y(t>=idx, PE.Jdom(i)), 'linewidth', lw(i),...
            'color', PD.cols(i,:), 'linestyle', lstyles{i})
        axis tight
        xlabel('Time [y]'), ylabel('DOC [mmolC/m³]') 
        ax1 = gca;

        subplot(1,3,2)
        hold on
        plot(t(t>=idx), y(t>=idx, PE.Jbac(i)), 'linewidth', lw(i),...
            'color', PD.cols(i,:), 'linestyle', lstyles{i})
        axis tight
        xlabel('Time [y]'), ylabel('Biomass [mmolC/m³]') 
        ax2 = gca;
    end
    ax1.YLim(1) = min(min(y(t>=idx, PE.Jdom)));
    ax2.YLim(1) = 0 ;
    ax1.Position(4) = ax1.Position(4)*0.8;
    ax2.Position(4) = ax2.Position(4)*0.8;
    l = legend('AA', 'SA', 'LL', 'NA', 'TC', 'NADW', 'AABW');
    set(findall(gcf, '-property', 'fontsize'), 'fontsize', 10)
    l.Position = [0.652,0.226,0.09574,0.69761];
    ax1.Position(2) = ax1.Position(2)+0.12;
    ax2.Position(2) = ax2.Position(2)+0.12;
    
    a = annotation('textbox', [0.112, 1, 0, 0], 'string', 'A', 'Fontsize', 10, 'Fontweight', 'bold');
    a.Position(2) = a.Position(2)+0.01;
    annotation('textbox', [0.395, 1, 0, 0], 'string', 'B', 'Fontsize', 10, 'Fontweight', 'bold')
    
end

%% Plot left column time and right column spatial end distribution

if any(strcmp(varargin, 'time_and_boxes'))
    
    lw = 1.5;
    
    % change default colors
    co = PD.cols;
    set(groot, 'defaultAxesColorOrder', co, 'defaultAxesLineStyleOrder','-|--|:')

    figure('color', 'white', 'position', [304,284,660,148])
    subplot(3,2,1)
    plot(ty, y(:, PE.Jdom), 'linewidth', lw)
    ylabel(sprintf('DOC\n[mmol/m³]')), axis([0 inf 0 inf])
    title('Model run over time')

    subplot(3,2,3)
    plot(ty, y(:, PE.Jbac), 'linewidth', lw)
    ylabel(sprintf('Bacteria\n[mmol/m³]')), axis([0 inf 0 inf])
    
    % Plot with colorscale, at end of model run
    DOM_conc = y(end, PE.Jdom);
    Bac_conc = y(end, PE.Jbac);
    
    subplot(3,2,2)
    set(gca, 'YDir', 'reverse', 'Color', 'white')
    patch(cx, cy, DOM_conc)
    ylabel('Depth [m]'), axis tight
    c = colorbar; ylabel(c, 'DOC [mmol/m³] ')
    colormap(parula(100))
    legend('off')
    ax = gca;
    ax.XTick = ax.XLim;
    ax.XTickLabel = {'South', 'North'};
    ax.YTick = [0, 1000, 2000, 3000];
    title('End distribution')
    
    subplot(3,2,4)
        set(gca, 'YDir', 'reverse', 'Color', 'white')
    patch(cx, cy, Bac_conc)
    ylabel('Depth [m]'), axis tight
    c = colorbar; ylabel(c, 'Biomass [mmol/m³] ')
    colormap(parula(100))
    legend('off')
    ax = gca;
    ax.XTick = ax.XLim;
    ax.XTickLabel = {'South', 'North'};
    ax.YTick = [0, 1000, 2000, 3000];
      
    set(findall(gcf, '-property', 'fontsize'), 'fontsize', 10)
end


%% Plot spatial distribution in boxes
% Assuming that box surface areas are rectangles, such that the smallest
% box is square
% Height of box 6 and 7 are mean heights
if any(strcmp(varargin, 'boxes'))    

    % Plot with colorscale, at end of model run
    DOM_conc = y(end, PE.Jdom);
    Bac_conc = y(end, PE.Jbac);
    
    if any(strcmp(varargin, 'DOC_only'))
        
        figure('color', 'white', 'position', [304,284,660,148])
        s = subplot(1,2,2);
        set(s, 'YDir', 'reverse', 'Color', 'white')
        patch(cx, cy, DOM_conc)
        % xlabel('Horizontal extent [m]')
        ylabel('Depth [m]'), axis tight
        c = colorbar; ylabel(c, 'DOC [mmol/m³] ')
        colormap(parula(100))
        % title('DOC after spin-up phase')
        legend('off')
        ax = gca;
        ax.XTick = ax.XLim;
        ax.XTickLabel = {'South', 'North'};
        ax.YTick = [0, 1000, 2000, 3000];
        
        set(findall(gcf, '-property', 'fontsize'), 'fontsize', 11.5)
        
    elseif any(strcmp(varargin, 'subplots'))
        
        % find subplot data
        ind  = find(strcmp(varargin, 'subplots'));
        values = varargin{ind+1};
        
        for i = 1:size(values,1)
            s = subplot(2,ceil(size(values,1)/2),i);
            set(s, 'YDir', 'reverse', 'Color', 'white')
            patch(cx, cy, values(i,:))
            colormap(parula(100))
            caxis([min(min(values)) max(max(values))])
            axis tight
            ax = gca;
            ax.XTick = ax.XLim;
            ax.XTickLabel = {'S', 'N'};
            if i == 1 || i == ceil(size(values,1)/2)+1
                ax.YTick = [0, 1000, 2000, 3000];
                ylabel('Depth [m]')
            elseif i == size(values,1)
                ax.YTick = [];
                c = colorbar; 
                c.Position = c.Position + [0.09 0 0 0];
                ylabel(c, sprintf('DOC\n[mmol/m³]'))
            else
                ax.YTick = [];
            end
        end   

    set(findall(gcf, '-property', 'fontsize'), 'fontsize', 11.5)
        
    elseif any(strcmp(varargin, 'all'))
        
        figure('color', 'white', 'position', [309,285,1066,159])
        subplot(1,3,1, 'YDir', 'reverse')
        patch(cx, cy, DOM_conc)
        ylabel('Depth [m]'), axis tight
        c = colorbar; ylabel(c, sprintf('DOC\n[mmol/m³] '))
        if length(c.Ticks) == 6
            c.Ticks([1,3,5]) = [];
        elseif length(c.Ticks) == 5
            c.Ticks([2,4]) =[];
        end
        colormap(parula(100))
        % title('DOC after spin-up phase')
        ax = gca;
        ax.XTick = ax.XLim;
        ax.XTickLabel = {'South', 'North'};
        ax.YTick = [0, 1000, 2000, 3000];

        subplot(1,3,2, 'YDir', 'reverse')
        patch(cx, cy, Bac_conc)
        %xlabel('Horizontal extent [m]')
        axis tight
        c = colorbar; ylabel(c, sprintf('Bacterial biomass\n[mmol/m³]'))
        colormap(parula(100))
        % title('Bacteria after spin-up phase')
        ax = gca;
        ax.XTick = ax.XLim;
        ax.XTickLabel = {'South', 'North'};
        ax.YTick = [];
        if length(c.Ticks) == 6
            c.Ticks([1,3,5]) = [];
        elseif length(c.Ticks) == 5
            c.Ticks([2,4]) =[];
        end
        
        set(findall(gcf, '-property', 'fontsize'), 'fontsize', 11.5)
        
    else
        
        figure('color', 'white', 'position', [304,284,660,148])
        subplot(1,2,1, 'YDir', 'reverse')
        patch(cx, cy, DOM_conc)
        %xlabel('Horizontal extent [m]')
        ylabel('Depth [m]'), axis tight
        c = colorbar; ylabel(c, sprintf('DOC\n[mmolC/m³]'))
        colormap(parula(100))
        caxis([22 65])
        c.Ticks = 30:10:70;
        % title('DOC after spin-up phase')
        ax = gca;
        ax.XTick = ax.XLim;
        ax.XTickLabel = {'South', 'North'};
        ax.YTick = [0, 1000, 2000, 3000];

        subplot(1,2,2, 'YDir', 'reverse')
        patch(cx, cy, Bac_conc)
        % xlabel('Horizontal extent [m]')
        axis tight
        c = colorbar; ylabel(c, sprintf('Bacterial biomass\n[mmolC/m³]'))
        c.Ticks = (0:0.5:2);
        colormap(parula(100))
        % title('Bacteria after spin-up phase')
        caxis([0 2])
        ax = gca;
        ax.XTick = ax.XLim;
        ax.XTickLabel = {'South', 'North'};
        ax.YTick = [];
        
        annotation('textbox', [0.06, 0.99, 0, 0], 'string', 'A', 'Fontsize', 11, 'Fontweight', 'bold')
        annotation('textbox', [0.52, 0.99, 0, 0], 'string', 'B', 'Fontsize', 11, 'Fontweight', 'bold')
        
        set(findall(gcf, '-property', 'fontsize'), 'fontsize', 11.5)
        
    end
    
    % Ratio deep to surface
    rsd1 = mean(DOM_conc(PO.Idp))/mean(DOM_conc(PO.Isfc));
    rsd2 = mean(DOM_conc(PO.Idp))/mean(DOM_conc(PO.Isfc));
    fprintf('\nmean DOC deep/surface = %1.2f (without TC = %1.2f)', rsd1, rsd2)
    
    set(findall(gcf, '-property', 'fontsize'), 'fontsize', 11.5)
    
end


end