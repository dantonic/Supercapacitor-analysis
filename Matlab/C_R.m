function [C, R, Rd, nl, Pd, E] = C_R(fname, inputFunctionHandle, I, cycleNo, noCycles, currentNo, noCurrents, LimUr, LimUc, pf1, pf2)
% C_R2 - raèuna R i C i plota grafove ako je plt == true

persistent lastCurrentNo;   % to determine when processing of next current started

n = inputFunctionHandle(fname,I);   % data
dataLen = size(n,1);


% Find voltage drop (discharge start)
idx = 0;
startIdx = find(n(:,1)>=0.5,1);     % start from 0.5s
endIdx = find(n(:,1)>=5,1);       % calculate self-discharge du/dt from interval[0.5,5]s
base_du_dt = (n(endIdx,2)-n(startIdx,2)) / (n(endIdx,1)-n(startIdx,1));
baseInterval = (n(endIdx,1)-n(startIdx,1)) / (endIdx - startIdx);
found = false;
while ~found
    for i = startIdx : dataLen-1
        if n(i+1,1)-n(i,1) < baseInterval / 2   % ignore points that are too close
            continue
        end
        du_dt = (n(i+1,2)-n(i,2)) / (n(i+1,1)-n(i,1));
        if du_dt < 5 * base_du_dt     % discharge, so find 'more negative'
            idx = i;
            Ur = n(idx,2);
            t0 = n(idx,1);          % discharge start index, voltage and time
            break
        end
    end
    if idx==0
        %error(['Unable to find voltage drop (cycle ' num2str(cycleNo) ')!']);
        idx = 1;
        Ur = n(idx,2);
        t0 = n(idx,1);          % discharge start index, voltage and time
		break
    end
    % check for the false drop (voltage should decrease after the step for at least 1s)
    endIdx = find(n(:,1)>=t0+1,1);
    if isempty(endIdx)  % if not found
        endIdx = dataLen-1;
    end
    for i = idx : endIdx
        if n(i+1,2) > n(i,2)    % voltage raise
            break
        end
    end
    if i == endIdx              % voltage continuously decreasing
        found = true;
    else
        startIdx = idx + 1;
    end
end

% Find end of discharge
for i = idx : dataLen-3
    if n(i+1,2) > n(i,2) && n(i+2,2) > n(i+1,2) && n(i+3,2) > n(i+2,2)   % voltage raise for at least three samples
        break
    end
end
dataLen = i;


% Interval for calculating ESR
if LimUr(1) == 1    % 3rd order polynomial approximation between t0 and LimUr(2) is used to calculate ESR
    idxRStart = idx+1;
    idxREnd = find(n(:,1)>=t0+LimUr(2),1);
    if isempty(idxREnd)  % if not found
        idxREnd = dataLen;
    end
    p = polyfit(n(idxRStart:idxREnd,1),n(idxRStart:idxREnd,2),3);   % 3rd order polynomial
else                % linear approximation between voltages LimUr(1) and LimUr(2) is used to calculate ESR
    idxRStart = find(n(:,2)<Ur*LimUr(1),1);
    if isempty(idxRStart)  % if not found
        error(['Resistance calculation start voltage too low (cycle ' num2str(cycleNo) ')!']);
    end
    idxREnd = find(n(:,2)<Ur*LimUr(2),1);
    if isempty(idxREnd)  % if not found
        idxREnd = dataLen;
    end
    p = polyfit(n(idxRStart:idxREnd,1),n(idxRStart:idxREnd,2),1);   % 1st order polynomial
end


% Interval for calculating C
if LimUc(1) == 1    % Uc1 is voltage drop at the discharge start time
    idxCStart = find(n(:,2)<Ur,1);   % find first data point having voltage < Ur
    
    if LimUc(2) == 0    % Uc2 is voltage of last recorded discharge curve sample
        idxCEnd = dataLen;
    else
        idxCEnd = find(n(:,2)<Ur*LimUc(2),1);
        if isempty(idxCEnd)  % if not found
            idxCEnd = dataLen;
        end
    end
    
else
    idxCStart = find(n(:,2)<Ur*LimUc(1),1);
    if isempty(idxCStart)  % if not found
        error(['Capacitance calculation start voltage too low (cycle ' num2str(cycleNo) ')!']);
    end
    
    idxCEnd = find(n(:,2)<Ur*LimUc(2),1);
    if isempty(idxCEnd)  % if not found
        idxCEnd = dataLen;
    end
end

UcStart = n(idxCStart,2);
UcEnd = n(idxCEnd,2);



% CAPACITANCE C = 2W / [(0.9ur)^2 - (0.7ur)^2], W - exchanged energy between UcStart and UcEnd
% Energy
W = 0;
for i = idxCStart+1 : idxCEnd
    W = W + abs(n(i,3)) * n(i,2) * (n(i,1)-n(i-1,1));    % dW = i * u * dt
end
C = 2 * W / (UcStart*UcStart - UcEnd*UcEnd);

% Self-discharge resistance
U0 = n(1,2);
Rd = - t0 / (C * log(Ur / U0));


% ESR
R = (Ur - polyval(p,t0)) / abs(n(idx+1,3)); % R = dU / I


% Maximum power density
Pd = 0.25 * Ur * Ur / R;


% Retrived energy
idxStart = find(n(:,2)<Ur,1);   % find first data point having voltage < Ur
idxEnd = dataLen;

E = 0;
for i = idxStart+1 : idxEnd
    E = E + abs(n(i,3)) * n(i,2) * (n(i,1)-n(i-1,1));    % dW = i * u * dt
end



% NONLINEARITY
% Ideal (linear) discharge curve
t1x = n(idx+1,1);
t1y = n(idx+1,2);

t2x = n(dataLen,1);
t2y = n(dataLen,2);

% Nonlinearity (mean square deviation from ideal discharge curve
k = (t2y - t1y) / (t2x - t1x);
nl = 0;
for i = idx + 1 : dataLen
    x = n(i,1);
    y = t1y + k * (x - t1x);
    du = y - n(i,2);
    nl = nl + du * du;
end
nl = nl / (dataLen - idx);
nl = sqrt(nl);

% PLOTTING
if (mod(cycleNo,pf1)) == 0
    % Discharge curve plot
    h = figure('Name', 'C and ESR', 'OuterPosition', [100,100,1000,600]);
    plot([t1x t2x],[t1y t2y],'r--');    % ideal (linear) discharge
    hold all
    plot(n(1:dataLen,1),n(1:dataLen,2),'b','LineWidth',2);
    
    % Limits
    xlim([0,n(dataLen,1)]);
    ylim([n(dataLen,2),1.05*n(1,2)]);
    
    % Markers
    plot([0 n(idxCStart,1)], [n(idxCStart,2) n(idxCStart,2)],'m')
    my = n(idxCEnd,2);
    yl = ylim;
    if my <= yl(1)
        my = yl(1)+0.01;
    end
    plot([0 n(idxCEnd,1)], [my my],'m')
    
    mEnd = n(idxREnd,1)+5;  % end of ESR markers
    if mEnd > n(dataLen,1)
        mEnd = n(dataLen,1);
    end
    plot([n(idxRStart,1) mEnd], [n(idxRStart,2) n(idxRStart,2)],'g')
    plot([n(idxREnd,1) mEnd], [n(idxREnd,2) n(idxREnd,2)],'g')
    
    plot([t0 t0], [polyval(p,t0) U0],':g','LineWidth',2)
    % discharge curve interpolation (linear or cubic)
    y = zeros(idxREnd-idx+1,1);
    for i = idx : idxREnd
        y(i-idx+1) = polyval(p,n(i,1));
    end
    plot(n(idx:idxREnd,1), y,':g','LineWidth',2)
    
    ha = gca;
    set(ha, 'Box', 'on');
    grid on
    xlabel('t[s]');
    ylabel('u[V]');
    
    % Oznaèavanje i spremanje slike
    h1 = annotation('textbox', [0.75,0.8,0.1,0.1],'String', ...
        ['Cycle: ' int2str(cycleNo) char(10) char(10) ...
        'C = ' num2str(C,'%6.3f') ' F' char(10) ...
        'ESR = ' num2str(R,'%6.3f') ' \Omega' char(10) ...
        'Rd = ' num2str(Rd,'%6.0f') ' \Omega']);
    set(h1,'FontSize',12,'BackgroundColor',[0.9 0.9 0.9],'FitBoxToText','on');
    saveas(h,strtok(fname,'.'),'emf')
    saveas(h,strtok(fname,'.'),'fig')
    close(h)
end

% Cumulative discharge plot
if pf2>0 && mod(cycleNo,pf2) == 0
    noColors = noCycles / pf2;
    col = (cycleNo / pf2 - 1) * 0.75 / noColors;  % R and G from 0 to 0.75
    
    g = findobj('type','figure','Name', '_discharge');
    if isempty(g)   % first plot, set up figure
        h = figure('Name', '_discharge', 'OuterPosition', [100,100,1000,600]);
        plot(n(1:dataLen,1),n(1:dataLen,2),'Color',[col col 1],'LineWidth',2,'DisplayName',int2str(cycleNo));
        hold all
        ha = gca;
        set(ha, 'Box', 'on');
        grid on
        xlabel('t[s]');
        ylabel('u[V]');
    else            % plot already exists
        h = gcf;
        plot(n(1:dataLen,1),n(1:dataLen,2),'Color',[col col 1],'LineWidth',2,'DisplayName',int2str(cycleNo));
        hold all
    end
    ylim([n(dataLen,2),1.05*n(1,2)]);
    
    hLin = findobj(gca,'Type','line','LineStyle','-');
    hLin = flipud(hLin);
    legend(hLin);
end


% Cumulative discharge plot by currents - for last cycle
if isempty(lastCurrentNo)
    lastCurrentNo = 0;
end

if noCurrents > 1 && currentNo ~= lastCurrentNo && cycleNo == noCycles
    lastCurrentNo = currentNo;
    noColors = noCurrents;
    col = (currentNo - 1) * 0.75 / noColors;  % R and G from 0 to 0.75
    
    g = findobj('type','figure','Name', '_discharge_I');
    if isempty(g)   % first plot, set up figure
        h = figure('Name', '_discharge_I', 'OuterPosition', [100,100,1000,600]);
        plot(n(1:dataLen,1),n(1:dataLen,2),'Color',[col col 1],'LineWidth',2,'DisplayName',[num2str(I*1000) ' mA']);
        hold all
        ha = gca;
        set(ha, 'Box', 'on');
        grid on
        xlabel('t[s]');
        ylabel('u[V]');
    else            % plot already exists
        h = gcf;
        plot(n(1:dataLen,1),n(1:dataLen,2),'Color',[col col 1],'LineWidth',2,'DisplayName',[num2str(I*1000) ' mA']);
        hold all
    end
    ylim([n(dataLen,2),1.05*n(1,2)]);
    
    hLin = findobj(gca,'Type','line','LineStyle','-');
    hLin = flipud(hLin);
    legend(hLin);
end

