function SCanalysis(Ur, I, m, S, text1, text2, text3, text4, text5, files, dataFunct, LimUr, LimUc, pf1,pf2)
%SCanalysis Calculate supercapacitor capacitance and equivalent serial resistance
%
%SCanalysis(I,Ur,m,S,'text1','text2','text3','text4','text5',files,[Ur1 Ur2],[Uc1 Uc2],pf1,pf2)
%
%Parameters:
%   Ur[V] - Nominal voltage (to which SC is charged)
%   I [mA] - Discharge current
%       a)	Scalar (e.g.21.5) Analysis for single current is performed
%       b)	Vector (e.g.[10 17.8 31.6 56.2 100]) Performs multi cycle - multi current analysis.
%   m[mg] - Mass of active material. Should be zero if not used.
%   S[cm2] - Area of one electrode. Should be zero if not used.
%   text1..text5 - Five comment lines that will be printed at plots
%   files - File name filter.Files that match the filter are sorted according to the file name.
%           Each file should contain data for one discharge cycle.
%           Care should be taken to ensure that sorting order corresponds to the cycles order,
%           e.g.name files as '17_04_25_Data_001.txt' ... '17_04_25_Data_100.txt'.In this case
%           'fnames' field should be set to *Data*.txt,where '*' matches any number of characters.
%           If vector of discharge currents is defined,files field should be defined as cell array,
%           having same number of elements as the I vector (and optionally one additional element),
%           where each element defines file name filter for corresponding current measurements,
%               e.g.: {‘I_15mA*.txt’ ‘I_22mA*.txt’}
%           If additional file is defined, it should contain repeated measurement for the first current.
%           Data from that file is shown with red cross at plots.
%   dataFunct - Handle of function for reading data from file. Function should accept two parameters,
%           file name and discharge current. If data contains measured current, discharge current
%           could be discarded. Function should return N-by-3 array, where first column contains
%           time in seconds, second voltage in volts and third current in amperes.
%   [Ur1 Ur2] - Resistance calculation start and end voltage.Voltage is specified as a percentage of
%           nominal voltage Ur,e.g.[0.9 0.7]. Straight-line approximation is applied to the discharge
%           curve between Ur1 and Ur2. Equivalent serial resistance is calculated from the voltage drop
%           at the discharge start time, which is the determined from the value of the straight line
%           at the discharge start time.
%           Special cases:
%               [1 t] - Ur1 is voltage at first sample after the discharge current was applied;
%                       Ur2 is voltage t seconds later.  3rd order polynomial approximation is used to
%                       determine the voltage drop at the discharge start time.
%   [Uc1 Uc2] - Capacitance calculation start and end voltage.Voltage is specified as a percentage
%           of nominal voltage Ur,e.g.[0.9 0.7].
%           Special cases:
%               [1 Uc2] - Uc1 is voltage at first sample after the discharge current was applied
%               [1 0] - Uc1 is voltage drop at the discharge start time,Uc2 is voltage of last
%                       recorded discharge curve sample
%   pf1 - Plot Frequency. Discharge curve plot will be generated each pF1 cycles.
%   pf2 - Plot Frequency. Cumulative discharge curve plot will be generated and discharge curve will be
%         plotted each pf2 cycles.
%
%Stored values (in file 'Results.mat'):
%   C - Capacitance
%   Cs_m - Specific capacitance per mass
%   Cs_a - Specific capacitance per area
%   R - Equivalent serial resistance
%   Rd - Discharge resistance,calculated from self-discharge curve during initial 5s rest period
%   nl - Discharge curve nonlinearity,calculated as mean square deviation from the ideal (linear)
%        discharge curve.
%   Pd - Specific power density, per mass.
%   E  - Specific energy, per mass.
%


% Attributes validation
validateattributes(I,{'double'},{'positive','row','nonempty'},'SCanalysis','I',1);
noCurrents = size(I,2);
validateattributes(Ur,{'double'},{'positive','scalar','nonempty'},'SCanalysis','Ur',2);
validateattributes(m,{'double'},{'nonnegative','scalar','nonempty'},'SCanalysis','m',3);
validateattributes(S,{'double'},{'nonnegative','scalar','nonempty'},'SCanalysis','S',4);

validateattributes(text1,{'char'},{},'SCanalysis','text1',5);
validateattributes(text2,{'char'},{},'SCanalysis','text2',6);
validateattributes(text3,{'char'},{},'SCanalysis','text3',7);
validateattributes(text4,{'char'},{},'SCanalysis','text4',8);
validateattributes(text5,{'char'},{},'SCanalysis','text5',9);

validateattributes(files,{'char','cell'},{'nonempty'},'SCanalysis','files',10);
if iscell(files)    % check number of elements and that each ellement is char
    noFileGroups = size(files,2);
    if noFileGroups~=noCurrents && noFileGroups~=noCurrents+1
        error('Expected input number 10, files, size should be equal to number of currents or number of currents + 1.');
    end
    for i = 1 : noFileGroups
        if ~ischar(files{i})
            error('Expected input number 10, files, all elements should be strings.');
        end
    end
else
    noFileGroups = 1;
end
validateattributes(dataFunct,{'function_handle'},{'nonempty'},'SCanalysis','dataFunct',11);

validateattributes(LimUr,{'double'},{'positive','size',[1,2]},'SCanalysis','LimUr',12);
if LimUr(1) < 1 && LimUr(1) <= LimUr(2) error('Expected input number 12, LimUr, first element should be greater than second.'); end
validateattributes(LimUc,{'double'},{'nonnegative','size',[1,2],'<=',1'},'SCanalysis','LimUc',13);
if LimUc(1) <= LimUc(2) error('Expected input number 13, LimUc, first element should be greater than second.'); end

validateattributes(pf1,{'double'},{'positive','scalar','nonempty'},'SCanalysis','pf1',14);
validateattributes(pf2,{'double'},{'nonnegative','scalar','nonempty'},'SCanalysis','pf2',15);


% Determine number of cycles
if iscell(files)
    f = dir(files{1});
else
    f = dir(files);
end
noCycles = size(f,1);
if noCycles == 0
    error('No data files selected, check current folder and input number 10, files.');
end

% Prealocate result arrays
C = zeros(noCycles,noFileGroups);
Cs_m = zeros(noCycles,noFileGroups);
Cs_a = zeros(noCycles,noFileGroups);
R = zeros(noCycles,noFileGroups);
Rd = zeros(noCycles,noFileGroups);
nl = zeros(noCycles,noFileGroups);
Pd = zeros(noCycles,noFileGroups);
E = zeros(noCycles,noFileGroups);

% Currents
for fileIdx = 1 : noFileGroups
    if fileIdx <= noCurrents	% if noFileGroups == noCurrents + 1, then measurement is repeated for the first current (repeatability) 
        I0 = I(fileIdx);
        fileSuffix = [num2str(I0) '_mA'];
        additionalCurrent = false;
    else
        I0 = I(1);
        fileSuffix = [num2str(I0) '_mA_2nd'];
        additionalCurrent = true;
    end
    % Replace '.' with '_'
    fileSuffix = strrep(fileSuffix,'.','_');
	
    % File list
    if iscell(files)
        f = dir(files{fileIdx});
    else
        f = dir(files);
    end
    if size(f,1) ~= noCycles
        error(['All currents should have the same number of cycles (' num2str(I0) 'mA: ' ...
            num2str(size(f,1)) ', should be ' num2str(noCycles)]);
    end
    fn = sort({f.name});    % sorted cell array of file names
    
    
    for n = 1 : noCycles
        [c1, r1, rd1, nl1, pd, e] = C_R(cell2mat(fn(n)),dataFunct,I0/1000,n,noCycles,fileIdx,noCurrents,LimUr,LimUc,pf1,pf2); % I !!!
        C(n,fileIdx) = c1;
        R(n,fileIdx) = r1;
        Rd(n,fileIdx) = rd1;
        nl(n,fileIdx) = nl1;
        Pd(n,fileIdx) = pd;
        E(n,fileIdx) = e;
    end
	if m > 0
		Cs_m(:,fileIdx) = C(:,fileIdx) ./ (m ./ 1000);
		Pd(:,fileIdx) = Pd(:,fileIdx) ./ (m ./ 1000);
		E(:,fileIdx) = E(:,fileIdx) ./ (m ./ 1000);
	end
    if S > 0
		Cs_a(:,fileIdx) = C(:,fileIdx) ./ S;
	end
    
    baseX = 0.7;
    baseY = 0.4;
    fontSize = 11;
    
    h = figure('Name', 'Capacitance', 'OuterPosition', [100,100,1000,600]);
    plot(C(:,fileIdx),'b','LineWidth',2);
    yl = ylim;
    yl(1) = 0; yl(2) = 1.05 * yl(2);
    ylim(yl);
    ha = gca;
    set(ha, 'Box', 'on');
    grid on
    xlabel('Cycle');
    ylabel('C[F]');
    annotate(baseX,baseY,fontSize,I0,Ur,m,S,text1,text2,text3,text4,text5);
    saveas(h,['_C_' fileSuffix],'emf')
    saveas(h,['_C_' fileSuffix],'fig')
    if noCurrents > 1 close(h); end
    
    h = figure('Name', 'ESR', 'OuterPosition', [100,100,1000,600]);
    plot(R(:,fileIdx),'b','LineWidth',2);
    yl = ylim;
    yl(1) = 0; yl(2) = 1.05 * yl(2);
    ylim(yl);
    ha = gca;
    set(ha, 'Box', 'on');
    grid on
    xlabel('Cycle]');
    ylabel('ESR[\Omega]');
    annotate(baseX,baseY,fontSize,I0,Ur,m,S,text1,text2,text3,text4,text5);
    saveas(h,['_ESR_' fileSuffix],'emf')
    saveas(h,['_ESR_' fileSuffix],'fig')
    if noCurrents > 1 close(h); end
    
    h = figure('Name', 'Self-discharge resistance', 'OuterPosition', [100,100,1000,600]);
    plot(Rd(:,fileIdx),'b','LineWidth',2);
    yl = ylim;
    yl(1) = 0; yl(2) = 1.05 * yl(2);
    ylim(yl);
    ha = gca;
    set(ha, 'Box', 'on');
    grid on
    xlabel('Cycle');
    ylabel('Rd[\Omega]');
    annotate(baseX,baseY,fontSize,I0,Ur,m,S,text1,text2,text3,text4,text5);
    saveas(h,['_Rd_' fileSuffix],'emf')
    saveas(h,['_Rd_' fileSuffix],'fig')
    if noCurrents > 1 close(h); end
    
    h = figure('Name', 'Discharge curve non-linearity', 'OuterPosition', [100,100,1000,600]);
    plot(nl(:,fileIdx),'b','LineWidth',2);
    yl = ylim;
    yl(1) = 0; yl(2) = 1.05 * yl(2);
    ylim(yl);
    ha = gca;
    set(ha, 'Box', 'on');
    grid on
    xlabel('Cycle');
    ylabel('MSE');
    annotate(baseX,baseY,fontSize,I0,Ur,m,S,text1,text2,text3,text4,text5);
    saveas(h,['_Nl_' fileSuffix],'emf')
    saveas(h,['_Nl_' fileSuffix],'fig')
    if noCurrents > 1 close(h); end
    
    if m > 0
        h = figure('Name', 'Specific capacitance (F/g)', 'OuterPosition', [100,100,1000,600]);
        plot(Cs_m(:,fileIdx),'b','LineWidth',2);
        yl = ylim;
        yl(1) = 0; yl(2) = 1.05 * yl(2);
        ylim(yl);
        ha = gca;
        set(ha, 'Box', 'on');
        grid on
        xlabel('Cycle');
        ylabel('Cs_m[F/g]');
        annotate(baseX,baseY,fontSize,I0,Ur,m,S,text1,text2,text3,text4,text5);
        saveas(h,['_Cs_m_' fileSuffix],'emf')
        saveas(h,['_Cs_m_' fileSuffix],'fig')
        if noCurrents > 1 close(h); end
        
        h = figure('Name', 'Maximum power density (W/g)', 'OuterPosition', [100,100,1000,600]);
        plot(Pd(:,fileIdx),'b','LineWidth',2);
        yl = ylim;
        yl(1) = 0; yl(2) = 1.05 * yl(2);
        ylim(yl);
        ha = gca;
        set(ha, 'Box', 'on');
        grid on
        xlabel('Cycle');
        ylabel('Pdm[W/g]');
        annotate(baseX,baseY,fontSize,I0,Ur,m,S,text1,text2,text3,text4,text5);
        saveas(h,['_Pdm_' fileSuffix],'emf')
        saveas(h,['_Pdm_' fileSuffix],'fig')
        if noCurrents > 1 close(h); end
        
        h = figure('Name', 'Energy density (J/g)', 'OuterPosition', [100,100,1000,600]);
        plot(E(:,fileIdx),'b','LineWidth',2);
        yl = ylim;
        yl(1) = 0; yl(2) = 1.05 * yl(2);
        ylim(yl);
        ha = gca;
        set(ha, 'Box', 'on');
        grid on
        xlabel('Cycle');
        ylabel('Ed[J/g]');
        annotate(baseX,baseY,fontSize,I0,Ur,m,S,text1,text2,text3,text4,text5);
        saveas(h,['_Ed_' fileSuffix],'emf')
        saveas(h,['_Ed_' fileSuffix],'fig')
        if noCurrents > 1 close(h); end
    end
    
    if S > 0
        h = figure('Name', 'Specific capacitance (/cm2)', 'OuterPosition', [100,100,1000,600]);
        plot(Cs_a(:,fileIdx),'b','LineWidth',2);
        yl = ylim;
        yl(1) = 0; yl(2) = 1.05 * yl(2);
        ylim(yl);
        ha = gca;
        set(ha, 'Box', 'on');
        grid on
        xlabel('Cycle');
        ylabel('Cs_a[F/cm^2]');
        annotate(baseX,baseY,fontSize,I0,Ur,m,S,text1,text2,text3,text4,text5);
        saveas(h,['_Cs_a_' fileSuffix],'emf')
        saveas(h,['_Cs_a_' fileSuffix],'fig')
        if noCurrents > 1 close(h); end
    end
    
    % Cumulative discharge curves
    g = findobj('type','figure','Name', '_discharge');
    if ~isempty(g)   % plot exists, save and close
        saveas(g,['_discharge_' fileSuffix],'emf')
        saveas(g,['_discharge_' fileSuffix],'fig')
        if noCurrents > 1 close(g); end
    end    
end


% CURRENT INFLUENCE ANALYSIS
if noCurrents > 1
    baseX = 0.7;
    baseY = 0.4;
    fontSize = 11;
    
    h = figure('Name', 'Capacitance', 'OuterPosition', [100,100,1000,600]);
    plot(I, C(noCycles,1:noCurrents),'b','LineWidth',2);
    if additionalCurrent
        hold on
        plot(I(1),C(noCycles,noFileGroups),'rx','LineWidth',2);
    end
    yl = ylim;
    yl(1) = 0; yl(2) = 1.05 * yl(2);
    ylim(yl);
    ha = gca;
    set(ha, 'Box', 'on');
    grid on
    xlabel('I[mA]');
    ylabel('C[F]');
    annotate(baseX,baseY,fontSize,0,Ur,m,S,text1,text2,text3,text4,text5);
    saveas(h,'_C','emf')
    saveas(h,'_C','fig')
    
    h = figure('Name', 'ESR', 'OuterPosition', [100,100,1000,600]);
    plot(I, R(noCycles,1:noCurrents),'b','LineWidth',2);
    if additionalCurrent
        hold on
        plot(I(1),R(noCycles,noFileGroups),'rx','LineWidth',2);
    end
    yl = ylim;
    yl(1) = 0; yl(2) = 1.05 * yl(2);
    ylim(yl);
    ha = gca;
    set(ha, 'Box', 'on');
    grid on
    xlabel('I[mA]');
    ylabel('ESR[\Omega]');
    annotate(baseX,baseY,fontSize,0,Ur,m,S,text1,text2,text3,text4,text5);
    saveas(h,'_ESR','emf')
    saveas(h,'_ESR','fig')
    
    h = figure('Name', 'Self-discharge resistance', 'OuterPosition', [100,100,1000,600]);
    plot(I, Rd(noCycles,1:noCurrents),'b','LineWidth',2);
    if additionalCurrent
        hold on
        plot(I(1),Rd(noCycles,noFileGroups),'rx','LineWidth',2);
    end
    yl = ylim;
    yl(1) = 0; yl(2) = 1.05 * yl(2);
    ylim(yl);
    ha = gca;
    set(ha, 'Box', 'on');
    grid on
    xlabel('I[mA]');
    ylabel('Rd[\Omega]');
    annotate(baseX,baseY,fontSize,0,Ur,m,S,text1,text2,text3,text4,text5);
    saveas(h,'_Rd','emf')
    saveas(h,'_Rd','fig')
    
    h = figure('Name', 'Discharge curve non-linearity', 'OuterPosition', [100,100,1000,600]);
    plot(I, nl(noCycles,1:noCurrents),'b','LineWidth',2);
    if additionalCurrent
        hold on
        plot(I(1),nl(noCycles,noFileGroups),'rx','LineWidth',2);
    end
    yl = ylim;
    yl(1) = 0; yl(2) = 1.05 * yl(2);
    ylim(yl);
    ha = gca;
    set(ha, 'Box', 'on');
    grid on
    xlabel('I[mA]');
    ylabel('MSE');
    annotate(baseX,baseY,fontSize,I0,Ur,m,S,text1,text2,text3,text4,text5);
    saveas(h,'_Nl','emf')
    saveas(h,'_Nl','fig')
    
    if m > 0
        h = figure('Name', 'Specific capacitance (F/g)', 'OuterPosition', [100,100,1000,600]);
        plot(I, Cs_m(noCycles,1:noCurrents),'b','LineWidth',2);
        if additionalCurrent
            hold on
            plot(I(1),Cs_m(noCycles,noFileGroups),'rx','LineWidth',2);
        end
        yl = ylim;
        yl(1) = 0; yl(2) = 1.05 * yl(2);
        ylim(yl);
        ha = gca;
        set(ha, 'Box', 'on');
        grid on
        xlabel('I[mA]');
        ylabel('Cs_m[F/g]');
        annotate(baseX,baseY,fontSize,I0,Ur,m,S,text1,text2,text3,text4,text5);
        saveas(h,'_Cs_m','emf')
        saveas(h,'_Cs_m','fig')
        
        h = figure('Name', 'Maximum power density (W/g)', 'OuterPosition', [100,100,1000,600]);
        plot(I, Pd(noCycles,1:noCurrents),'b','LineWidth',2);
        if additionalCurrent
            hold on
            plot(I(1),Pd(noCycles,noFileGroups),'rx','LineWidth',2);
        end
        yl = ylim;
        yl(1) = 0; yl(2) = 1.05 * yl(2);
        ylim(yl);
        ha = gca;
        set(ha, 'Box', 'on');
        grid on
        xlabel('I[mA]');
        ylabel('Pdm[W/g]');
        annotate(baseX,baseY,fontSize,I0,Ur,m,S,text1,text2,text3,text4,text5);
        saveas(h,'_Pdm','emf')
        saveas(h,'_Pdm','fig')
        
        h = figure('Name', 'Energy density (J/g)', 'OuterPosition', [100,100,1000,600]);
        plot(I, E(noCycles,1:noCurrents),'b','LineWidth',2);
        if additionalCurrent
            hold on
            plot(I(1),E(noCycles,noFileGroups),'rx','LineWidth',2);
        end
        yl = ylim;
        yl(1) = 0; yl(2) = 1.05 * yl(2);
        ylim(yl);
        ha = gca;
        set(ha, 'Box', 'on');
        grid on
        xlabel('I[mA]');
        ylabel('Ed[J/g]');
        annotate(baseX,baseY,fontSize,I0,Ur,m,S,text1,text2,text3,text4,text5);
        saveas(h,'_Ed','emf')
        saveas(h,'_Ed','fig')
    end
    
    if S > 0
        h = figure('Name', 'Specific capacitance (/cm2)', 'OuterPosition', [100,100,1000,600]);
        plot(I, Cs_a(noCycles,1:noCurrents),'b','LineWidth',2);
        if additionalCurrent
            hold on
            plot(I(1),Cs_a(noCycles,noFileGroups),'rx','LineWidth',2);
        end
        yl = ylim;
        yl(1) = 0; yl(2) = 1.05 * yl(2);
        ylim(yl);
        ha = gca;
        set(ha, 'Box', 'on');
        grid on
        xlabel('I[mA]');
        ylabel('Cs_a[F/cm^2]');
        annotate(baseX,baseY,fontSize,I0,Ur,m,S,text1,text2,text3,text4,text5);
        saveas(h,'_Cs_a','emf')
        saveas(h,'_Cs_a','fig')
    end
    
    % Cumulative discharge curves - by current
    g = findobj('type','figure','Name', '_discharge_I');
    if ~isempty(g)   % plot exists, save
        saveas(g,'_discharge_I_','emf')
        saveas(g,'_discharge_I_','fig')
    end    
    
end
save ('Results.mat','C','Cs_m','Cs_a','R','Rd','nl','Pd','E')

function annotate(baseX,baseY,fontSize,I,U,m,S,t1,t2,t3,t4,t5)
Y = baseY;
txt = '';
mStxt = '';
if m > 0
    mStxt = ['m = ' num2str(m,'%4.1f') ' mg  '];
end
if S > 0
    mStxt = [mStxt 'S = ' num2str(S,'%4.1f') ' cm2'];
end
if m > 0 || S > 0
    txt = mStxt;
end
if I == 0
    currTxt = '';
else
    currTxt = ['I = ' num2str(I,'%6.1f') ' mA'];
end
txt = [txt char(10) t1 char(10) t2 char(10) t3 char(10) t4 char(10) t5 char(10) ...
    'U = ' num2str(U,'%3.1f') ' V  ' currTxt];

h1 = annotation('textbox', [baseX,baseY,0.1,0.1],'String', txt);
set(h1,'FontSize',fontSize,'BackgroundColor',[0.9 0.9 0.9],'FitBoxToText','on');



