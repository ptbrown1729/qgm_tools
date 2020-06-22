Folders = [21,22,23,24,25,26];
BinEdges = [0,5,10,15,20,25,30,35];
BinCenters = 0.5*(BinEdges(2:end)+BinEdges(1:end-1));
NBins = length(BinEdges)-1;
RampTs = [200,300,400,500,600,700]; %ns

%correlation lengths
CorrLen = zeros(length(Folders),NBins);
Err = zeros(length(Folders),NBins);
%another bin

%domain sizes
DomSize = zeros(length(Folders),1);
DomPts = zeros(length(Folders),1);
DomErr = zeros(length(Folders),1);

AllContainers = {};

for ii = 1:length(Folders)
    df = DataFolder({2017 06 15 Folders(ii) 1 1},[],BinEdges);
    for jj = 1:NBins
        [Cl,Cle] = df.showCorrelationRange(jj);
        CorrLen(ii,jj) = Cl;
        Err(ii,jj) = Cle;
    end
    
    BinContainer = df.domainAnalysis(df.DistGrid,df.BinEdges);
    AllContainers{ii} = BinContainer;
    [DomainSizes,DomainSizesStd,DomainPts,~] = df.showDomainsVsBin(BinContainer);
    DomSize(ii) = DomainSizes(1);
    DomErr(ii) = DomainSizesStd(1);
    DomPts(ii) = DomainPts(1);
end

Exponents = zeros(NBins,1);
Uncs = zeros(NBins,1);
for kk = 1:NBins
%fit exponents
    InitP = [0.5,1];
    FixedP = [0,0];
    [Fp,Pn,FFH,SE] = fit1D(log(RampTs),log(CorrLen(:,kk)),[],{'line1D'},InitP,FixedP);
    Exponents(kk) = Fp(1);
    Uncs(kk) = SE(1);
    
    figure;
    loglog(1./RampTs,CorrLen(:,kk),'bo')
    hold on;
    loglog(1./RampTs,exp(FFH(log(RampTs))))
    grid on;
    xlabel('Rate')
    ylabel('Corr Length')
    title(sprintf('Corr Length Vs. Rate, bin %d',kk))

    figure;
    errorbar(RampTs,CorrLen(:,kk),Err(:,kk),'bo')
    hold on;
    InterpT = linspace(0,max(RampTs),300);
    plot(InterpT,exp(FFH(log(InterpT))),'b');
    % plot(InterpT,exp(line1D([0.5,1],log(InterpT))),'r');
    grid on;
    xlabel('Ramp Time (ns)')
    ylabel('Correlation Length (sites)')
    title(sprintf('Corr Length vs. Tquench Exponent = %0.2f +/- %0.2f, Bin %d',Fp(1),SE(1),kk));
end

figure;
errorbar(BinCenters,Exponents,Uncs,'bo');
grid on;
ylim([0,1])
xlabel('Bin Center')
ylabel('Exponent')

figure;
errorbar(df.Occs_AzAvg,Exponents,Uncs,'bo');
grid on;
ylim([0,1])
xlabel('Average Occupation')
ylabel('Exponent')
