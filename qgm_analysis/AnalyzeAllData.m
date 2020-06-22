MI_Dataset_List = {{2016 10 28 28 1 2};{2016 11 1 31 1 2};{2016 11 1 37 1 2};{2016 10 31 34 1 2};{2016 10 31 23 1 2};{2016 11 1 42 1 2};{2016 10 31 53 1 2}};%;{}};
SzSz_Blow1_List = {{2016 10 28 24 1 2};{2016 11 1 27 1 2};{2016 11 1 33 1 2};{2016 10 31 30 1 2};{2016 10 31 19 1 2};{2016 11 1 38 1 2};{2016 10 31 49 1 2}};%;{2016 11 11 24 1 2}};
SzSz_Blow2_List = {{2016 10 28 26 1 2};{2016 11 1 29 1 2};{2016 11 1 35 1 2};{2016 10 31 32 1 2};{2016 10 31 21 1 2};{2016 11 1 40 1 2};{2016 10 31 51 1 2}};%;{2016 11 11 23 1 2}};
SxSx_Blow1_List = {{2016 10 28 27 1 2};{2016 11 1 30 1 2};{2016 11 1 36 1 2};{2016 10 31 33 1 2};{2016 10 31 22 1 2};{2016 11 1 41 1 2};{2016 10 31 52 1 2}};%;{}};
SxSx_Blow2_List = {{2016 10 28 25 1 2};{2016 11 1 28 1 2};{2016 11 1 34 1 2};{2016 10 31 31 1 2};{2016 10 31 20 1 2};{2016 11 1 39 1 2};{2016 10 31 50 1 2}};%;{}};

%BinEdges = [0,sqrt(linspace(8^2,30^2,25))];
BinEdges = sqrt([0,linspace(8^2,30^2,30)]);
SavePath = fullfile('Data',sprintf('Data_%s',datestr(now,'yyyy-mm-dd;HH;MM')));
SaveName = 'Data.mat';
SaveFile = fullfile(SavePath,SaveName);

NSets = length(MI_Dataset_List);
PNominal = [0.0225,0.1810,0.3387,0.4783,0.6533,0.7687,0.9415]; %,0.3];
PUnc = [0.0417,0.0405,0.0437,0.0430,0.0506,0.0325,0.0242]; %,0];

UoverTNominal = [8*ones(1,NSets),14];
%CorrDSets = [];
DSummary = DataSummaryV2();
DSummary.PolarizationUnc = PUnc;

for ii = 1:NSets
    DSetNames = {MI_Dataset_List{ii};SzSz_Blow2_List{ii};SzSz_Blow1_List{ii};SxSx_Blow2_List{ii};SxSx_Blow1_List{ii}};
    CurrentCorrDSet = CorrDataSetV2(DSetNames,BinEdges);
    CurrentCorrDSet.shrinkToSave();
    CurrentCorrDSet.SavePath = SavePath;
    CurrentCorrDSet.GlobalPolarization = PNominal(ii);
    CurrentCorrDSet.UoverTs = UoverTNominal(ii);
    %CorrDSets = [CorrDSets,CurrentCorrDSet];
    CurrentCorrDSet.saveRadProfilesToTxt;

    DSummary.append(CurrentCorrDSet);
    
    %Create new name for saving
    VarName = sprintf('DSet%d',ii);
    eval(sprintf('%s = CurrentCorrDSet',VarName));
    
    if exist(SaveFile,'file')
        save(SaveFile,VarName,'-append');
    else
        save(SaveFile,VarName);
    end
    
    clear CurrentCorrDSet
    eval(sprintf('clear %s',VarName));
end

save(SaveFile,'DSummary','-append');

DSummary.SavePath = SavePath;
DSummary.showCorrVsP();
DSummary.showCorrMat();
DSummary.generateCorrMatTextFiles();
DSummary.showMaxMomentCorr();

DSummary.saveNNCorrs();
DSummary.saveMaxMomentCorr();


%%
%produce a few interesting plots
load Data2
DStack = [DSet1,DSet2,DSet3,DSet4,DSet5,DSet6,DSet7];
%colormap
colormap(jet(length(DStack)));
cmap = colormap;

%NN Spin Correlator vs. Density
figure('name','Spin Corr NN Vs. Density')
for ii = 1:length(DStack)
errorbar(DStack(ii).MI.Occs_AzAvg,squeeze(DStack(ii).SzSz_Corr(5,6,:)),squeeze(DStack(ii).SzSz_CorrUnc(5,6,:)),'-o','color',cmap(ii,:))
hold on;
end
hold off;
grid on;
xlabel('Singles Density/Moment')
ylabel('C^z_s(0,1)')

%NN spin correlator normalized to density Vs. Density
figure('name','Spin Corr NN Normalized to Density Vs. Density')
for ii = 1:length(DStack)
%errorbar(DStack(ii).MI.Occs_AzAvg,squeeze(DStack(ii).SzSz_Corr(5,6,:)),squeeze(DStack(ii).SzSz_CorrUnc(5,6,:)),'-o','color',cmap(ii,:))
%plot(DStack(ii).MI.Occs_AzAvg,squeeze(sum(sum(DStack(ii).SzSz_Corr(5,5,:),1),2))./DStack(ii).MI.Occs_AzAvg,'-o','color',cmap(ii,:))
plot(DStack(ii).MI.Occs_AzAvg,squeeze(sum(sum(DStack(ii).SzSz_Corr(5,5,:),1),2)),'-o','color',cmap(ii,:))

hold on;
end
hold off;
grid on;
xlabel('Singles Density/Moment')
ylabel('C^z_s(0,1)')
ylim([0,2.5])

%Spin susceptibility from fluctuation dissipation theorem
figure('name','Spin Susceptibility Vs. Density')
for ii = 1:length(DStack)
%errorbar(DStack(ii).MI.Occs_AzAvg,DStack(ii).SpinSusceptibility,DStack(ii).SpinSusceptibilityUnc,'-o','color',cmap(ii,:))
plot(DStack(ii).MI.Occs_AzAvg,DStack(ii).SpinSusceptibility./DStack(ii).MI.Occs_AzAvg,'-o','color',cmap(ii,:))
hold on;
end
hold off;
grid on;
xlabel('Singles Density/Moment')
ylabel('\chi')

%central density vs. polarization
figure('name','Singles Density Vs. Polarization')
errorbar(DSummary.Polarization,DSummary.MI_CentralDensity,DSummary.MI_CentralDensityUnc,'r-o');
hold on;
errorbar(DSummary.Polarization,DSummary.Sz_Up_CentralDensity,DSummary.Sz_Up_CentralDensityUnc,'b-o');
errorbar(DSummary.Polarization,DSummary.Sz_Down_CentralDensity,DSummary.Sz_Down_CentralDensityUnc,'g-o');
plot(DSummary.Polarization,DSummary.Sz_Up_CentralDensity+DSummary.Sz_Down_CentralDensity,'m-o');
grid on;
xlabel('Polarization')
ylabel('Central Singles Density')

%%
%analyze U/t=15 susceptibility data.
MI_Dataset_List = {{2016 11 14 23 1 2}};
SzSz_Blow1_List = {{2016 11 14 22 1 2}}; %{{2016 11 11 24 1 2}};
SzSz_Blow2_List = {{2016 11 14 21 1 2}}; %{{2016 11 11 23 1 2}};
SxSx_Blow1_List = {{}};
SxSx_Blow2_List = {{}};

% BinEdges = sqrt(linspace(0,34^2,15));
BinEdges = [0,9,11,13,16,18,20,22,24,26,27,28,29,30,31,32,33,34,35];

ii = 1;
DSetNames = {MI_Dataset_List{ii};SzSz_Blow2_List{ii};SzSz_Blow1_List{ii};SxSx_Blow2_List{ii};SxSx_Blow1_List{ii}};
CurrentCorrDSet = CorrDataSetV2(DSetNames,BinEdges);

CurrentCorrDSet.GlobalPolarization = 0.29; %+/- 0.3
CurrentCorrDSet.SavePath = fullfile('Data',sprintf('Data_%s',datestr(now,'yyyy-mm-dd;HH;MM')));
CurrentCorrDSet.saveRadProfilesToSingleTxt();

figure('name','Radial Profiles');
errorbar(CurrentCorrDSet.Up_Sz.BinAvg,CurrentCorrDSet.Up_Sz.Occs_AzAvg,CurrentCorrDSet.Up_Sz.Occs_AzAvgUnc,'ro-')
hold on;
errorbar(CurrentCorrDSet.Down_Sz.BinAvg,CurrentCorrDSet.Down_Sz.Occs_AzAvg,CurrentCorrDSet.Down_Sz.Occs_AzAvgUnc,'bo-')
errorbar(CurrentCorrDSet.Down_Sz.BinAvg,CurrentCorrDSet.Sz_LocalPolarization,CurrentCorrDSet.Sz_LocalPolarizationUnc,'go-')

figure('name','Pol Vs Density');
errorbar(CurrentCorrDSet.Up_Sz.Occs_AzAvg+CurrentCorrDSet.Down_Sz.Occs_AzAvg,CurrentCorrDSet.Sz_LocalPolarization,CurrentCorrDSet.Sz_LocalPolarizationUnc,'ro-')

%%
%comparison to debayan susceptibility.


%analyze U/t=15 susceptibility data.
MI_Dataset_List = {{2016 11 14 23 1 2}};
SzSz_Blow1_List = {{2016 11 14 22 1 2}}; %{{2016 11 11 24 1 2}};
SzSz_Blow2_List = {{2016 11 14 21 1 2}}; %{{2016 11 11 23 1 2}};
SxSx_Blow1_List = {{}};
SxSx_Blow2_List = {{}};

% BinEdges = sqrt(linspace(0,30^2,30));
BinEdges = [0,8,12,16,20,24,28];
ii = 1;
DSetNames = {MI_Dataset_List{ii};SzSz_Blow2_List{ii};SzSz_Blow1_List{ii};SxSx_Blow2_List{ii};SxSx_Blow1_List{ii}};
CurrentCorrDSet = CorrDataSetV2(DSetNames,BinEdges);

%CurrentCorrDSet.SzSz_Corr = 2*(CurrentCorrDSet.Up_Sz.Density_Corr_AzAvg + CurrentCorrDSet.Down_Sz.Density_Corr_AzAvg);
Susceptibility = squeeze(sum(sum(CurrentCorrDSet.SzSz_Corr,1),2));
%Density = CurrentCorrDSet.Up_Sz.Occs_AzAvg + CurrentCorrDSet.Down_Sz.Occs_AzAvg;
Density = CurrentCorrDSet.MI.Occs_AzAvg;

h = figure('name','Vs. Position');
errorbar(CurrentCorrDSet.BinCenters,CurrentCorrDSet.Sz_LocalPolarization,CurrentCorrDSet.Sz_LocalPolarizationUnc);
hold on;
errorbar(CurrentCorrDSet.Up_Sz.BinAvg,CurrentCorrDSet.Up_Sz.Occs_AzAvg,CurrentCorrDSet.Up_Sz.Occs_AzAvgUnc);
errorbar(CurrentCorrDSet.Down_Sz.BinAvg,CurrentCorrDSet.Down_Sz.Occs_AzAvg,CurrentCorrDSet.Down_Sz.Occs_AzAvgUnc);
errorbar(CurrentCorrDSet.Down_Sz.BinAvg,CurrentCorrDSet.Down_Sz.Occs_AzAvg+CurrentCorrDSet.Up_Sz.Occs_AzAvg,CurrentCorrDSet.Down_Sz.Occs_AzAvgUnc);
plot(CurrentCorrDSet.Up_Sz.BinAvg,Susceptibility,'b-o')
plot(CurrentCorrDSet.Up_Sz.BinAvg,Susceptibility./Density,'r-o')
grid on;
legend({'SpinUp','SpinDown','MI','Susceptibility','Susceptibility/Density'})
xlim([0,25]);

%comparision to debayan...
DebDensity=[0.8946,0.8523,0.7032,0.4915,0.2669,0.0857];
DebSusceptibility=[0.3067,0.2753,0.245,0.2019,0.115,0.0336];


h = figure('name','Susceptibility Vs. Density');
plot(Density,Susceptibility,'b-o');
hold on;
plot(DebDensity,DebSusceptibility,'r-o')
legend({'PB Susceptibility','Deb Susceptibility'})
xlabel('Density');


% h = figure('name','Susceptibility/Density Vs. Density');
% plot(Density,Susceptibility./Density,'b-o');
% grid on;
% xlabel('Singles Density')
% ylabel('Susceptibility/Density')