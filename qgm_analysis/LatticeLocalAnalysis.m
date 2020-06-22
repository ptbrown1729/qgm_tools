% BinEndPts = sqrt(linspace(0,40^2,40));
BinEndPts = sqrt(linspace(0,50^2,40));
% Frqs = [320:3:359,360:1.5:420];
Frqs = [300:3:410];

% Folder2 = {2017 04 04 17,1,1};
Folder2 = {2017 04 04 22 1 1};
df2 = DataFolder(Folder2,[],BinEndPts);

DistGrid = df2.DistGrid;
BinAvg = df2.BinAvg;

% Folder3 = {2017 04 04 18,1,1};
Folder3 = {2017 04 04 23 1 1};
df3 = DataFolder(Folder3,[],BinEndPts,'external',DistGrid);

% Folders4 = {2017 04 04 19,1,1};
Folders4 = {2017 04 04 24 1 1};
df4 = DataFolder(Folder3,[],BinEndPts,'external',DistGrid);

AzAvg_Stack = 1/3*(df2.Occs_AzAvgStack+df3.Occs_AzAvgStack+df4.Occs_AzAvgStack);
AzAvg_StackUnc = 1/3*sqrt((df2.Occs_AzAvgStackUnc).^2+(df3.Occs_AzAvgStackUnc).^2+(df4.Occs_AzAvgStackUnc).^2);

%%
%display

BinIndex = 1;
figure;
errorbar(Frqs,AzAvg_Stack(BinIndex,:),AzAvg_StackUnc(BinIndex,:));
hold on;
grid on;
title(sprintf('Bin Index %d, Mean Distance = %0.2f',BinIndex,BinAvg(BinIndex)));
% grid on;
% plot(Frqs,df2.Occs_AzAvgStack(BinIndex,:));
% plot(Frqs,df3.Occs_AzAvgStack(BinIndex,:));

%%
%fit each bin
FurthestBinToFit = 25;
ShowFits = 1;

InterpPts = linspace(min(Frqs),max(Frqs),300);
Centers = [];
Widths = [];

InitP = [[347,8,-0.8,1],[375,4,-0.5,0],[385,4,-0.8,0]];
for ii = 1:FurthestBinToFit
%     FirstDeriv = diff(smooth(AzAvg_Stack(ii,:)));
%     FirstDerivPts = 0.5*(Frqs(2:end)+Frqs(1:end-1));
% 
%     
%     SecondDeriv = diff(smooth(FirstDeriv));
%     SecondDerivPts = 0.5*(FirstDerivPts(2:end)+FirstDerivPts(1:end-1));
    
    
    
    FixedP = [[0,0,0,0],[0,0,0,1],[0,0,0,1]];
    LowerLims = [[min(Frqs),0,-1.1,-inf],[min(Frqs),0,-1.1,-inf],[min(Frqs),0,-1.1,-inf]];
    UpperLims = [[max(Frqs),max(Frqs),0,inf],[max(Frqs),max(Frqs),0,inf],[max(Frqs),max(Frqs),0,inf]];
    [Fp,Pn,Fh,SE] = fit1D(Frqs,AzAvg_Stack(ii,:),[],{'lorentzian1D','lorentzian1D','lorentzian1D'},InitP,FixedP,LowerLims,UpperLims);
    Centers = cat(1,Centers,Fp([1,5,9]));
    Widths = cat(1,Widths,Fp([2,6,10]));
    if ShowFits
        figure('name',sprintf('Bin %d',ii));
        errorbar(Frqs,AzAvg_Stack(ii,:),AzAvg_StackUnc(ii,:));
        hold on;
        plot(InterpPts,Fh(InterpPts));
        grid on;
    end
    
    InitP([1,2,5,6,9,10]) = Fp([1,2,5,6,9,10]);
    
end

figure('name',sprintf('Peak fitting,up to bin %d',FurthestBinToFit));
plot(BinAvg(1:size(Centers,1)),Centers);
ylabel('Peak Frq (KHz)');
xlabel('Distance (Lattice Sites)');
grid on;


%%
%Fit lattice depths
RetroAtten = 0.47;
LattServos = [0,10];
InitParams = [4,RetroAtten];
FixedParams = [0,1];

LattDepths = [];

for jj = 1:size(Centers,1)
    
    Peaks1 = [0,Centers(jj,1)]; Peaks2 = [0,Centers(jj,2)]; Peaks3 = [0,Centers(jj,3)];
    Width1 = [1,Widths(jj,1)]; Width2 = [1,Widths(jj,2)]; Width3 = [1,Widths(jj,3)];
%     FitParams = Lattice_Interpolated(LattServos,Peaks1,Peaks2,Peaks3,Width1,Width2,Width3,InitParams,FixedParams);
    FitParams = Fit_Lattice_Depth(LattServos,Peaks1,Peaks2,Peaks3,Width1,Width2,Width3,InitParams,FixedParams,[1,0,1]);
    LattDepths = cat(1,LattDepths,LattServos(2)*FitParams(1));
end

figure('name','Lattice Depth Vs. Position')
plot(BinAvg(1:length(LattDepths)),LattDepths);
ylabel('Lattice Depth (Er)');
xlabel('Distance (Lattice Sites)')
grid on;
title(sprintf('Lattice Depth Vs. Radius, Efield Attenuation = %0.2f',RetroAtten));

%%
%Fit potential...
fn = @(P,R) P(1) - 0.5*P(2).^2*R.^2;
fitfn = @(P) fn(P,BinAvg(1:length(LattDepths)))-LattDepths;

Vo = 58;
Omega = 0.2;

TrapFitParams = lsqnonlin(fitfn,[58,0.2]);
mLi = 9.9883e-27;
ErecHz = 14.66e3;
LattSpace = 1064e-9/sqrt(2);
hbar = 1.054e-34;
OmegaRealHz = TrapFitParams(2)/sqrt(mLi/abs(ErecHz*2*pi*hbar))/LattSpace;

Rinterp = linspace(min(BinAvg(1:length(LattDepths))),max(BinAvg(1:length(LattDepths))),100);
figure('name','Lattice Depth Vs. Position')
plot(BinAvg(1:length(LattDepths)),LattDepths,'r.');
hold on;
plot(Rinterp,fn(TrapFitParams,Rinterp),'b');
ylabel('Lattice Depth (Er)');
xlabel('Distance (Lattice Sites)')
grid on;
title(sprintf('Lattice Depth Vs. Radius, Efield Attenuation = %0.2f \n Vo = %0.1f Er, Omega = (2pi) %0.1fHz',RetroAtten,TrapFitParams(1),OmegaRealHz/(2*pi)));


%%
%U/ts 

%load data
FilePath = 'LatticeBandData_Interpolation.txt';
NBands = 8;

%Read data.
Data = dlmread(FilePath,',',1,0);
Depths = Data(:,1);
EFieldRetroAttenuation = Data(:,2);

NDepths = length(unique(Depths));
DepthsGrid = reshape(Depths,[length(Depths)/NDepths,NDepths]);
RetroAttenGrid = reshape(EFieldRetroAttenuation,[length(Depths)/NDepths,NDepths]);
DataGrid = reshape(Data,[length(Depths)/NDepths,NDepths,size(Data,2)]);

%interpolate
txFn = @(Depth,Atten) interp2(DepthsGrid,RetroAttenGrid,DataGrid(:,:,7),Depth,Atten);
tyFn = @(Depth,Atten) interp2(DepthsGrid,RetroAttenGrid,DataGrid(:,:,8),Depth,Atten);
tdiagFn = @(Depth,Atten) interp2(DepthsGrid,RetroAttenGrid,DataGrid(:,:,9),Depth,Atten);
uFn = @(Depth,Atten) interp2(DepthsGrid,RetroAttenGrid,DataGrid(:,:,10),Depth,Atten);

UoverTfn = @(Depth,Atten) uFn(Depth,Atten)./(0.5*(txFn(Depth,Atten)+tyFn(Depth,Atten)));

figure('name','Lattice Depth Vs. Position')
plot(BinAvg(BinAvg<31),UoverTfn(LattDepths(BinAvg<31)*1.12/10,RetroAtten),'r.');
hold on;
plot(Rinterp,UoverTfn(fn(TrapFitParams,Rinterp)*1.12/10,RetroAtten),'b');
ylabel('U/t');
xlabel('Distance (Lattice Sites)')
grid on;
title(sprintf('U over t Vs. Radius, Efield Attenuation = %0.2f, Inferred for ServoV = 1.12V',RetroAtten));