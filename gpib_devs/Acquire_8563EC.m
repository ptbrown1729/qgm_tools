if ~isempty(instrfind)
    fclose(g);
    delete(g)
    clear g
end

%open connection
BoardIndex = 0;
PrimaryAddress = 21;
g = gpib('ni',BoardIndex,PrimaryAddress);
g.InputBufferSize=50000;
g.Timeout=120;
fopen(g);
g.EOSCharCode = 'LF';


%read data
ID = readDat(g,'ID?');
FirmWareRev = readDat(g,'REV?');
%measurement type, units, etc.
RefLvl = readDat(g,'RL?',1); %dBm
ResolutionBW = readDat(g,'RB?',1);
VideoBW = readDat(g,'VB?',1);
AverageOn = readDat(g,'VAVG?');

%marker info
MarkerAmp = readDat(g,'MKA?',1);
MarkerFrq = readDat(g,'MKF?',1);
Marker3dBBW = readDat(g,'MKBW?',1);

%frequency info
Span = readDat(g,'SP?',1); %Hz
StartFrq = readDat(g,'FA?',1);
CtrFrq = readDat(g,'CF?',1);
StepSize = readDat(g,'SS?',1);
FrqData = StartFrq:(Span/600):(StartFrq+Span);
FrqUnits = 'Hz';

%Spectrum
SpecData = readDat(g,'TRA?',1);
SpecUnits = 'dBm';

% stopasync(g)
closeDev(g);

Fname = sprintf('8563EC_SpecAnalyzer_Data_%s.mat',datestr(now,'yyyy_mm_dd_hh;MM;ss'));
save(Fname,'ID','FirmWareRev','RefLvl','ResolutionBW','VideoBW','AverageOn','MarkerAmp'...
    ,'MarkerFrq','Marker3dBBW','Span','StartFrq','CtrFrq','StepSize','FrqData','FrqUnits','SpecData','SpecUnits');


