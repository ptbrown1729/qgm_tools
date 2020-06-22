if ~isempty(instrfind)
    fclose(g);
    delete(g)
    clear g
end

%open connection
BoardIndex = 0;
PrimaryAddress = 11;
g = gpib('ni',BoardIndex,PrimaryAddress);
g.InputBufferSize=50000;
g.Timeout=120;
fopen(g);
g.EOSCharCode = 'LF';


%read data
ID = readDat(g,'*IDN?');

 %measurement type, units, etc.
MeasIndex = readDat(g,'MEAS?-1',1);
AllowedMeasTypes = {'Spectrum','PSD','Time Record','Octave Analysis'};
MeasType = AllowedMeasTypes{MeasIndex+1};

DispIndex = readDat(g,'DISP?-1',1);
AllowedDispTypes = {'Log Magnitude','Linear Magnitude','Real Part','Imaginary Part','Phase'};
DispType = AllowedDispTypes{DispIndex+1};

MeasUnitIndex = readDat(g,'UNIT?-1',1);
AllowedUnits = {'Volts Pk','Volts RMS','dBV','dBVrms'}; %01,2,3,4
MeasUnits = AllowedUnits{MeasUnitIndex+1};

UnitTypeIndex = readDat(g,'VOEU?-1',1);
AllowedUnitTypes = {'Volts','EUs'};
UnitType = AllowedUnitTypes{UnitTypeIndex+1};

WindowFnIndex = readDat(g,'WNDO?-1',1);
AllowedWindowFns = {'Uniform','Flattop','Hanning','Blackman-Harris'}; %0,1,2,3
WindowFn = AllowedWindowFns{WindowFnIndex+1};

% AnalysisMd = readDat(g,'ANAM?-1',1);
%read marker position
MarkerPos = readDat(g,'MRKX?-1',1);
MarkerVal = readDat(g,'MRKY?-1',1);

% input range information
arange_index = readDat(g, 'ARNG?', 1);
AllowedRanges = {'Manual', 'Auto'};
range_mode = AllowedRanges{arange_index+1};

input_range = readDat(g, 'IRNG?', 1);

%information about averaging 
AverageOn = readDat(g,'AVGO?',1);
NumAvgs = readDat(g,'NAVG?',1);

AvgTypeIndex = readDat(g,'AVGT?',1);
AllowedAvgTypes = {'RMS','Vector','Peak Hold'};
AvgType = AllowedAvgTypes{AvgTypeIndex+1};

Avg_ModeIndex = readDat(g,'AVGM?',1);
AllowedAvgModes = {'Linear','Exponential'};
AvgMode = AllowedAvgModes{Avg_ModeIndex+1};

AvgOverlap = readDat(g,'OVLP?',1); %percentage
%frequency information
SpanIndex = readDat(g,'SPAN?',1);
AllowedSpans = [191e-3,382e-3,763e-3,...
    1.5,3.1,6.1,12.2,24.4,48.75,97.5,...
    195,390,780,1.56e3,3.125e3,6.25e3,...
    12.5e3,25e3,50e3,100e3]; %corresponding to numbers 0 to 19
Span = AllowedSpans(SpanIndex+1);
NumBins = 400; %#0-#399
LineWidth = Span/NumBins;
StartFrq = readDat(g,'STRF?',1);
CtrFrq = readDat(g','CTRF?',1);
FrqData = StartFrq:LineWidth:(StartFrq+Span-LineWidth);

%read spectrum
SpecData = readDat(g,'SPEC?-1',1);

% stopasync(g)
closeDev(g);

fname = sprintf('%s_SR770_data.mat', datestr(now, 'yyyy_mm_dd_hh;MM;ss'));
save(fname,'ID','MeasType','DispType','MeasUnits','UnitType','WindowFn',...
    'range_mode', 'input_range', 'MarkerPos','MarkerVal','AverageOn','NumAvgs','AvgType','AvgMode',...
    'AvgOverlap','Span','NumBins','LineWidth','StartFrq','CtrFrq','FrqData','SpecData');


