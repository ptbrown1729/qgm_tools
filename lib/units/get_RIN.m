function RIN = get_RIN(val,unit,BW_Hz,Vdc_Volts)
%Possible units for val = 'dBV', 'dBm', 'Vrms'
%Band Width must be in Hz
%DC voltage must be in Volts

if strcmp(unit,'dBV')
    %this is not standard usage for dBV, but is how the SR7700 FFT Network
    %analyzer defines dBV...it defines dBv = 20*log10(V/1volt). And
    %analagously for dBVrms
    Vrms = 10.^(val/20); %This is in Volts RMS
    RIN = 10*log10((Vrms/Vdc_Volts).^2/BW_Hz);
elseif strcmp(unit,'dBm')
    Vrms = sqrt(50e-3*10.^(val/10));
    RIN = 10*log10((Vrms/Vdc_Volts).^2/BW_Hz);
elseif strcmp(unit,'Vrms')
    Vrms = val;
    %in dB/Hz
    RIN = 10*log10((Vrms/Vdc_Volts).^2/BW_Hz);
elseif strcmp(unit,'Vrms/sqrt(Hz)')
    RIN = 20*log10(val/Vdc_Volts);
else
    Vrms = 1;
    RIN = 0;
end
