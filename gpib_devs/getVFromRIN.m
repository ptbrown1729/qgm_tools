function val = getVFromRIN(RIN, Vdc, BW)
%assume input unit is dBc/sqrt(Hz)

%output units = Vrms
%dBc = 10*log10((Vrms/Vdc)^2/Bw)

val = sqrt(10^(RIN / 10) * BW) * Vdc;

end