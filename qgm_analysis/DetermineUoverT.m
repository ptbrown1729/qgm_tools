%Field from RF spec

Field = 616.06;
FieldUnc = 0.12;

PeakShift = 75.8697-75.86525; %MHz

U = @(Field) PeakShift*a3D(Field,2)/(a3D(Field,2)-a3D(Field,3));

UBest = U(Field);
UMax = U(Field+FieldUnc);
UMin = U(Field-FieldUnc);

UErr = mean([UMax-UBest,UBest-UMin]);

%Band structure calc
LattDepth10V = 59;
RetroEFieldAtten = 0.48;
LattDepth = 1.5;
[t_Er,tx_Er,ty_Er,t_diag,U_BandStruct,UOverT_BandStruct] =  Lattice_2D_Asymmetric_Interpolated(LattDepth10V/10*LattDepth,RetroEFieldAtten);

abohr = 5.2918e-11;
Er = 14.66e3;
U_BandStruct = U_BandStruct*a3D(Field,2)/(-888*abohr)*Er;
UOverT_BandStruct = UOverT_BandStruct*a3D(Field,2)/(-888*abohr);

fprintf('t = %0.1f Hz \n',t_Er*Er);
fprintf('Band Struct U = %0.1f KHz \n',U_BandStruct/1e3);
fprintf('Band Struct U/t = %0.2f \n',UOverT_BandStruct);
fprintf('Spectroscopy U = %0.2f +/- %0.2f KHz \n',UBest*1e3,UErr*1e3);
fprintf('USpec/t = %0.2f \n\n',UBest*1e6/(t_Er*14.66e3));