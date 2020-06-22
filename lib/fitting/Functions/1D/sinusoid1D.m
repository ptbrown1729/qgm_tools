function [V, PNames, ArgStr] = sinusoid1D(P, X)
% sinusoid1D(P, X), functional form Amp * sin(2 * pi * F * X + Phase) + Bg
%
% P: vector of paramters, P = [F,A,Phase,Bg]
%  
% X: array of locations of arbitrary size
%

PNames = {'F','Amp','Phase','Bg'};
ArgStr = 'Amp*sin(2*pi*F*X+Phase)+Bg';

if isempty(P) || isempty(X)
    V = 0;
else
    F = P(1); 
    Amp = P(2); 
    Phase = P(3); 
    Bg = P(4);
    
    V = Amp * sin(2 * pi * F * X + Phase) + Bg;
end
end



