function [V, PNames] = thomasFermi2D(P,X,Y)
%P = [CxTF,CyTF,RxTF,RyTF,ATF,ThetaTF,Bg]
%2D Thomas-Fermi Function.
PNames = {'CxTF','CyTF','RxTF','RyTF','ATF','ThetaTF','Bg'};
ArgStr = '';

if isempty(P) || isempty(X) || isempty(Y)
    V = 0;
else
    CxTF = P(1);
    CyTF = P(2);
    RxTF = abs(P(3));
    RyTF = abs(P(4));
    ATF = abs(P(5));
    ThetaTF = P(6);
    Bg = P(7);
    TFPart = ATF * (1- ((X - CxTF) .* cos(ThetaTF) - (Y - CyTF) .* sin(ThetaTF)).^2 ./ RxTF.^2 -...
                       ((Y - CyTF) .* cos(ThetaTF) + (X - CxTF) .* sin(ThetaTF)).^2 ./ RyTF.^2);
    V= TFPart .* (TFPart > 0) + Bg;
    
end

end

