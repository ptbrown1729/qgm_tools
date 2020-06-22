function [V,PNames,ArgStr] = rf2DSpectra(P,X)
%P = [A,Eb,Eb',Sign,Bg,FreePeak]. Sign parameter decides which way peak faces.
%RFSpectra for 2D. For reference see PRL 108, 045302.
PNames = {'FreePeakCenter','Amp','Eb','EbPrime','Sign','Bg'};
ArgStr = 'Amp*heaviside((X-XcFree)*Sign-Eb)./(X-XcFree).^2./(log((abs((X-XcFree)*Sign-Eb)./EbP)).^2+pi^2)+Bg';
%not sure how will handle pi...might need to replace with numbers in
%string...

if isempty(P)||isempty(X)
    V = 0;
else
    XcFree = P(1); Amp = P(2); Eb = P(3); EbP = P(4); Sign = sign(P(5)); Bg = P(6);
    %V = A*heaviside((X-XcFree)*Sign-Eb)./(X-XcFree).^2./(log((abs((X-XcFree)*Sign-Eb)./EbP)).^2+pi^2) + Bg;%.*(log(abs(Eb./EbP)).^2./(log((abs((X-XcFree)*Sign-Eb)./EbP)).^2+pi^2)) + Bg;
    V = Amp*heaviside((X-XcFree)*Sign-Eb)./(X-XcFree).^2./(log((abs((X-XcFree)*Sign-Eb)./EbP)).^2+pi^2)+Bg;
    V(isnan(V)) = 0;
end
end
