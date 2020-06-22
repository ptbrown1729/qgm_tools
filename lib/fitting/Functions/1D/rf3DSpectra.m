function [V,PNames,ArgStr] = rf3DSpectra(P,X)
%P = [AMol,Eb,Eb',Sign,Bg,FreePeak,HWHM,AAtom]. Sign parameter decides which way peak faces.
PNames = {'AMol','Eb','EbPrime','Sign','Bg','FreePeak','HWHM','AAtom'};
ArgStr = '';

if isempty(P)||isempty(X)
    V = 0;
else
    AMol = P(1); Eb = P(2); EbPrime = P(3); Sign = sign(P(4)); Bg = P(5);FreePeak = P(6); HWHM = P(7); AAtom = P(8);
    V = Bg+AMol*X.^-2.*sqrt(abs((X-FreePeak)*Sign-Eb)).*((X-FreePeak)*Sign-Eb+EbPrime).^-1+lorentzian1D([FreePeak,HWHM,AAtom,0],X,0);
    %V = Bg+AMol*X.^-2.*sqrt(abs((X)*Sign-Eb)).*((X)*Sign-Eb+EbPrime).^-1+lorentzian1D([FreePeak,HWHM,AAtom,0],X,0);
    V(isnan(V)) = 0;
    V(isinf(V)) = 0;
end
end
