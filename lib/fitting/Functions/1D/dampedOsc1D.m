function [V,PNames,ArgStr] = dampedOsc1D(P,X)
%[V,PNames,ArgStr] = dampedOsc1D(P,X)
%This function is the solution to the problem d^2n/dt^2 + (2/DecayLen)*dn/dt + Omega_o^2*n = 0
%P = [InitPos,InitVel,Omega0,DecayLen,Cx,Bg]
%Omega_mod = sqrt(1/DecayLen^2 - Omega0^2)
%It seems to be more common to define decay len with a different factor of
%four in the above...i.e. so there is no factor of 2 in initial
%differential equation.
%InitPos is related to the initial velocity. It is the coefficient of the symmetric term, which is
% a cos or cosh term. Similarly, InitVel is the initial velocity. The
% coefficient of the asymmetric term is (InitPos/DecayLen +
% InitVel)/OmegaEff
PNames = {'InitPos','InitVel','Omega0','DecayLen','Cx','Bg'};
ArgStr = '';

if isempty(P)||isempty(X)
    V = zeros(size(X));
else
    InitPos = P(1); InitVel = P(2); Omega0 = P(3); DecayLen = P(4); Cx = P(5); Bg = P(6);
    
    EffFrq = sqrt(1/DecayLen^2-Omega0^2);
    %1/DecayT^2>Omega0^2 = Overdamped. Get cosh and sinh
    %1/DecayT^2<Omega0^2 = underamped. get cos and sin
    %1/DecayT^2=Omega0^2 = critically damped. Also handle this case.
    if DecayLen==0
        V = zeros(size(X));
    elseif 1/DecayLen^2~=Omega0^2
        %this term should be real modulo machine precesion
        SymmPart = 0.5*InitPos*real(exp(EffFrq*(X-Cx))+exp(-EffFrq*(X-Cx)));
        
        %this term is either purely real or purely imaginary. If it is
        %purely imaginary, we really want that term /1i. This is because
        %sin(x) and sinh(x) differ by a factor of 1/1i in their
        %definitions.
        AsymmPart = 0.5*(InitPos/DecayLen+InitVel)/EffFrq*(exp(EffFrq*(X-Cx))-exp(-EffFrq*(X-Cx)));
        if abs(imag(AsymmPart))>abs(real(AsymmPart))
            AsymmPart = AsymmPart/1i;
        end
        AsymmPart = real(AsymmPart);
        
        V = (SymmPart+AsymmPart).*exp(-(X-Cx)/DecayLen) + Bg;
    else
        %expression is different for the critically damped case.
        V = (InitPos + (InitPos/DecayLen+InitVel)*X).*exp(-(X-Cx)/DecayLen)+ Bg;
    end
end

%most common way I see Nans sneak in is InitPos = 0, exp(EffFrq*X) -> inf.
%Then 0*inf = nan. Typically these should be zeroes.
V(isnan(V)) = 0;

end



