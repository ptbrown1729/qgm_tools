function [EntropyHarmonicTrap] = EntropyHarmonicTrap(T,mu,omega)

EntropyHarmonicTrap = kB*(...
    3*(kB*T/(hbar*omega)^2)*polylog(3,-exp(mu/T))-...
    mu*kB*T/(hbar*omega)^2*polylog(2,-exp(mu/T)));

end