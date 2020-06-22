function [EntropyHarmonicTrap] = Entropy2DHarmonicTrap(T,mu,omega)

warning('this function is not the correct format (Entropy2DHarmonicTrap)');

kB=1.38e-23;
hbar = 1.055e-34;
EntropyHarmonicTrap = kB*(...
    3*(kB*T/(hbar*omega)^2)*polylog(3,-exp(mu/T))-...
    mu*kB*T/(hbar*omega)^2*polylog(2,-exp(mu/T)));

end


