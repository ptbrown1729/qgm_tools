function vpp = dBm2vpp(dBm)
    vpp = sqrt(10^(dBm/10)*1e-3*50)*2*sqrt(2);
end