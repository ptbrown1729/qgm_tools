function dbm = vpp2dBm(vpp)
    dbm = 10*log10((vpp/2/sqrt(2))^2/50/1e-3);
end