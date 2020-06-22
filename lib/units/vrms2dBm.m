function dbm = vrms2dBm(vrms)
    dbm = 10*log10((vrms)^2/50/1e-3);
end