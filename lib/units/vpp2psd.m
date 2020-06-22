function psd = vpp2psd(vpp,bandwidth)
    %psd in Vrms/sqrt(Hz)
    
    psd = (vpp/2/sqrt(2))/sqrt(bandwidth);
    psd = psd/sqrt(2); %don't understand why there should be this factor of two here, but it is used by the Sr770 FFT Network Analyzer.
end