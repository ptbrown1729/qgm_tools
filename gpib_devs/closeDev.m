function closeDev(GPIB_Dev)

    clrdevice(GPIB_Dev) %clear hardware buffer of device
    fclose(GPIB_Dev); %close device
    delete(GPIB_Dev); %delete handle object
    clear GPIB_Dev; %remove the variable
end