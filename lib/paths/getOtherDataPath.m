function [Path] = getOtherDataPath(Date)
    RootOtherDirectory = fullfile('//','128.112.86.75','lithium','OTHER DATA');
    %DateStr1 = datestr(Date,'yyyy/mm');
    %DateStr2 = datestr(Date,'-mmm/dd');
    %DateStr = strcat(DateStr1,DateStr2);
    DateStr = fullfile(datestr(Date,'yyyy'),sprintf('%s-%s',datestr(Date,'mm'),datestr(Date,'mmm')),datestr(Date,'dd'));
    Path = fullfile(RootOtherDirectory,DateStr);
end