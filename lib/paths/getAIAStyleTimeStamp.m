function [ AIAStyleDateString ] = getAIAStyleTimeStamp( DateString )
%DATEFORMATTOAIA Summary of this function goes here
%   Detailed explanation goes here

AIAStyleDateString = datestr(datenum(DateString),'yyyy-mm-dd-HH;MM;SS');

end

