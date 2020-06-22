function [temp] = getTempSingle(sigma, time_of_flight, pixel_area, Constants)
% estimate temperature from a single expansion time
%TOF in ms.
mLi = Constants.mass;
kb = Constants.K;

time_of_flight = time_of_flight * 1e-3; 
temp = mLi * pixel_area * (sigma / time_of_flight)^2 / kb;

end

