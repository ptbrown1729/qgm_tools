function temperature=getTemp(file_name, t_add, Constants)
%Takes the name of a log file and an amount of time that should be added to
%the get the appropriate time of flight to compute temperature. Returns the
%temperature.

%Get constants
K = Constants.K;

%Run settings file.
[pathstr, ~, ~] = fileparts(file_name);
settings_path = pathstr(1:end-9);
run(fullfile(settings_path, 'settings'))

%Extract values from log file
[vals, keywords] = parselog(file_name);
%find proper values. Note that if names are changed in readimg or in Cicero
%these will have to be changed.
tindex = find(strcmp(keywords, 'TimeOfFlight'));
sxindex = find(strcmp(keywords, 'Sx'));
syindex = find(strcmp(keywords, 'Sy'));

tof = str2double(vals(:, tindex)) + t_add;
sx = str2double(vals(:, sxindex));
sy = str2double(vals(:, syindex));

%fit using square root function
params = polyfit(tof.^2, sy.^2, 1);

%plot data and fit
plot(tof, sy, '.');
hold on;
%plot(T,sy,'.');
%hold on;
temps = linspace(tof(1), tof(length(tof)), 25);
line = sqrt(params(1) * temps .^ 2 + params(2));
plot(temps, line);

%extract temperature
temperature=mLi * (1000 * pixel_size / magnification)^2 * params(1) / K;