function [power] = HP437B_ReadPowerFreq(HP437B, freq)

% HP437B Cal Factor Lookup Table
CalFac = zeros(19,2);
CalFac(:,1) = [0.1e9,0.5e9,2e9,3e9,4e9,5e9,6e9,7e9,8e9,9e9,10e9,11e9,12e9,13e9,14e9,15e9,16e9,17e9,18e9];
CalFac(:,2) = [100,99.5,99.4,99.1,99.1,98.3,98.6,98.5,98.5,98.6,98.1,98.7,96.7,98.5,97.1,97.6,97.7,96.7,96.5];

% Example syntax
% sprintf('KB%.1f%%',99.5)

check = find(CalFac == freq); % Check if the desired frequency has an associated Cal Factor
x = size(check);

if x(1) == 0 % If it doesn't, do a normal power fetch
    power = str2num(query(HP437B, 'GET'));
else % Otherwise, input the desired calibration factor
    percentage = CalFac(check, 2);
    string = sprintf('KB%.1f%%',percentage);
    power = str2num(query(HP437B, string));
end

% str2num(query(HP437B, 'KB99.5%'));
% power = str2num(query(HP437B, 'FR2GZ'));

end

