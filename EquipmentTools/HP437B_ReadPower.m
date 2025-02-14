function [power] = HP437B_ReadPower(HP437B)
%HP437B_READPOWER Summary of this function goes here
%   Detailed explanation goes here
power = str2num(query(HP437B, 'GET'));
end

