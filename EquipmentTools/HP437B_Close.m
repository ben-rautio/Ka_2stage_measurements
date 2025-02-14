function [  ] = HP437B_Close(HP437B)
%HP437B_SETUP Summary of this function goes here
%   Detailed explanation goes here
fclose(HP437B)
delete(HP437B)
clear HP437B
end

