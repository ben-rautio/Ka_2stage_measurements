function [ ] = NRP67T_Close( NRP67T )
% NRP67T_Close Turns off and disconnects from NRP67T Power Supply
%   No detailed explanation.

fclose(NRP67T)
delete(NRP67T)
clear NRP67T
end