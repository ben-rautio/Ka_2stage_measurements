function [ ] = E4419B_Close( E4419B )
%E4419B_Close Disconnects from E4419B Power Meter. 
%   No detailed explanation
fclose(E4419B)
delete(E4419B)
clear E4419B
end

