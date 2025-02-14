function [ ] = N6705A_Close( N6705A )
%N6705A_Close Turns off and disconnects from N6705A Power Supply
%   No detailed explanation.

%fprintf(N6705A,'OUTP OFF,(@1:4)');
fclose(N6705A)
delete(N6705A)
clear N6705A
end