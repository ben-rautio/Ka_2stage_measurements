function [ ] = FSW_Close( FSW )
%N6705A_Close Turns off and disconnects from N6705A Power Supply
%   No detailed explanation.

%fprintf(N6705A,'OUTP OFF,(@1:4)');
fclose(FSW)
delete(FSW)
clear FSW
end