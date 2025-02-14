function [ ] = FSQ26_Close( FSQ26 )
%FSQ26_Clear Disconnects from the FSQ26 Vector Signal Analyzer
%   No detailed explanation

fclose(FSQ26); 
delete(FSQ26)
clear FSQ26
end

