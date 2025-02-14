function [ ] = SMJ100A_Close( SMJ100A )
%SMJ100A_Close Disconnects from the SMJ100A Vector Signal Generator
%   No detailed explanation
[status, Result] = rs_send_query( SMJ100A, 'OUTP:STAT OFF; *OPC?' );
if( ~status ); clear;  return; end
fclose(SMJ100A); 
delete(SMJ100A)
clear SMJ100A
end

