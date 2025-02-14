function [ ] = SMJ100A_OnOff( SMJ100A, ChangeTo )
%SMJ100A_OnOff Turns the SMJ100A On and Off
%   Detailed explanation goes here



if strcmp(ChangeTo, 'On') | strcmp(ChangeTo, 'on')
   [status, Result] = rs_send_query( SMJ100A, 'OUTP:STAT ON; *OPC?' );
   if( ~status ); clear;  return; end
elseif strcmp(ChangeTo, 'Off') | strcmp(ChangeTo, 'off')
   [status, Result] = rs_send_query( SMJ100A, 'OUTP:STAT OFF; *OPC?' );
if( ~status ); clear;  return; end
else
    disp('Unknown ChangeTo for SMJ100A_OnOff')
end
end

