function [nrp] = NRP_50SN_setup()
% nrp = visa('keysight',"USB::0x0aad::0x0145::102483::INSTR") % original
% nrp = visa('keysight',"USB0::0x0AAD::0x0138::102270::0::INSTR")
nrp = visa('keysight',"USB0::0x0AAD::0x0162::101101::0::INSTR")

fopen(nrp)
fprintf(nrp, 'INIT:CONT OFF');
pause(9e-1)
fprintf(nrp, 'SENS:FUNC "POW:AVG"');
pause(9e-1)
fprintf(nrp, 'SENS:AVER:COUN:AUTO OFF');
pause(9e-1)
fprintf(nrp, 'SENS:AVER:COUN 16');
pause(9e-1)
fprintf(nrp, 'SENS:AVER:STAT ON');
pause(9e-1)
fprintf(nrp, 'SENS:AVER:TCON REP');
pause(9e-1)
fprintf(nrp, 'FORMAT ASCII');
end