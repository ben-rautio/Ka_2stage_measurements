% nrp = visa('keysight',"USB::0x0aad::0x0145::102545::INSTR")
% nrp = visa('keysight',"USB0::0x0AAD::0x0138::102270::0::INSTR") 
nrp = visa('keysight',"USB0::0x0AAD::0x0145::102483::0::INSTR") % the new one

fopen(nrp)

fprintf(nrp, 'CAL:ZERO:AUTO ONCE')

pause(5)

fclose(nrp)
delete(nrp)
clear nrp