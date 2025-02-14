function [ E36232A ] = E36232A_Setup_USB(Top)

if Top == 1
    E36232A = visa('keysight',"USB0::0x2A8D::0x3002::MY59001597::0::INSTR")
else
    E36232A = visa('keysight',"USB0::0x2A8D::0x3002::MY59001580::0::INSTR")
end
fopen(E36232A)

end