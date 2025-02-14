function [] = E36232A_Close_USB(E36232A)
fclose(E36232A)
delete(E36232A)
clear E36232A
end