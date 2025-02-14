function [] = NRP_Close(nrp)
fclose(nrp)
delete(nrp)
clear nrp
end