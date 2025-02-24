pdel_ADS = DriveupADSdataNG902stage.pdel(DriveupADSdataNG902stage.freq == 40);
PAE_ADS = DriveupADSdataNG902stage.pae(DriveupADSdataNG902stage.freq == 40);
Gain_ADS = DriveupADSdataNG902stage.gain(DriveupADSdataNG902stage.freq == 40);

figure
hold on; box on; grid on; grid minor; set(gcf,'color','w')
yyaxis left
ylim([0 20])
plot(pdel_ADS,Gain_ADS)
ylabel('Gain (dB)')
% ylim([4 13])
xlim([20 30])
yyaxis right
ylim([10 30])
plot(pdel_ADS,PAE_ADS)
ylabel('PAE (%)')
xlabel('Pout (dBm)')
% legend('38 GHz','38.2 GHz','38.4 GHz','38.6 GHz', ...
%     '38.8 GHz','39 GHz','39.2 GHz','39.4 GHz','39.6 GHz','39.8 GHz','40 GHz'); % or however you're doing it% ylim([0 40])
hold off