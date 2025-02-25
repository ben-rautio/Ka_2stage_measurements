%% load other WS, calculate everything from measured data

Attenuation_scaled = repmat(output_off_avg - 0.17, length(Pdes), 1); % Output Attenuation
InputOffsets_scaled = repmat(inp_off_avg + 0.1, length(Pdes), 1); % Input Offsets

Pout1 = FundTone_Out + Attenuation_scaled; % Use whichever measured power you want. Or all of them.
Pout2 = FundTone_OutFreq + Attenuation_scaled; % Use whichever measured power you want. Or all of them.
Pout3 = FundTone_OutFreq2 + Attenuation_scaled; % Use whichever measured power you want. Or all of them.

Pout_watts1 = 10.^(Pout1./10)./1000;
Pout_watts2 = 10.^(Pout2./10)./1000;
Pout_watts3 = 10.^(Pout3./10)./1000;

Pin1 = FundTone_In + InputOffsets_scaled;
Pin2 = FundTone_InFreq + InputOffsets_scaled;
Pin3 = FundTone_InFreq2 + InputOffsets_scaled;

Pin_watts1 = 10.^(Pin1./10)./1000;
Pin_watts2 = 10.^(Pin2./10)./1000;
Pin_watts3 = 10.^(Pin3./10)./1000;

Pdc_Drain1 = Vd1_eng.*Id1_eng;
Pdc_Gate1 = Vg1_eng.*Ig1_eng;
Pdc_tot1 = Pdc_Drain1 + Pdc_Gate1;

Pdc_Drain2 = Vd2_eng.*Id2_eng;
Pdc_Gate2 = Vg2_eng.*Ig2_eng;
Pdc_tot2 = Pdc_Drain2 + Pdc_Gate2;

Pdc_tot = Pdc_tot1 + Pdc_tot2;

DE1 = (Pout_watts1./Pdc_tot).*100;
PAE1 = ((Pout_watts1-Pin_watts1)./Pdc_tot).*100;
DE2 = (Pout_watts2./Pdc_tot).*100;
PAE2 = ((Pout_watts2-Pin_watts2)./Pdc_tot).*100;
DE3 = (Pout_watts3./Pdc_tot).*100;
PAE3 = ((Pout_watts3-Pin_watts3)./Pdc_tot).*100;

Gain1 = Pout1 - Pin1;
Gain2 = Pout2 - Pin2;
Gain3 = Pout3 - Pin3;

%% plot
close all

pdel_ADS = DriveupADSdataNG902stage.pdel(DriveupADSdataNG902stage.freq == 40);
PAE_ADS = DriveupADSdataNG902stage.pae(DriveupADSdataNG902stage.freq == 40);
Gain_ADS = DriveupADSdataNG902stage.gain(DriveupADSdataNG902stage.freq == 40);

figure
hold on; box on; grid on; grid minor; set(gcf,'color','w')
yyaxis left
ylim([0 20])
plot(pdel_ADS,Gain_ADS)
plot(Pout1(:,11),Gain1(:,11))
plot(Pout2(:,11),Gain2(:,11),'Color','r')
plot(Pout3(:,11),Gain3(:,11),'Color','g')
ylabel('Gain (dB)')
% ylim([4 13])
xlim([20 30])
yyaxis right
ylim([10 30])
plot(pdel_ADS,PAE_ADS)
plot(Pout1(:,11),PAE1(:,11))
plot(Pout2(:,11),PAE2(:,11),'Color','r')
plot(Pout3(:,11),PAE3(:,11),'Color','g')
ylabel('PAE (%)')
xlabel('Pout (dBm)')
legend("Sim Gain", "Meas Gain1", "Meas Gain2","Meas Gain3","Sim PAE", "Meas PAE1", "Meas PAE2","Meas PAE3")
title("CW Measurements vs. Simulated")
hold off