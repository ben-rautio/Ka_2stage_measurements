%% load other WS, calculate everything from measured data
%cw datasets
FundTone_out_cw_s = load("driveup_38_40GHz.mat","FundTone_Out");
FundTone_out_cw = FundTone_out_cw_s.FundTone_Out;
FundTone_outFreq_cw_s = load("driveup_38_40GHz.mat","FundTone_OutFreq");
FundTone_outFreq_cw = FundTone_outFreq_cw_s.FundTone_OutFreq;
FundTone_outFreq2_cw_s = load("driveup_38_40GHz.mat","FundTone_OutFreq2");
FundTone_outFreq2_cw = FundTone_outFreq2_cw_s.FundTone_OutFreq2;
FundTone_In_cw_s = load("driveup_38_40GHz.mat","FundTone_In");
FundTone_In_cw = FundTone_In_cw_s.FundTone_In;
FundTone_InFreq_cw_s = load("driveup_38_40GHz.mat","FundTone_InFreq");
FundTone_InFreq_cw = FundTone_InFreq_cw_s.FundTone_InFreq;
FundTone_InFreq2_cw_s = load("driveup_38_40GHz.mat","FundTone_InFreq2");
FundTone_InFreq2_cw = FundTone_InFreq2_cw_s.FundTone_InFreq2;
Vg1_eng_cw_s = load("driveup_38_40GHz.mat","Vg1_eng");
Vg1_eng_cw = Vg1_eng_cw_s.Vg1_eng;
Vg2_eng_cw_s = load("driveup_38_40GHz.mat","Vg2_eng");
Vg2_eng_cw = Vg2_eng_cw_s.Vg2_eng;
Vd1_eng_cw_s = load("driveup_38_40GHz.mat","Vd1_eng");
Vd1_eng_cw = Vd1_eng_cw_s.Vd1_eng;
Vd2_eng_cw_s = load("driveup_38_40GHz.mat","Vd2_eng");
Vd2_eng_cw = Vd2_eng_cw_s.Vd2_eng;
Ig1_eng_cw_s = load("driveup_38_40GHz.mat","Ig1_eng");
Ig1_eng_cw = Ig1_eng_cw_s.Ig1_eng;
Ig2_eng_cw_s = load("driveup_38_40GHz.mat","Ig2_eng");
Ig2_eng_cw = Ig2_eng_cw_s.Ig2_eng;
Id1_eng_cw_s = load("driveup_38_40GHz.mat","Id1_eng");
Id1_eng_cw = Id1_eng_cw_s.Id1_eng;
Id2_eng_cw_s = load("driveup_38_40GHz.mat","Id2_eng");
Id2_eng_cw = Id2_eng_cw_s.Id2_eng;
%calibration offsets
avg_inp_off_s = load("MIWAVE_input_cal.mat","avg_inp_off");
avg_inp_off = avg_inp_off_s.avg_inp_off;
avg_output_atten_s = load("outputCal_data1_raw.mat","avg_output_atten");
avg_output_atten = avg_output_atten_s.avg_output_atten;
avg_GSG_atten_s = load("GSG_cal.mat","avg_gsg_loss");
avg_gsg_atten = avg_GSG_atten_s.avg_gsg_loss;
%calculate input offset to get input power
inp_offset = avg_inp_off - avg_gsg_atten; %add this to inp power meter to get inp power
%calculate output attenuation to get output power
output_atten = avg_output_atten + avg_gsg_atten; %add this to output power meter to get output power
%avg the last 5-10 samples to get a single value vs power
inp_off_avg = mean(inp_offset(size(inp_offset,2)-5:end,:),1);
output_off_avg = mean(output_atten(size(output_atten,2)-5:end,:),1);
%take into account approximate adapter losses
output_adapter_loss = 0.17;
input_adapter_loss = 0.1;
output_off_avg = output_off_avg - output_adapter_loss;
inp_off_avg = inp_off_avg + input_adapter_loss;
%load Pdes for replication in next section
Pdes_cw_s = load("driveup_38_40GHz.mat","Pdes");
Pdes_cw = Pdes_cw_s.Pdes;
%freq tested
freq_cw_s = load("driveup_38_40GHz.mat","Fo");
freq_cw = freq_cw_s.Fo;
%% calculate everything
Attenuation_scaled = repmat(output_off_avg, length(Pdes_cw), 1); % Output Attenuation
InputOffsets_scaled = repmat(inp_off_avg, length(Pdes_cw), 1); % Input Offsets

Pout1 = FundTone_out_cw + Attenuation_scaled; % Use whichever measured power you want. Or all of them.
Pout2 = FundTone_outFreq_cw + Attenuation_scaled; % Use whichever measured power you want. Or all of them.
Pout3 = FundTone_outFreq2_cw + Attenuation_scaled; % Use whichever measured power you want. Or all of them.

Pout_watts1 = 10.^(Pout1./10)./1000;
Pout_watts2 = 10.^(Pout2./10)./1000;
Pout_watts3 = 10.^(Pout3./10)./1000;

Pin1 = FundTone_In_cw + InputOffsets_scaled;
Pin2 = FundTone_InFreq_cw + InputOffsets_scaled;
Pin3 = FundTone_InFreq2_cw + InputOffsets_scaled;

Pin_watts1 = 10.^(Pin1./10)./1000;
Pin_watts2 = 10.^(Pin2./10)./1000;
Pin_watts3 = 10.^(Pin3./10)./1000;

Pdc_Drain1 = Vd1_eng_cw.*Id1_eng_cw;
Pdc_Gate1 = Vg1_eng_cw.*Ig1_eng_cw;
Pdc_tot1 = Pdc_Drain1 + Pdc_Gate1;

Pdc_Drain2 = Vd2_eng_cw.*Id2_eng_cw;
Pdc_Gate2 = Vg2_eng_cw.*Ig2_eng_cw;
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
tiledlayout("flow")
set(0,'DefaultLineLineWidth', 2)
DriveupADSdataNG902stage_s = load("Driveup_CW_ADS.mat","DriveupADSdataNG902stage");
DriveupADSdataNG902stage = DriveupADSdataNG902stage_s.DriveupADSdataNG902stage;
pdel_ADS = DriveupADSdataNG902stage.pdel(DriveupADSdataNG902stage.freq == 40);
PAE_ADS = DriveupADSdataNG902stage.pae(DriveupADSdataNG902stage.freq == 40);
Gain_ADS = DriveupADSdataNG902stage.gain(DriveupADSdataNG902stage.freq == 40);

nexttile
hold on; box on; grid on; grid minor; set(gcf,'color','w')
yyaxis left
ylim([0 20])
plot(pdel_ADS,Gain_ADS,'LineStyle','--','DisplayName','Sim Gain')
plot(Pout1(:,11),Gain1(:,11),'LineStyle','-','DisplayName','Meas. Gain')
% for i = 1:length(freq_cw)-1
%     plot(Pout1(:,i),Gain1(:,i),'Color','r','DisplayName',"Gain at "+num2str(freq_cw(i)/1e9)+" GHz")
% end
ylabel('Gain (dB)')
% ylim([4 13])
xlim([16 30])
yyaxis right
ylim([0 30])
yticks(0:3:30)
plot(pdel_ADS,PAE_ADS,'LineStyle','--','DisplayName','Sim PAE')
plot(Pout1(:,11),PAE1(:,11),'LineStyle','-','DisplayName','Meas. PAE')
% for i = 1:length(freq_cw)-1
%     plot(Pout1(:,i),PAE1(:,i),'Color','b','DisplayName',"PAE at "+num2str(freq_cw(i)/1e9)+" GHz")
% end
ylabel('PAE (%)')
xlabel('Pout (dBm)')
% legend("Sim Gain", "Meas Gain1", "Meas Gain2","Meas Gain3","Sim PAE", "Meas PAE1", "Meas PAE2","Meas PAE3")
legend
title("CW Measurements vs. Simulated")
hold off
%% S-parameters plotting
pth = "sparams_2stg.s2p";
pth1 = "NG90_2stage_sparams.s2p";

sp = sparameters(pth);
params_meas = sp.Parameters;
s11m = params_meas(1,1,:);
s21m = params_meas(2,1,:);
s22m = params_meas(2,2,:);

sp1 = sparameters(pth1);
params_sim = sp1.Parameters;
s11s = params_sim(1,1,:);
s21s = params_sim(2,1,:);
s22s = params_sim(2,2,:);

nexttile
hold on
plot(sp.Frequencies,20.*log10(abs(s11m(:))),'r','DisplayName','Meas. $S_{11}$','LineWidth',2)
plot(sp.Frequencies,20.*log10(abs(s21m(:))),'g','DisplayName','Meas. $S_{21}$','LineWidth',2)
plot(sp.Frequencies,20.*log10(abs(s22m(:))),'b','DisplayName','Meas. $S_{22}$','LineWidth',2)
plot(sp1.Frequencies,20.*log10(abs(s11s(:))),'--r','DisplayName','Sim $S_{11}$','LineWidth',2)
plot(sp1.Frequencies,20.*log10(abs(s21s(:))),'--g','DisplayName','Sim $S_{21}$','LineWidth',2)
plot(sp1.Frequencies,20.*log10(abs(s22s(:))),'--b','DisplayName','Sim $S_{22}$','LineWidth',2)
yticks(-25:5:20)
ylim([-25,20])
xticks(30e9:1e9:40e9);
xticklabels({'30' '31' '32' '33' '34' '35' '36' '37' '38' '39' '40'})
xlabel("Frequency (GHz)")
ylabel("|S-Parameters| (dB)")
grid on
h=legend;
set(h,'Interpreter','latex')
title("S-Parameters Meas. and Sim")
hold off

%% Modulated Measurements
%plot Pout, PAEavg, ACPR, EVM, Gain(?) over freq
%keep Pav ~constant at good power level
% can include direct dpd if it makes things look better
%get all the data from the modulated meas workspace
Pout1_mod_s = load("Modulated_Meas_WS.mat","Pout1");
Pout1_mod = Pout1_mod_s.Pout1;
Pout1_mod_dpd_s = load("Modulated_Meas_WS.mat","Pout1_dpd");
Pout1_mod_dpd = Pout1_mod_dpd_s.Pout1_dpd;
Gain1_mod_s = load("Modulated_Meas_WS.mat","Gain1");
Gain1_mod = Gain1_mod_s.Gain1;
Gain1_mod_dpd_s = load("Modulated_Meas_WS.mat","Gain1_dpd");
Gain1_mod_dpd = Gain1_mod_dpd_s.Gain1_dpd;
freq_mod_s = load("Modulated_Meas_WS.mat","Fo");
freq_mod = freq_mod_s.Fo;
PAE1_mod_s = load("Modulated_Meas_WS.mat","PAE1");
PAE1_mod = PAE1_mod_s.PAE1;
PAE1_mod_dpd_s = load("Modulated_Meas_WS.mat","PAE1_dpd");
PAE1_mod_dpd = PAE1_mod_dpd_s.PAE1_dpd;

ACPR_L_S = load("Modulated_Meas_WS.mat","adjpow_lower");
ACPR_L = ACPR_L_S.adjpow_lower;
ACPR_L_DPD_S = load("Modulated_Meas_WS.mat","adjpow_lower_dpd");
ACPR_L_DPD = ACPR_L_DPD_S.adjpow_lower_dpd;

ACPR_U_S = load("Modulated_Meas_WS.mat","adjpow_upper");
ACPR_U = ACPR_U_S.adjpow_upper;
ACPR_U_DPD_S = load("Modulated_Meas_WS.mat","adjpow_upper_dpd");
ACPR_U_DPD = ACPR_U_DPD_S.adjpow_upper_dpd;

EVM_s = load("Modulated_Meas_WS.mat","EVM");
EVM = EVM_s.EVM;
EVM_dpd_s = load("Modulated_Meas_WS.mat","EVM_dpd");
EVM_dpd = EVM_dpd_s.EVM_dpd;

figure
tiledlayout("flow")
nexttile
hold on
yyaxis left
for i = 1:1:numel(freq_mod)
    plot(Pout1_mod(:,i),Gain1_mod(:,i),'LineStyle','-','DisplayName',num2str(freq_mod(i)/1e9)+" GHz NDPD")
    plot(Pout1_mod_dpd(:,i),Gain1_mod_dpd(:,i),'LineStyle','--','DisplayName',num2str(freq_mod(i)/1e9)+" GHz DPD")
end
xlabel('Pout (dBm)')
ylabel('Gain (dB)')
yyaxis right
for i = 1:1:numel(freq_mod)
    plot(Pout1_mod(:,i),PAE1_mod(:,i),'LineStyle','-','DisplayName',num2str(freq_mod(i)/1e9)+" GHz NDPD")
    plot(Pout1_mod_dpd(:,i),PAE1_mod_dpd(:,i),'LineStyle','--','DisplayName',num2str(freq_mod(i)/1e9)+" GHz DPD")
end
ylabel('Efficiency (%)')
yticks(1:0.5:4.5)
ylim([1 4.5])
grid on
legend
title("Modulated Gain, PAE")
hold off

%ACPR
nexttile
hold on
for i = 1:1:numel(freq_mod)
    plot(Pout1_mod(:,i),ACPR_L(:,i),'LineStyle','-','DisplayName',num2str(freq_mod(i)/1e9)+" GHz NDPD")
    plot(Pout1_mod_dpd(:,i),ACPR_L_DPD(:,i),'LineStyle','--','DisplayName',num2str(freq_mod(i)/1e9)+" GHz DPD")
end
xlabel('Pout (dBm)')
ylabel('ACPR Lower (dBc)')
grid on
legend
hold off

nexttile
hold on
for i = 1:1:numel(freq_mod)
    plot(Pout1_mod(:,i),ACPR_U(:,i),'LineStyle','-','DisplayName',num2str(freq_mod(i)/1e9)+" GHz NDPD")
    plot(Pout1_mod_dpd(:,i),ACPR_U_DPD(:,i),'LineStyle','--','DisplayName',num2str(freq_mod(i)/1e9)+" GHz DPD")
end
xlabel('Pout (dBm)')
ylabel('ACPR Upper (dBc)')
grid on
legend
hold off

%EVM
nexttile
hold on
for i = 1:1:numel(freq_mod)
    plot(Pout1_mod(:,i),EVM(:,i),'LineStyle','-','DisplayName',num2str(freq_mod(i)/1e9)+" GHz NDPD")
    plot(Pout1_mod_dpd(:,i),EVM_dpd(:,i),'LineStyle','--','DisplayName',num2str(freq_mod(i)/1e9)+" GHz DPD")
end
xlabel('Pout (dBm)')
ylabel('EVM (%)')
grid on
legend
hold off

%next, choose a pout that has decent efficiency and linearity and plot over
%frequency