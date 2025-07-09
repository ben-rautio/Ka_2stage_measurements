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
output_off_avg = output_off_avg(output_off_avg~=0);
output_off_avg = output_off_avg - output_adapter_loss;
inp_off_avg = inp_off_avg(inp_off_avg~=0);
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

FundTone_out_cw = FundTone_out_cw(:,1:numel(freq_cw));
FundTone_outFreq_cw = FundTone_outFreq_cw(:,1:numel(freq_cw));
FundTone_outFreq2_cw = FundTone_outFreq2_cw(:,1:numel(freq_cw));

Pout1 = FundTone_out_cw + Attenuation_scaled; % Use whichever measured power you want. Or all of them.
Pout2 = FundTone_outFreq_cw + Attenuation_scaled; % Use whichever measured power you want. Or all of them.
Pout3 = FundTone_outFreq2_cw + Attenuation_scaled; % Use whichever measured power you want. Or all of them.

Pout_watts1 = 10.^(Pout1./10)./1000;
Pout_watts2 = 10.^(Pout2./10)./1000;
Pout_watts3 = 10.^(Pout3./10)./1000;

FundTone_In_cw = FundTone_In_cw(:,1:numel(freq_cw));
FundTone_InFreq_cw = FundTone_InFreq_cw(:,1:numel(freq_cw));
FundTone_InFreq2_cw = FundTone_InFreq2_cw(:,1:numel(freq_cw));

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
Pdc_tot = Pdc_tot(:,1:numel(freq_cw));

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
pbaspect([1 0.5 0.5])
plot(pdel_ADS,Gain_ADS,'LineStyle','--','DisplayName','Sim Gain')
plot(Pout1(:,11),Gain1(:,11),'LineStyle','-','DisplayName','Meas. Gain')
% for i = 1:length(freq_cw)-1
%     plot(Pout1(:,i),Gain1(:,i),'Color','r','DisplayName',"Gain at "+num2str(freq_cw(i)/1e9)+" GHz")
% end
ylabel('Gain (dB)')
yticks(0:2:20)
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
legend("Location","best")
title("CW Measurements vs. Simulated (40 GHz)")
saveas(gcf,fullfile("figures","pdf","CW_Measurements_vs_Simulated(40GHz).pdf"))
saveas(gcf,fullfile("figures","svg","CW_Measurements_vs_Simulated(40GHz).svg"))
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

figure
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
h=legend("Location","best");
set(h,'Interpreter','latex')
title("S-Parameters Meas. and Sim")
hold off
pbaspect([1 0.5 0.5])
saveas(gcf,fullfile("figures","pdf","Sparams_vs_Simulated.pdf"))
saveas(gcf,fullfile("figures","svg","Sparams_vs_Simulated(40GHz).svg"))

%% Modulated Measurements
close all
%get all the data from the modulated meas workspace
freq_mod = load("PowerMeter_modulated_3-5-2025.mat","Fo").Fo;
Pdes_mod = load("PowerMeter_modulated_3-5-2025.mat","Pdes").Pdes;

FundTone_Out_mod = load("PowerMeter_modulated_3-5-2025.mat","FundTone_Out").FundTone_Out;
CW_output_offset_mod = repmat(output_off_avg, length(Pdes_mod), 1);
Pout1_mod = FundTone_Out_mod + CW_output_offset_mod;
Pout1_watts_mod = 10.^(Pout1_mod./10)./1000;

FundTone_In_mod = load("PowerMeter_modulated_3-5-2025.mat","FundTone_In").FundTone_In;
CW_input_offset_mod = repmat(inp_off_avg, length(Pdes_mod), 1);
Pin1_mod = FundTone_In_mod + CW_input_offset_mod;
Mod_Pin1_watts = 10.^(Pin1_mod./10)./1000;

Mod_PDC = load("PowerMeter_modulated_3-5-2025.mat","Pdc_tot").Pdc_tot;

DE1_mod = (Pout1_watts_mod./Mod_PDC).*100;
PAE1_mod = ((Pout1_watts_mod-Mod_Pin1_watts)./Mod_PDC).*100;
Gain1_mod = Pout1_mod - Pin1_mod;

% Pout1_mod = load("PowerMeter_modulated_3-5-2025.mat","Pout1").Pout1;
% Gain1_mod = load("PowerMeter_modulated_3-5-2025.mat","Gain1").Gain1;
% PAE1_mod = load("PowerMeter_modulated_3-5-2025.mat","PAE1").PAE1;
% Pout1_mod = Pout1_mod_s.Pout1;
% Pout1_mod_dpd_s = load("Modulated_Meas_WS.mat","Pout1_dpd");
% Pout1_mod_dpd = Pout1_mod_dpd_s.Pout1_dpd;
% Gain1_mod_s = load("Modulated_Meas_WS.mat","Gain1");
% Gain1_mod = Gain1_mod_s.Gain1;
% Gain1_mod_dpd_s = load("Modulated_Meas_WS.mat","Gain1_dpd");
% Gain1_mod_dpd = Gain1_mod_dpd_s.Gain1_dpd;
% freq_mod_s = load("Modulated_Meas_WS.mat","Fo");
% freq_mod = freq_mod_s.Fo;
% PAE1_mod_s = load("Modulated_Meas_WS.mat","PAE1");
% PAE1_mod = PAE1_mod_s.PAE1;
% PAE1_mod_dpd_s = load("Modulated_Meas_WS.mat","PAE1_dpd");
% PAE1_mod_dpd = PAE1_mod_dpd_s.PAE1_dpd;

ACPR_L_S = load("FSW_modulated_3-5-2025.mat","adjpow_lower");
ACPR_L = ACPR_L_S.adjpow_lower;
% ACPR_L_DPD_S = load("Modulated_Meas_WS.mat","adjpow_lower_dpd");
% ACPR_L_DPD = ACPR_L_DPD_S.adjpow_lower_dpd;

ACPR_U_S = load("FSW_modulated_3-5-2025.mat","adjpow_upper");
ACPR_U = ACPR_U_S.adjpow_upper;
% ACPR_U_DPD_S = load("Modulated_Meas_WS.mat","adjpow_upper_dpd");
% ACPR_U_DPD = ACPR_U_DPD_S.adjpow_upper_dpd;

EVM_s = load("FSW_modulated_3-5-2025.mat","EVM");
EVM = EVM_s.EVM;
% EVM_dpd_s = load("Modulated_Meas_WS.mat","EVM_dpd");
% EVM_dpd = EVM_dpd_s.EVM_dpd;

CFAC_out = load("FSW_modulated_3-5-2025.mat","CFAC_out").CFAC_out;

% figure
% tiledlayout("flow")

for i = 1:1:numel(freq_mod)
    %saving figures
    % nexttile
    figure
    hold on
    set(gcf,'Color','w')
    box on
    title("Gain, Avg. PAE, Output CF at "+num2str(freq_mod(i)/1e9)+ " GHz")
    fname = "Modulated_Gain_PAE_"+num2str(freq_mod(i)/1e9)+ "_GHz";
    yyaxis left
    ylim([0 25])
    yticks(0:2.5:25)
    xlim([21 29])
    xlabel('Pout (dBm)')
    ylabel('Gain (dB), Avg. PAE (%)')
    dname1 = "Gain";
    plot(Pout1_mod(:,i),Gain1_mod(:,i),'LineStyle','-','DisplayName',dname1, 'LineWidth',2,'Color',[0 0.4470 0.7410])
    dname2 = "Avg. PAE";
    plot(Pout1_mod(:,i),PAE1_mod(:,i),'LineStyle','-','DisplayName',dname2,'LineWidth',2,'Color',[0.8500 0.3250 0.0980])
    ax = gca;
    ax.YColor = "k";

    yyaxis right
    plot(Pout1_mod(:,i),CFAC_out(:,i),'LineStyle','-','DisplayName',"Output CF",'LineWidth',2,'Color',[0.4940 0.1840 0.5560])
    grid on
    legend("Location","best")
    ylim([0 10])
    yticks(0:1:10)
    ylabel('Output CF (dB)')
    ax = gca;
    ax.YColor = "k";
    pbaspect([1 0.5 0.5])
    hold off
    saveas(gcf,fullfile("figures","pdf",fname+".pdf"))
    saveas(gcf,fullfile("figures","svg",fname+".svg"))
end

%ACPR
% figure
% tiledlayout("flow")

for i = 1:1:numel(freq_mod)
    % nexttile
    figure
    hold on
    set(gcf,'Color','w')
    box on
    title("ACPR "+num2str(freq_mod(i)/1e9)+" GHz")
    fname = "ACPR_"+num2str(freq_mod(i)/1e9)+"_GHz";
    plot(Pout1_mod(:,i),ACPR_L(:,i),'Color','r','DisplayName',"ACPR Lower "+num2str(freq_mod(i)/1e9)+" GHz NDPD")
    plot(Pout1_mod(:,i),ACPR_U(:,i),'Color','b','DisplayName',"ACPR Upper "+num2str(freq_mod(i)/1e9)+" GHz NDPD")
    xlabel('Pout (dBm)')
    ylabel('ACPR Lower (dBc)')
    xlim([18 30])
    ylim([-40 -15])
    yticks(-40:5:-15)
    grid on
    legend("Location","best")
    pbaspect([1 0.5 0.5])
    hold off
    saveas(gcf,fullfile("figures","pdf",fname+".pdf"))
    saveas(gcf,fullfile("figures","svg",fname+".svg"))
end

%EVM
% figure
% tiledlayout("flow")

% for i = 1:1:numel(freq_mod)
%     % nexttile
%     figure
%     hold on
%     title("EVM "+num2str(freq_mod(i)/1e9)+" GHz")
%     fname = "EVM_"+num2str(freq_mod(i)/1e9)+"GHz";
%     plot(Pout1_mod(:,i),EVM(:,i),'LineStyle','-','DisplayName',"EVM "+num2str(freq_mod(i)/1e9)+" GHz NDPD")
%     xlabel('Pout (dBm)')
%     ylabel('EVM (%)')
%     xlim([18 30])
%     ylim([0 25])
%     yticks(0:5:25)
%     grid on
%     legend("Location","best")
%     hold off
%     pbaspect([1 0.5 0.5])
%     saveas(gcf,fullfile("figures","pdf",fname+".pdf"))
%     saveas(gcf,fullfile("figures","svg",fname+".svg"))
% end

%CFAC
% figure
% tiledlayout("flow")
% 
% for i = 1:1:numel(freq_mod)
%     % nexttile
%     figure
%     hold on
%     title("CFAC Out "+num2str(freq_mod(i)/1e9)+" GHz")
%     plot(Pout1_mod(:,i),CFAC_out(:,i),'LineStyle','-','DisplayName',"CFAC Out "+num2str(freq_mod(i)/1e9)+" GHz NDPD")
%     xlabel('Pout (dBm)')
%     ylabel('CFAC Out (dB)')
%     xlim([18 30])
%     ylim([0 10])
%     yticks(0:1:10)
%     grid on
%     legend
%     pbaspect([1 0.5 0.5])
%     hold off
% end

%% plot Coupled power vs. SMW power for modulated and CW measurements for sanity check
% close all
% Pcoupled_cw_s = load("driveup_38_40GHz.mat","FundTone_In");
% Pcoupled_cw = Pcoupled_cw_s.FundTone_In;
% SMW_power_cw_s = load("driveup_38_40GHz.mat","Pavs");
% SMW_power_cw = SMW_power_cw_s.Pavs;
% Pcoupled_cw = Pcoupled_cw(:,1:size(SMW_power_cw,2));
% Pcoupled_mod_s = load("PowerMeter_modulated_3-5-2025.mat","FundTone_In");
% Pcoupled_mod = Pcoupled_mod_s.FundTone_In;
% SMW_power_mod_s = load("PowerMeter_modulated_3-5-2025.mat","Pavs");
% SMW_power_mod = SMW_power_mod_s.Pavs;
% 
% OutputCoupled_cw = load("driveup_38_40GHz.mat","FundTone_Out").FundTone_Out;
% OutputCoupled_cw = OutputCoupled_cw(:,1:size(SMW_power_cw,2));
% 
% figure
% tiledlayout("flow")
% for i = 1:length(freq_cw)
%     nexttile
%     hold on
%     xlabel("SMW Power (dBm)")
%     ylabel("Input Coupled Power")
%     title("Coupled Power at " + num2str(freq_cw(i) / 1e9) + " GHz")
%     plot(SMW_power_cw(:,i),Pcoupled_cw(:,i),'LineStyle','-','DisplayName','CW')
%     plot(SMW_power_mod(:,i),Pcoupled_mod(:,i),'LineStyle','--','DisplayName','Mod')
%     grid on
%     legend
%     pbaspect([1 0.5 0.5])
%     hold off
% end
% 
% figure
% tiledlayout("flow")
% for i = 1:length(freq_cw)
%     nexttile
%     hold on
%     xlabel("SMW Power (dBm)")
%     ylabel("Output Power")
%     title("Output Power at " + num2str(freq_cw(i) / 1e9) + " GHz")
%     plot(SMW_power_cw(:,i),OutputCoupled_cw(:,i),'LineStyle','-','DisplayName','CW')
%     plot(SMW_power_mod(:,i),FundTone_Out_mod(:,i),'LineStyle','--','DisplayName','Mod')
%     ylim([-10 5])
%     yticks(-10:1:5)
%     grid on
%     legend
%     pbaspect([1 0.5 0.5])
%     hold off
% end

%% plot linearity over freq for a given Pout
close all
pout_des = 24;
EVM_freq = zeros(1, numel(freq_mod));
ACPR_L_freq = zeros(1, numel(freq_mod));
ACPR_U_freq = zeros(1, numel(freq_mod));
PAE_freq = zeros(1, numel(freq_mod));
PAE_freq_CW = zeros(1, numel(freq_mod));
Gain_freq = zeros(1, numel(freq_mod));
Pout_freq = zeros(1, numel(freq_mod));
tol = 0.3;
for i = 1:numel(freq_mod)
    col = Pout1_mod(:,i);
    col_P = abs((col - pout_des)) < tol;
    if sum(col_P) > 1
        diff_vector = abs((col - pout_des));
        smallest = min(diff_vector);
        col_P = abs(diff_vector - smallest) < 1e-6;
        col_ACPR_L = ACPR_L(:,i);
        ACPR_L_freq(i) = col_ACPR_L(col_P);
        col_ACPR_U = ACPR_U(:,i);
        ACPR_U_freq(i) = col_ACPR_U(col_P);
        col_EVM = EVM(:,i);
        EVM_freq(i) = col_EVM(col_P);
        col_PAE = PAE1_mod(:,i);
        PAE_freq(i) = col_PAE(col_P);
        col_Gain = Gain1_mod(:,i);
        Gain_freq(i) = col_Gain(col_P);
        col_Pout = Pout1_mod(:,i);
        Pout_freq(i) = col_Pout(col_P);
    elseif sum(col_P) == 1
        col_ACPR_L = ACPR_L(:,i);
        ACPR_L_freq(i) = col_ACPR_L(col_P);
        col_ACPR_U = ACPR_U(:,i);
        ACPR_U_freq(i) = col_ACPR_U(col_P);
        col_EVM = EVM(:,i);
        EVM_freq(i) = col_EVM(col_P);
        col_PAE = PAE1_mod(:,i);
        PAE_freq(i) = col_PAE(col_P);
        col_Gain = Gain1_mod(:,i);
        Gain_freq(i) = col_Gain(col_P);
        col_Pout = Pout1_mod(:,i);
        Pout_freq(i) = col_Pout(col_P);
    else
        % tol=tol+0.1;
        error("no Pout within allowed tolerance")
    end
end

%plot it
figure
% tiledlayout("flow")
% nexttile
hold on
title("EVM at " + num2str(pout_des) + " dBm Output Power")
fname = "EVMvsFreq_at_" + num2str(pout_des) + "_dBm_Pout";
plot(freq_mod ./ 1e9,EVM_freq)
grid on
xlabel("Frequency (GHz)")
ylabel("EVM (%)")
set(gcf,'Color','w')
box on
pbaspect([1 0.5 0.5])
hold off
saveas(gcf,fullfile("figures","pdf",fname+".pdf"))
saveas(gcf,fullfile("figures","svg",fname+".svg"))

% nexttile
figure
hold on
title("ACPR at " + num2str(pout_des) + " dBm Output Power")
fname = "ACPRvsFreq_at_" + num2str(pout_des) + "_dBm_Pout";
plot(freq_mod ./ 1e9,ACPR_U_freq,'DisplayName','ACPR Upper')
plot(freq_mod ./ 1e9,ACPR_L_freq,'DisplayName','ACPR Lower')
grid on
xlabel("Frequency (GHz)")
ylabel("ACPR (dBc)")
legend("Location","best")
set(gcf,'Color','w')
box on
pbaspect([1 0.5 0.5])
hold off
saveas(gcf,fullfile("figures","pdf",fname+".pdf"))
saveas(gcf,fullfile("figures","svg",fname+".svg"))

% nexttile
figure
hold on
title("Gain, PAE at " + num2str(pout_des) + " dBm Output Power")
fname = "Gain_PAEvsFreq_at_" + num2str(pout_des) + "_dBm_Pout";
% yyaxis left
plot(freq_mod ./ 1e9,Gain_freq,'DisplayName','Gain','LineWidth',2,'Color','r')
ylim([5 20]);
yticks(5:1:20);
ylabel("Gain (dB), PAE (%)")
% yyaxis right
plot(freq_mod ./ 1e9,PAE_freq,'DisplayName','PAE','LineWidth',2,'Color','b')
% ylim([8 20]);
% yticks(8:2:20);
% yticks
grid on
xlabel("Frequency (GHz)")
% ylabel("PAE (%)")
legend("Location","best")
set(gcf,'Color','w')
box on
pbaspect([1 0.5 0.5])
hold off
saveas(gcf,fullfile("figures","pdf",fname+".pdf"))
saveas(gcf,fullfile("figures","svg",fname+".svg"))
%% Do the same for CW measurements

close all
pout_des = 24;
% EVM_freq = zeros(1, numel(freq_mod));
% ACPR_L_freq = zeros(1, numel(freq_mod));
% ACPR_U_freq = zeros(1, numel(freq_mod));
% PAE_freq = zeros(1, numel(freq_mod));
PAE_freq = zeros(1, numel(freq_cw));
Gain_freq = zeros(1, numel(freq_cw));
Pout_freq = zeros(1, numel(freq_cw));
tol = 0.42;
for i = 1:numel(freq_cw)
    col = Pout1(:,i);
    col_P = abs((col - pout_des)) < tol;
    if sum(col_P) > 1
        diff_vector = abs((col - pout_des));
        smallest = min(diff_vector);
        col_P = abs(diff_vector - smallest) < 1e-6;
        % col_ACPR_L = ACPR_L(:,i);
        % ACPR_L_freq(i) = col_ACPR_L(col_P);
        % col_ACPR_U = ACPR_U(:,i);
        % ACPR_U_freq(i) = col_ACPR_U(col_P);
        % col_EVM = EVM(:,i);
        % EVM_freq(i) = col_EVM(col_P);
        col_PAE = PAE1(:,i);
        PAE_freq(i) = col_PAE(col_P);
        col_Gain = Gain1(:,i);
        Gain_freq(i) = col_Gain(col_P);
        col_Pout = Pout1(:,i);
        Pout_freq(i) = col_Pout(col_P);
    elseif sum(col_P) == 1
        % col_ACPR_L = ACPR_L(:,i);
        % ACPR_L_freq(i) = col_ACPR_L(col_P);
        % col_ACPR_U = ACPR_U(:,i);
        % ACPR_U_freq(i) = col_ACPR_U(col_P);
        % col_EVM = EVM(:,i);
        % EVM_freq(i) = col_EVM(col_P);
        col_PAE = PAE1(:,i);
        PAE_freq(i) = col_PAE(col_P);
        col_Gain = Gain1(:,i);
        Gain_freq(i) = col_Gain(col_P);
        col_Pout = Pout1(:,i);
        Pout_freq(i) = col_Pout(col_P);
    else
        % tol=tol+0.1;
        error("no Pout within allowed tolerance")
    end
end

%plot it

% nexttile
figure
hold on
title("CW Gain, PAE at " + num2str(pout_des) + " dBm Output Power")
fname = "Gain_PAEvsFreq_at_" + num2str(pout_des) + "_dBm_Pout_CW";
plot(freq_cw ./ 1e9,Gain_freq,'DisplayName','Gain','LineWidth',2,'Color','r')
ylim([9 21]);
yticks(9:1:21);
ylabel("Gain (dB), PAE (%)")
plot(freq_cw ./ 1e9,PAE_freq,'DisplayName','PAE','LineWidth',2,'Color','b')
grid on
xlabel("Frequency (GHz)")
legend("Location","best")
set(gcf,'Color','w')
box on
pbaspect([1 0.5 0.5])
hold off
saveas(gcf,fullfile("figures","pdf",fname+".pdf"))
saveas(gcf,fullfile("figures","svg",fname+".svg"))

%% output coupler power


%apply CW offset to Modulated measurement readings

% output powers
% Mod_Tone_Out1_s = load("Modulated_Meas_WS.mat","FundTone_Out");
% Mod_Tone_Out1 = Mod_Tone_Out1_s.FundTone_Out;
% Mod_Tone_Out2_s = load("Modulated_Meas_WS.mat","FundTone_OutFreq");
% Mod_Tone_Out2 = Mod_Tone_Out2_s.FundTone_OutFreq;
% Mod_Tone_Out1_dpd_s = load("Modulated_Meas_WS.mat","FundTone_Out_dpd");
% Mod_Tone_Out1_dpd = Mod_Tone_Out1_dpd_s.FundTone_Out_dpd;
% Mod_Tone_Out2_dpd_s = load("Modulated_Meas_WS.mat","FundTone_OutFreq_dpd");
% Mod_Tone_Out2_dpd = Mod_Tone_Out2_dpd_s.FundTone_OutFreq_dpd;
%DC powers
% Mod_PDC_s = load("Modulated_Meas_WS.mat","Pdc_tot");
% Mod_PDC = Mod_PDC_s.Pdc_tot;
%input powers
%no dpd
% Mod_Pin1 = load("Modulated_Meas_WS.mat","Pin1").Pin1;
% Mod_Pin1_watts = load("Modulated_Meas_WS.mat","Pin1_watts").Pin1_watts;
% Mod_Pin2 = load("Modulated_Meas_WS.mat","Pin2").Pin2;
% Mod_Pin2_watts = load("Modulated_Meas_WS.mat","Pin2_watts").Pin2_watts;
% %dpd
% Mod_Pin1_dpd = load("Modulated_Meas_WS.mat","Pin1_dpd").Pin1_dpd;
% Mod_Pin1_watts_dpd = load("Modulated_Meas_WS.mat","Pin1_watts_dpd").Pin1_watts_dpd;
% Mod_Pin2_dpd = load("Modulated_Meas_WS.mat","Pin2_dpd").Pin2_dpd;
% Mod_Pin2_watts_dpd = load("Modulated_Meas_WS.mat","Pin2_watts_dpd").Pin2_watts_dpd;

%CW output atten
% CW_offset_s = load("Modulated_Meas_WS.mat","Attenuation_scaled");
% CW_offset = CW_offset_s.Attenuation_scaled;
% %no dpd
% Pout1_mod = Mod_Tone_Out1 + CW_offset;
% Pout1_watts_mod = 10.^(Pout1_mod./10)./1000;
% Pout2_mod = Mod_Tone_Out2 +  CW_offset;
% Pout2_watts_mod = 10.^(Pout2_mod./10)./1000;
% %dpd
% Pout1_dpd_mod = Mod_Tone_Out1_dpd + CW_offset;
% Pout1_watts_dpd_mod = 10.^(Pout1_dpd_mod./10)./1000;
% Pout2_dpd_mod = Mod_Tone_Out2_dpd + CW_offset;
% Pout2_watts_dpd_mod = 10.^(Pout2_dpd_mod./10)./1000;

%calc effic and gain with cw offset
%no dpd
% DE1_mod = (Pout1_watts_mod./Mod_PDC).*100;
% PAE1_mod = ((Pout1_watts_mod-Mod_Pin1_watts)./Mod_PDC).*100;
% Gain1_mod = Pout1_mod - Mod_Pin1;
% 
% DE2_mod = (Pout2_watts_mod./Mod_PDC).*100;
% PAE2_mod = ((Pout2_watts_mod-Mod_Pin2_watts)./Mod_PDC).*100;
% Gain2_mod = Pout2_mod - Mod_Pin2;
% %dpd
% DE1_mod_dpd = (Pout1_watts_dpd_mod./Mod_PDC).*100;
% PAE1_mod_dpd = ((Pout1_watts_dpd_mod-Mod_Pin1_watts_dpd)./Mod_PDC).*100;
% Gain1_mod_dpd = Pout1_dpd_mod - Mod_Pin1_dpd;
% 
% DE2_mod_dpd = (Pout2_watts_dpd_mod./Mod_PDC).*100;
% PAE2_mod_dpd = ((Pout2_watts_dpd_mod-Mod_Pin2_watts_dpd)./Mod_PDC).*100;
% Gain2_mod_dpd = Pout2_dpd_mod - Mod_Pin2_dpd;
%plot all this stuff

% figure
% tiledlayout("flow")
% for i = 1:1:numel(freq_cw)
%     nexttile
%     title("Mod Gain, PAE at " + num2str(freq_cw(i) / 1e9) + " GHz")
%     hold on
%     yyaxis left
%     xlabel("Pout (dBm)")
%     ylabel("Gain")
%     plot(Pout1_mod(:,i),Gain1_mod(:,i),'LineStyle','-','DisplayName','Gain NDPD')
%     plot(Pout1_dpd_mod(:,i),Gain1_mod_dpd(:,i),'LineStyle','--','DisplayName','Gain DPD')
%     % yticks(8:1:14)
%     % ylim([8 14])
%     yyaxis right
%     % yticks(1:0.5:4)
%     % ylim([1 4])
%     plot(Pout1_mod(:,i),PAE1_mod(:,i),'LineStyle','-','DisplayName','PAE NDPD')
%     plot(Pout1_dpd_mod(:,i),PAE1_mod_dpd(:,i),'LineStyle','--','DisplayName','PAE DPD')
%     ylabel("PAE (%)")
%     grid on
%     legend
%     hold off
% end

%% drain efficiency of stage 2
% close all
% Pdc_tot2=Pdc_tot2(:,1:numel(freq_cw));
% DE_s2_v1 = (Pout_watts1./Pdc_tot2(:,1:numel(freq_cw))).*100;
% %plot at 40 GHz, over output power
% hold on; box on; grid on; grid minor; set(gcf,'color','w')
% yyaxis left
% ylim([0 20])
% % plot(pdel_ADS,Gain_ADS,'LineStyle','--','DisplayName','Sim Gain')
% plot(Pout1(:,11),Gain1(:,11),'LineStyle','-','DisplayName','Meas. Gain')
% % for i = 1:length(freq_cw)-1
% %     plot(Pout1(:,i),Gain1(:,i),'Color','r','DisplayName',"Gain at "+num2str(freq_cw(i)/1e9)+" GHz")
% % end
% ylabel('Gain (dB)')
% % ylim([4 13])
% xlim([16 30])
% yyaxis right
% ylim([0 36])
% yticks(0:3:36)
% % plot(pdel_ADS,PAE_ADS,'LineStyle','--','DisplayName','Sim PAE')
% plot(Pout1(:,11),DE_s2_v1(:,11),'LineStyle','-','DisplayName','Meas. DE')
% % for i = 1:length(freq_cw)-1
% %     plot(Pout1(:,i),PAE1(:,i),'Color','b','DisplayName',"PAE at "+num2str(freq_cw(i)/1e9)+" GHz")
% % end
% ylabel('DE (%)')
% xlabel('Pout (dBm)')
% % legend("Sim Gain", "Meas Gain1", "Meas Gain2","Meas Gain3","Sim PAE", "Meas PAE1", "Meas PAE2","Meas PAE3")
% legend("Location","best")
% title("DE of 2nd Stage vs. Output Power (40 GHz)")
% pbaspect([1 0.5 0.5])
% saveas(gcf,fullfile("figures","pdf","DE_2ndStage_vs_pout(40GHz).pdf"))
% saveas(gcf,fullfile("figures","svg","DE_2ndStage_vs_pout(40GHz).svg"))
% hold off