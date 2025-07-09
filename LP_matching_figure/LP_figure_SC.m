close all
clc

LP_targets = readtable("LP_Targets_36_42_GHz.csv");
f_range = 36:1:42;
LP_Z_str = strings(1,numel(f_range));
for i = 1:numel(f_range)
    inter1 = LP_targets(:,i).Variables;
    LP_Z_str(i) = inter1{1};
end
LP_Z = zeros(1,numel(LP_Z_str));
for i = 1:numel(LP_Z_str)
    str1 = split(LP_Z_str(i)," ");
    str2 = str1(3);
    imagVal = erase(str2,"j");
    imagVal = imagVal+"i";
    if str1(2) == "+"
        cmplxZ = str2double(str1(1)) + str2double(imagVal);
    else
        cmplxZ = str2double(str1(1)) - str2double(imagVal);
    end
    LP_Z(i) = cmplxZ;
end
%%

OMN_S11 = sparameters("OMN_traj.s2p");
OMN_S11 = OMN_S11.Parameters(1,1,:);
OMN_S11 = OMN_S11(:);

%%
%finally, plot these impedances on SC
close all
plotsc(LP_Z,'Domain','Z','Z0',50,'Marker',"o",'MarkerSize',15,'LineStyle','none','LineWidth',4)
hold on
plotsc(OMN_S11,'Domain','G','Z0',50,'LineStyle','-','LineWidth',6)
legend('LP','OMN')
% pbaspect([1 0.5 1])

% fH = gcf;
% aH = fH.Children;
% pbaspect(aH(2),[0.5 1 1])
% daspect(aH(2),[0.5 1 1])
