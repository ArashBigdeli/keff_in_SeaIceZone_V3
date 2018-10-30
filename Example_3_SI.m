close all;clc,clear

%% Setting inputs
i_c = 0:0.5:100;
U = 5;
%% setting ks
k_w14 = 0.24.*(0.251.*U.^2).*(1-i_c/100);
for i = 1:size(i_c,2)
    k__i(i,1) = keff_SIZ(0.03.*U,U,'I',i_c(i),0,0);
    k__i(i,2) = keff_SIZ(0.02.*U,U,'I',i_c(i),0,0);
    k__i(i,3) = keff_SIZ(0.016.*U,U,'I',i_c(i),0,0);
    k__i(i,4) = keff_SIZ(0.01*U,U,'I',i_c(i),0,0);
    k__i(i,5) = keff_SIZ(0.002.*U,U,'I',i_c(i),0,0);
    k__i(i,6) = keff_SIZ(0.*U,U,'I',i_c(i),0,0);
end

%% plotting
hold on
pl1=plot (i_c,k_w14./k_w14(1),':k','LineWidth',2);
pl2=plot (i_c,k__i(:,1)./k_w14(1),'-r','LineWidth',2);
pl3=plot (i_c,k__i(:,2)./k_w14(1),'-m','LineWidth',2);
pl4=plot (i_c,k__i(:,4)./k_w14(1),'-g','LineWidth',2);
pl5=plot (i_c,k__i(:,6)./k_w14(1),'-b','LineWidth',2);

legend([pl1,pl2,pl3,pl4,pl5],sprintf('Linear relation'),'U_{ice}/U_{10} = 0.03','U_{ice}/U_{10} = 0.02','U_{ice}/U_{10} = 0.01','U_{ice}/U_{10} = 0','location','northeast')
legend boxon;
ylim([-0.4 2.2])
xlim([0 100])

xlabel('Ice cover (%)');
ylabel('K_{eff} / K_{open}','fontsize',10)
title('(\itb)')
