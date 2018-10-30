clf;close all;clear;

%% Setting inputs and compute other gas exchange formulations

C = 0.24;
my_temp = 0;
my_u2 = 0:1:20;
k_n00 = 0.24.*(0.1.*my_u2 + 0.22.*my_u2.^2);
k_h06 = 0.24.*(0.254.*my_u2.^2);
k_s07 = 0.24.*(0.27.*my_u2.^2);
k_w14 = 0.24.*(0.251.*my_u2.^2);

%% Calculating Keff from the main code , at wave ages "W" of 2 to 

for i = 1:size(my_u2,2)

    KL_20(i,1) = keff_SIZ(0,my_u2(i),'W',20,my_temp,my_temp);
    KL_28(i,1) = keff_SIZ(0,my_u2(i),'W',28,my_temp,my_temp);
    KL_32(i,1) = keff_SIZ(0,my_u2(i),'W',32,my_temp,my_temp);
    KL_34(i,1) = keff_SIZ(0,my_u2(i),'W',34,my_temp,my_temp);
    
end

hold on;
p1=plot(my_u2,KL_20,'b','LineWidth',2);
p2=plot(my_u2,KL_28,'c','LineWidth',2);
p3= plot(my_u2,KL_32,'m','LineWidth',2);
p4=plot(my_u2,KL_34,'r','LineWidth',2);
p6=plot(my_u2,k_n00,':k','LineWidth',2);
p7=plot(my_u2,k_h06,'-.k','LineWidth',2);
p8=plot(my_u2,k_s07,'--k','LineWidth',2);
p9=plot(my_u2,k_w14','-k','LineWidth',2);

xlim([0 15])
legend([p1 p2 p3 p4 p6 p7 p8 p9],'cp/u* = 20','cp/u* = 28','cp/u* = 32','cp/u* = 34','N00','H06','S07','W14','Location','NorthWest');
xlabel('U_{10} (ms^{-1})')
ylabel('K (m d^{-1})')
ylim([0 12])
box on
