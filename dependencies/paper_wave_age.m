clf;close all;clear;

scrsz = get(groot,'ScreenSize');
set(groot, 'defaultAxesFontSize', 12)
set(groot,'defaultTextFontSize',12)
set(groot,'defaultLegendFontSize',12)

figure('Position',[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2.5],'visible', 'on')
positionVector1 = [0.06, 0.13, 0.25, 0.8];
positionVector2 = [0.38, 0.13, 0.25, 0.8];
positionVector3 = [0.7, 0.13, 0.25, 0.8];

my_temp = 0;
my_u2 = 0:1:20;
c_d = (1.03E-3 + 0.04E-3*my_u2.^1.48)./(my_u2.^0.21); % COARE3 .Edson et.al (2013)

u_f=my_u2'.*sqrt(c_d');

k_n00 = (0.1.*my_u2 + 0.22.*my_u2.^2);
k_h06 = (0.254.*my_u2.^2);
k_s07 = (0.27.*my_u2.^2);
k_w14 = (0.251.*my_u2.^2);
ind1 = (my_u2 < 3.6); ind2 = ((my_u2 >= 3.6) & my_u2 <= 13);
ind3 = (my_u2 > 13);
k_lis= zeros(size(my_u2));
k_liss_first_part = 0.17.*my_u2;
k_lis(ind1) = 0.17.*my_u2(ind1);
k_lis(ind2) = 2.85.*my_u2(ind2) - 9.65;
k_lis(ind3) = 5.9.*my_u2(ind3) - 49.3;

for i = 1:size(my_u2,2);
    
    KL_2(i,1) = (1/0.24).*keff_SIZ(0,my_u2(i),'W',2,my_temp,my_temp);
    KL_4(i,1) = (1/0.24).*keff_SIZ(0,my_u2(i),'W',4,my_temp,my_temp);
    KL_8(i,1) = (1/0.24).*keff_SIZ(0,my_u2(i),'W',8,my_temp,my_temp);
    KL_6(i,1) = (1/0.24).*keff_SIZ(0,my_u2(i),'W',6,my_temp,my_temp);
    KL_12(i,1) = (1/0.24).*keff_SIZ(0,my_u2(i),'W',12,my_temp,my_temp);
    KL_14(i,1) = (1/0.24).*keff_SIZ(0,my_u2(i),'W',14,my_temp,my_temp);
    KL_17(i,1) = (1/0.24).*keff_SIZ(0,my_u2(i),'W',17,my_temp,my_temp);
    KL_20(i,1) = (1/0.24).*keff_SIZ(0,my_u2(i),'W',20,my_temp,my_temp);
    KL_23(i,1) = (1/0.24).*keff_SIZ(0,my_u2(i),'W',23,my_temp,my_temp);
    KL_26(i,1) = (1/0.24).*keff_SIZ(0,my_u2(i),'W',26,my_temp,my_temp);
    KL_28(i,1) = (1/0.24).*keff_SIZ(0,my_u2(i),'W',28,my_temp,my_temp);
    KL_32(i,1) = (1/0.24).*keff_SIZ(0,my_u2(i),'W',32,my_temp,my_temp);
    KL_34(i,1) = (1/0.24).*keff_SIZ(0,my_u2(i),'W',34,my_temp,my_temp);
    KL_38(i,1) = (1/0.24).*keff_SIZ(0,my_u2(i),'W',38,my_temp,my_temp);
    KL_40(i,1) = (1/0.24).*keff_SIZ(0,my_u2(i),'W',40,my_temp,my_temp);
    KL_75(i,1) = (1/0.24).*keff_SIZ(0,my_u2(i),'W',75,my_temp,my_temp);
    
    
end


subplot ('Position',positionVector1)
hold on
my_u1 = 0:0.1:20;
my_fetch = 100:10:20E3;
 for i = 1:length(my_u1)
for j = 1:length(my_fetch)
[my_epsr(i,j),my_wa(i,j)]= Mb_coef_cal(my_u1(i),'f',my_fetch(j));
end
 end
pcolor(my_u1,my_fetch./1000,my_wa')
colormap('jet')
xlabel('U_{10} (ms^{-1})')
ylabel('Fetch (km)')
xlim([1 12])
caxis([8 34]);
shading flat
h=colorbar('horizental');
ylabel(h,'wave age')
title('(\ita)')

subplot ('Position',positionVector2)
hold on
C = 0.24;

p1=plot(my_u2,C.*KL_20,'b','LineWidth',2);
% plot(my_u2,KL_26,'g','LineWidth',2);
p2=plot(my_u2,C.*KL_28,'c','LineWidth',2);
p3= plot(my_u2,C.*KL_32,'m','LineWidth',2);
p4=plot(my_u2,C.*KL_34,'r','LineWidth',2);
p5=plot(my_u2,C.*k_lis','^-k','LineWidth',1,'Markersize',5);
p6=plot(my_u2,C.*k_n00,':k','LineWidth',2);
p7=plot(my_u2,C.*k_h06,'-.k','LineWidth',2);
p8=plot(my_u2,C.*k_s07,'--k','LineWidth',2);
p9=plot(my_u2,C.*k_w14','-k','LineWidth',2);

xlim([0 15])
legend([p1 p2 p3 p4 p5 p6 p7 p8 p9],'cp/u* = 20','cp/u* = 28','cp/u* = 32','cp/u* = 34','LM','N00','H06','S07','W14','Location','NorthWest');
xlabel('U_{10} (ms^{-1})')
ylabel('K (m d^{-1})')
ylim([0 12])
title('(\itb)')
box on
% ,'cp/u* = 26'
% %
load('../data/Zappa/zappa2007.mat');
subplot ('Position',positionVector3)
hold on
plot(my_u2,C.*KL_12,'b','LineWidth',2);
plot(my_u2,C.*KL_40,'g','LineWidth',2);
plot(my_u2,C.*KL_75,'--r','LineWidth',2);
plot(my_u2,C.*k_w14','-k','LineWidth',2);
ylabel('K (m d^{-1})')
plot(Zappa_u,C.*Zappa_k,'^k','Markersize',6);
title('(\itc)')
% colorbar
xlabel('U_{10} (ms^{-1})')
xlim([0 12])
ylim([0 12])

% set(gca,'YTickLabel',[]);
% xlabel('U_{10} (ms^{-1})')
box on
legend('cp/u* = 12','cp/u* = 40','cp/u* = 75','W14','Zappa 2007','Location','NorthWest');
% ~ 12 to 40+
