clf;clear;close all;

%% Setting the inputs
my_temp = 0;
my_u2 = 0:0.5:20;
k_w14 = 0.24.*(0.251.*my_u2.^2);

%% Running the code for I_W_F = 'f' and fetch limit of 5, 20, 50 km 
for i = 1:length(my_u2)
    KL_05000(i,1) = keff_SIZ(0,my_u2(i),'f',5E3,my_temp,my_temp);
    KL_20000(i,1) = keff_SIZ(0,my_u2(i),'f',20E3,my_temp,my_temp);
    KL_50000(i,1) = keff_SIZ(0,my_u2(i),'f',50E3,my_temp,my_temp);
end

%% plotting
hold on
p1=plot(my_u2,KL_05000,'b','LineWidth',2);
p2=plot(my_u2,KL_20000,'m','LineWidth',2);
p3=plot(my_u2,KL_50000,'r','LineWidth',2);
p4=plot(my_u2,k_w14,'--k','LineWidth',2);
legend([p1 p2 p3 p4 ],'X= 5000 m','X = 20000 m', 'X = 50000 m','W14','location','northwest')
xlabel('U_{10} (ms^{-1})')
ylabel('K (m d^{-1})')
xlim([0 12]);
ylim([0 4])
box on
