function [eps_ratio,wave_age_nominal] = Mb_coef_cal(Uwind,I_F_W,f,varargin)
% function [eps_ratio,wave_age] = Mb_coef_cal(Uwind,I_F,f,wave_age)
% Inputs :
% Uwind (m/s)
% I_F_W : flag for SI(I) , Fetch(F) or Waveage(W)
%0<f<1 for I_F_W = "I" **
% f in (m) for I_F_W = "f"
% f = cp/u* for I_F_W = W
% Optional : Depth

g = 9.8;

if nargin > 3
    depth_p = varargin{1};
else
    depth_p = 1000;
end
C_D_EC=0.001*(1.1+0.035*Uwind);
C_D = (1.03E-3 + 0.04E-3*Uwind.^1.48)./(Uwind.^0.21);
if strcmpi(I_F_W,'w')
    x1=0;
end
if strcmpi(I_F_W,'F')
    x1 = f;
elseif strcmpi (I_F_W,'I')
    f = f;
    x_hat = (162.*f.^(-0.49));
    x1 = Uwind.^2 .* x_hat ./g; %Smith & Thomson (2016)
end

use_EC = 1; % cosatal protection
use_SP = 0; % shore protection man
use_R10 = 0; % Romero 2010 and KLEISS 2010
use_JONSWAP = 0; % Hasselmann et al. 1973

u_f=Uwind.*sqrt(C_D);

if use_EC == 1 % Costal engineering manual
    u_f_ec=Uwind.*sqrt(C_D_EC);
%     u_f_ec = 1*u_f; % averaging
    if (u_f_ec./g).*7.51E-1.*((g.*x1)./(u_f_ec.^2)).^(1/3) < 2.398E2.*u_f_ec./g; %CEM
        Ts = (u_f_ec./g).*7.51E-1.*(g.*x1./u_f_ec.^2).^(1/3);
    else
        Ts = 2.398E2.*u_f_ec./g;
    end
    if isinf(x1)
        Ts = 2.398E2.*u_f_ec./g;
    end
end

if use_SP == 1 %use shore protection manual
    
    U_A = 0.71 * Uwind^(1.23);
    term1 = 0.833 * ( (g*depth_p)/(U_A^2) )^(3/8);
    T_sh = 7.54 * (U_A/g) *tanh(term1)  * ...
        tanh( (0.0379*(g*x1/Uwind^2)^(1/3) ) / tanh (term1) ); % T shallow
    T0_dp = 2.857E-1 * (g*x1/(U_A^2))^(1/3) * (U_A/g); % T Deep
    T0_fd = 8.134 * U_A/g; % T fully developed
    
    if T0_dp >= T0_fd % check range for deep water
        T_dp = T0_fd;
    else
        T_dp = T0_dp;
    end
    L_estim_dp = (g/pi/2)*T_dp^2 * sqrt(tanh(4*(pi^2)/(T_dp^2).*depth_p/g));
    if depth_p >= 2*L_estim_dp
        Ts = T_dp;
    else
        Ts = T_sh;
    end
    if isinf(x1)
        Ts = T0_fd;
    end
end

if use_R10 == 1 % use Romero 2010
    
x = g.*u_f^(-2)*f;
f_p_hat = 0.28* x ^ -0.25; % dimionsionless peak Freq based on KLEISS 2010
f_p = f_p_hat*g/u_f;
T_p = 1/f_p;
Ts = T_p;
wave_age_foo = (0.57)*(x)^(0.25);
end

if use_JONSWAP ==1
   x_bar = f*g/(Uwind^2);
   f_bar = 3.5 * x_bar^-0.33;
   freq_j = f_bar/(Uwind/g);
   Ts = 1/freq_j;
end

if strcmpi (I_F_W,'I')
%     x1 = Uwind.^2 .* (162.*f.^(-0.49))./9.81; %Smith & Thomson (2016)
freq_hat = 4.0.*(162.*f.^(-0.49))^(-0.33);
fs_i = (g./Uwind).*freq_hat;
Ts_i = 1/fs_i;
Ts = Ts_i.*1.45;
end

% fetch to Wave_age
L = (g/pi/2)*Ts^2 * sqrt(tanh(4*(pi^2)/(Ts^2).*depth_p/g)); % estimate wave length
cp = L/Ts;
wave_age_nominal = cp./u_f;
% wave_age_nominal = min(wave_age_nominal,32); %  upperlimit at 32

if strcmpi(I_F_W,'W')
    wave_age_nominal = f;
end

% cp2cm = 1/1.45;
cp2cm = 1;
wave_age = wave_age_nominal.*cp2cm;%their formulation is not on peak freq
eps_ratio = 7.8E-4.*((wave_age).^3.15);
eps_ratio = max(eps_ratio,1); % the limit is 1 for 9.9670 = cm/u*, first number they have is 2
% sprintf('wave_age = %.2f eps_ratio =%0.2f',wave_age,eps_ratio)

end
%#ok<*BDSCA>
% end