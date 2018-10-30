function [keff,eps_f,eps_iw,eps_aw] = keff_SIZ(Uice,Uwind,I_F_W,f,waterT,airT,varargin)

%function  [keff,eps_f,eps_iw,eps_aw] = keff_SIZ(Uice,Uwind,I_F_W,f,waterT,airT,varargin)
%
%
% 
% INPUT:
%       Uice = ice velocity (m/s)
%       Uwind = wind speed (not vector) in m/s
%       I_F_W = Ice or Fetch or WaveAge flag 
%       "I" = Sic , "F" = Fetch , 'W' = Waveage. Set f accordingly
%       I_F_W == "I" => f : sic = sea ice concentration (%)
%       I_F_W == "F" => f : fetch = effective fetch in (m)
%       I_F_W == "W" => f : = wave age = cp/u*
%       waterT = water temp in deg C
%       airT = air temp in deg C.
%       
%       Optional Inputs
%       depth , default is 1000 ; Shallow or deep water range in wave speed
%       Humidity ,default calculated based on temperature
%       Mixed layer depth (m) , default 30m
%       ice thickness in (m) defualt 2m
%
% OUTPUT:
%
% keff (m/d): The effective gas transfer velocity 
%
%
% DEPENDENCIES:
%   These routines require the Gibbs Seawater toolbox:
%   http://www.teos-10.org/
%
% AUTHOR:  
%          Arash Bigdeli (arash.big@gmail.com)
%          Brice Loose (osoberlice@gmail.com)    
% REFERENCE:

%	Bigdeli et al., (2017), "Wave attenuation and gas exchange
%	velocity in marginal sea ice zone." 
%	https://doi.org/10.1002/2017JC013380


%       Loose et al., (2014), "A parameter model of gas exchange for the
%       seasonal sea ice zone", Ocean Sci., 10(4), doi:10.5194/os-10-1-2014.
%
% DISCLAIMER:
%    This software is provided "as is" without warranty of any kind.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%addpath('dependencies');

%%% Do a check if we are on Open Water/Fetch limit/ Sea Ice zone
if strcmpi(I_F_W,'F') 
    sic = 0;
    fetch_eff = f; % set fetch and SI as zero
elseif strcmpi(I_F_W,'I')
    sic = f; % set value from f into fetch
elseif strcmpi(I_F_W,'W')
    sic = 0;
    wave_age = f; % set value from into wave age
end

if nargin > 6
    depth_p = varargin{1};
else
    depth_p = 1000;
end

if nargin > 7
    hum = varargin{2};
else
    hum = 273.16./8.3141./satvap(waterT(:))./1e2*.06;
end
if nargin > 8
    %- Set the Mixed-layer depth (m) if it has not been specified.
    mld = varargin{3};
else
    mld = 300;
end
if nargin > 9
    %- Set the salinity if it has not been specified.
    waterS = varargin{4};
else
    waterS = repmat(32,size(waterT));
end
if nargin > 10
    %- Set the ice thickness if not specified
    z_ice = varargin{5};
else
%     z_ice = 3;
z_ice = 8;
% z_ice = max(sic,1) ; % better estimate, thicker ice at higher %
% z_ice = 0.01;
end

%- SIC (as decimal)
% Check if SI still in vector form 
if numel(sic) > 1
    error('Sea ice % no longer accepts vectors')
end
Ai = sic./100;
airT = airT(:);
waterT = waterT(:);
waterS = waterS(:);
Uice = Uice(:);
Uwind = Uwind(:);

%- convert to Kelvin
waterT = waterT+273.16;
airT = airT+273.16;

%--- COEFICIENTS
g = 9.8;
c_d = (1.03E-3 + 0.04E-3*Uwind^1.48)./(Uwind^0.21); % COARE3 .Edson et.al (2013)
A_charnock = 0.015;      %- Charnocks coeficient
vis = 1.787e-6; %- Kinematic viscosity m2/s
Sc = 660;      %- Schmidt Number (unitless)
n = 0.5;        %- Schmidt Number exponent

% z = 0.01;
z0 = 0.01;       %- Roughness length
acp  = 1.005e3; %- Specific heat of air (J/kg/K)
Lh = 2.26e6;    %- Latent heat of melting (J/kg)
rho_a = 1.013e5/287.058./(airT);
rho_w = gsw_rho_t_exact(waterS, waterT-273.16,0);


%% --------------------------------------------------------------------- %%
%----------------- (1) Air Water SHEAR CALCULATION -----------------------%
%-------------------------------------------------------------------------%
tau_aw = (1-.14)*rho_a.*c_d.*(Uwind.^2);
us_a_aw = sqrt(tau_aw./rho_a);
us_w_aw = sqrt(tau_aw./rho_w);
% Sutherland and Melville (2014) microbreaking coeficient for epsilon

if strcmpi(I_F_W,'F')
    [MB_coef] = Mb_coef_cal(Uwind,'F',fetch_eff,depth_p);
elseif  (strcmpi (I_F_W,'I') && Ai >= 0.01)
    [MB_coef] = Mb_coef_cal(Uwind,'I',Ai,depth_p);
elseif (strcmpi (I_F_W,'I') && Ai < 0.01)
    [MB_coef] = Mb_coef_cal(Uwind,'W',32,depth_p);
    Ai = 0;
elseif strcmpi (I_F_W,'W')
    [MB_coef] = Mb_coef_cal(Uwind,'W',wave_age);
end

z_st_formula = 'n'; % This Formula is EXPERIMENTAL ONLY, keep the value at 'n'
if z_st_formula == 'y'
    %%% due to computational cost we should pre guess eps
    % eps_v = logspace(-14,0,1E6); Accuracy increase with 1E6
    accuracy_inc = 1E5;
    eps_v = logspace(-14,1,accuracy_inc);
    
    kol_l_v = ((vis.^3)./eps_v).^(1/4);
    % Tunning_coef = (5.736.*1E4); % original
    Tunning_coef = (9.*1E3);
    
    z_0w = A_charnock/g .* us_a_aw.^2 + 0.11 .* vis./us_a_aw; %- based on Smith 88
    % z = Kol_l/z0 * L_problem indepedent
    z_aw_v = Tunning_coef.*((vis.^2/g).^(1/3)) .* kol_l_v./z_0w;
    eps_aw_v = MB_coef'.*(us_w_aw.^3)./0.41./z_aw_v.*(1-Ai);
    
    [~,ind_of_match] = min(abs(eps_aw_v-eps_v));% iteration on eps
    eps_aw = eps_v(ind_of_match); % Found eps match
    eps_aw = eps_aw.*(1-Ai);
    z_aw = z_aw_v(ind_of_match); % Corresponding depth
    z = z_aw;
    z_st_out = z/rho_w;
else
    load ('z_st_ref.mat');
    z_aw = interp1(u_ref,z_st,Uwind,'pchip');
    eps_aw = MB_coef'.*(us_w_aw.^3)./0.41./(z_aw).*(1-Ai);
end
z = z_aw;
%% --------------------------------------------------------------------- %%
%----------------- (2) BUOYANCY FLUX CALCULATION -------------------------%
%-------------------------------------------------------------------------%

%- heat capacity and expansion/contraction coefs.
wcp = gsw_cp_t_exact(waterS,waterT-273.16,0);
alpha = gsw_alpha_wrt_t_exact(waterS,waterT-273.16,0);
beta = gsw_beta_const_t_exact(waterS,waterT-273.16,0);
sens_hf = rho_a*acp*c_d.*abs(Uwind).*(waterT-airT);
Jb_sens = 9.8*alpha.*sens_hf/rho_w/wcp;

%- precalculate u*0 because we need for ice-water
us = Uice.*.41/log(z_ice/z0);

%--- Andreas and Murphy (1986), "Bulk transfer coefficients for heat and
%- momentum over leads and polynyas". JPO, V. 16, 1875-1883.
%- Qs = Qsat(ts).  Then Qr is "ref" specific humidity at same height as
%- Uwind.
%- Units are kg-w/m3-air
%- pv = nrt --- n/v(mol/m3)=rt/p.  n/v*.06kg/mol
Qs = 273.16./8.3141./satvap(waterT-273.15)./1e2*.06;

%- Buoyancy flux due to bulk evaporative Latent heat flux.
JqLf = Lh.*c_d.*abs(Uwind).*(Qs-hum).*rho_a;
Jb_Lf = g.*alpha.*JqLf./rho_a./acp;

%- Buoyancy flux due to salinification from evaporation
Jb_salt = g.*beta.*waterS.*(JqLf./Lh)./1e3;


%- Buoyancy flux from ice melt/freeze;
%- Based on interface fluxes using the 3 equation solution (IOBL, 6.9).
%- <w'T'>_0 = w_0*Q_L + q
%- <w'S'>_0 = w_0*(Sice-S0)
%-
%- Reference: McPhee (2010) Air-ice-Ocean Interaction, Springer.
a.us0 = us;
a.alpha_h = .0057;
a.alpha_S = 4e-4;
a.Tml = waterT-273.16;
a.Sml = waterS;
a.Sice = 6;     %- ice salinity
a.wperc = 1e-7;

%- calculate ice temperature as half way between water and air temp.
Tice = nanmean([gsw_t_freezing(waterS,0)+273.16 airT],2);
%Tice = nanmean([gsw_t_freezing(waterS,0)+273.16 airT],2);
%- Kice ~ Kfresh + beta*Sice/Tice
Kice = 2.04 + 0.117*a.Sice./Tice;

a.Hup_ice = -Kice.*(Tice-waterT);

b = S_solve(a);

%- Output b.wb0 is already in W/kg or m2/s3.
Jb_ice = b.wb0;


%- Positive Flux out of the ocean, corresponds to a loss of buoyancy
%- This is why first terms have "negative" signs in front.
Jb = nansum([-(1-Ai).*Jb_sens,-(1-Ai).*Jb_Lf, (1-Ai).*Jb_salt, +Ai.*Jb_ice],2);


%% --------------------------------------------------------------------- %%
%----------------- (3) ICE SHEAR CALCULATION ---------------------------------%
%-------------------------------------------------------------------------%

%- CALC tau_skin-iw the ice/water skin friction using the Rossby similarity
%- theorem:  Vo/u* = 1/k(log(u*/fzo) - A + iB)
%--------------------------------------------
%- Reference: McPhee (2010) Air-ice-Ocean Interaction, Springer.

% now set roughness to 1 cm, Coriolis parameter for 75N
kappa = .41;
f=1.4e-4;
Rostar=b.us0./f./z0;

% compute the Obukhov length
%L=b.us0.^3./kappa./b.wb0;
L=b.us0.^3./kappa./Jb;
Lo = L;
L(L<1) = 1;
L = 1;

% mu* is ratio of planetary to Obukhov scales
mustar=b.us0./f./L;

%- When mustar, stability parameter is negative, stratification is not
%- affecting surface shear.
mustar(mustar<0)=0;

[Und,etastar,A,B]=U0(Rostar,mustar);

Vstable=b.us0.*Und./etastar;
Vsr = real(Vstable).*kappa./(log(Rostar)-A);
Vsi = imag(Vstable).*kappa./(log(Rostar)-A-B);
us_iw = sqrt(Vsr.^2+Vsi.^2);
tau_iw = us_iw.^2.*rho_w;


%- CALC FORM DRAG FROM ESTIMATE OF SEA ICE COVER AND DRIFT VELOCITY -%
%-- This takes into account the recognition that the flow size
%- distribution is a pmf w.r.t. sea ice concentration.

%- Calculate L from SIC:
[floe_dim, dia, PMF] = floe_dimension(sic);
floe_dim = floe_dim.*1e3;


%- Form drag
tau_f =  0.5*rho_w(:).*(z_ice./floe_dim(:)).*Uice.^2;
%- Friction velocity from form drag
us_f = sqrt(tau_f./rho_w);

% REPLACE NANs in US_IW with 0
us_iw(isnan(us_iw)) = 0;


%- FINAL SHEAR-BASED CALCULATIONS ---%
% us_tot = sqrt(Ai.*(us_f.^2 + us_iw.^2) + ((MB_coef'.^(1/3)).*us_aw).^2.*(1-Ai));

%- make a logical index for Obukhov convection criteria
convect = mld(:)./z./Lo;
use = convect<1;

% eps_f = us_f.^3./0.41/z_ice .*Ai; % def
eps_f = us_f.^3./0.41/z_ice .*Ai; % test
% eps_f = 0 

eps_iw = us_iw.^3/0.41/z_ice .*Ai;
    
eps = eps_f + eps_iw + eps_aw;

%- Calculate Epsilon from heat flux and shear.
%- Reference:
eps_Jb = abs(.58*.84*Jb.*use);
% eps_sh = .84*1.76*eps;
% eps_jb = 0 ; % no more convection staright.
eps_sh = eps; % only its effects on shear
%- calculate k from shear (visc. (m2/s), eps (m2/s3)=> m/s
%- multiply by 86400 to get m/d.
k = .419*Sc^(-n).*(vis.*(nansum([eps_Jb,eps_sh],2))).^.25*86400;
keff = (1-Ai).*k;

%%%%DEBUG
% sprintf('eps = %.3f eps_f =%0.3f eps_iw = %.3f  eps_aw = %.3f',eps.*1E6,eps_f.*1E6,eps_iw.*1E6,eps_aw.*1E6)

end



