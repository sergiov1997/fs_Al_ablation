function [R,alpha]=reflectivity_Al_Redka_Winter(Te,Tl,lambda)


function y0 = y0_func(x)

% x = Dichteverhältniss p/p0

%Parameter
a=1.425;
b=1.359;
c_ab=(a-b)/(b+1);

y0=(1+c_ab).*x.^(2*a+1)./(1+c_ab.*x.^(a+1));

end

function t = t_func(T_e,x)

% t=normalized temperature (1)
% T_e = elektronen Temp in kK, x = Dichteverhältniss p/p0, p0=2.735 g/cm^3

k_B = 8.6173303e-5 * 10^3; % (eV/kK)
E_F0 = 11.1; % (eV) Fermi-Energie

t = (6.*k_B.*T_e)./(E_F0.*x.^(2/3));

end

function nu_ee = el_el_col_freq (T_e, x)

% T_e (kK), x (1) = Dichteverhältniss rho/rho0

%Parameter
eta_ee = 0.75;
a0 = 0.49188e15; %(1/s)
a1 = 1.2;
a2 = 0.8;
p0 = 2/3;
p1 = 0.75;
p2 = 2.25;

%  electron-electron collision frequency (1/s)

nu_ee = eta_ee.* x.^p0.*a0.*t_func(T_e,x).^2.*(1 + a1.*t_func(T_e,x).^p1)./(1 + a2.*t_func(T_e, x).^p2);

end

function nu_ei_sol = el_ion_sol_col_freq (T_i, x)

% T_i (kK) Gitter Temp, x (1) = Dichteverhältniss rho/rho0

%Kontsanten
nu_rt_sol = 1.167022e14; %(1/s) nu_sol at room temp
eta_sol = 0.5175;
T_rt = 0.293; % room temperature (kK)
rho_0=2.735; %(g/cm^3)
rho_rt=2.71223; %(g/cm^3)
x_rt=rho_rt/rho_0; %  Dichteverhältniss rho_rt/rho_0, rho__rt = density at room temperature


%  electron-ionen-solid collision frequency (1/s)
nu_ei_sol = nu_rt_sol + eta_sol.* nu_rt_sol.*(x_rt./x).^(2/3).*y0_func(x_rt)./y0_func(x).*T_i./T_rt;

end

function nu_tot_sol = tot_sol_col_freq (T_e, T_i, x)

% T_e (kK) elektronen Temp, T_i (kK) Gitter Temp, x (1) Dichteverhältniss rho/rho0


% total solid collision frequency (1/s)

nu_tot_sol = el_el_col_freq(T_e,x) + el_ion_sol_col_freq(T_i,x);

end

function epsilon = dielectric_function_solid (T_e, T_i, x, omega)

% T_e electron Temp (kK), T_i Gitter Temp (kK), x (1) Dichteverhältniss
% rho/rho_0,  omega (1/s) excitation angular frequency 


% Konstanten Drude Term
omega_p0 = 16.4868e15; % (1/s) Plasma Frequenzy

% Konstanten CPs Terms
A1 = 2.17020801; %(1)
Omega_p1 = 2.25028e15; % (1/s)
Gamma_p1 = 0.268726e15; % (1/s)
Phi_p1 = -0.52279535; %(1)

A2 = 3.95505075; %(1)
Omega_p2 = 2.26699e15; % (1/s)
Gamma_p2 = 1.16111e15; % (1/s)
Phi_p2 = 0.12358051; %(1)


% Drude Term
Drude = 1 - x.*omega_p0.^2 ./ (omega.*(omega + 1i.*tot_sol_col_freq(T_e, T_i, x)));

% CPs Terms
epsilon_1 = A1.*Omega_p1.*(exp(1i*Phi_p1) ./ (Omega_p1-omega-1i.*Gamma_p1) + exp(-1i.*Omega_p1) ./ (Omega_p1+omega+1i.*Gamma_p1));
epsilon_2 = A2.*Omega_p2.*(exp(1i*Phi_p2) ./ (Omega_p2-omega-1i.*Gamma_p2) + exp(-1i.*Omega_p2) ./ (Omega_p2+omega+1i.*Gamma_p2));

% dielectric function 
epsilon = Drude + epsilon_1 + epsilon_2;

end



c=3e8;
omega=(2*pi*c/(lambda*(1e-9)));


eps=dielectric_function_solid (1e-3*Te,1e-3*Tl,1,omega);
nn=sqrt(eps);

R=abs((nn-1)/(nn+1))^2;

alpha=2*omega*imag(nn)/c;

end