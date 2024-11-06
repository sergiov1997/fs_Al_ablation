function  [R,alpha]=reflectivity_Al_ours_w_RK_freq(lambda,T_e,T_l)
 %%% this funciton will calculate the reflectivity and absorption
 %%% coefficient of Al from our preexistint DCP model but considering
 %%% Wintter and Redka colision frequency 

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


gamma=tot_sol_col_freq (1e-3*T_e, 1e-3*T_l, 1); % colision frequency according to Winter-Redka model for stdr density x=1;

c=3e8; %speed of light in m/s
w=(2*pi*c/(lambda*(1e-9)));   %angular frecuency of the laser in rad/s

%Now we define the required parameters;

epsilon_inf=1; %first term of eqn (2)

w_D=2.059e16; %in rad/s
A_1=5.2306; A_2=5.2704; 
Omega_1=2.2694e15; Omega_2=2.4668e15;  %in rad/s
Gamma_1=3.2867e14; Gamma_2=1.773e15;  %in rad/s
Phi_1=-0.51202; Phi_2=0.42503; 

MIDDLE_TERM=w_D.^2./(w.^2+1i*gamma*w);   %Second term of right hand side of eqn (2)

FIRST_RIGHT_SUM=A_1*Omega_1*((exp(1i*Phi_1)./(Omega_1-w-1i*Gamma_1))+(exp(-1i*Phi_1)./(Omega_1+w+1i*Gamma_1)));
SECOND_RIGHT_SUM=A_2*Omega_2*((exp(1i*Phi_2)./(Omega_2-w-1i*Gamma_2))+(exp(-1i*Phi_2)./(Omega_2+w+1i*Gamma_2))); %terms as defined in the sum of eqn (2)


RIGHT_TERM=FIRST_RIGHT_SUM+SECOND_RIGHT_SUM;  %Last term of eqn (2)

EPSILON=epsilon_inf-MIDDLE_TERM+RIGHT_TERM; %Finally, we apply eqn (2)
epsilon_1=real(EPSILON);
%epsilon_2=imag(EPSILON);  %Real and imaginary part of the permitivity
mod_EPS=abs(EPSILON);

nn=sqrt(EPSILON); %complex refractive index

%n=sqrt((epsilon_1+mod_EPS)/2); %real physical refractive index
k=sqrt((-epsilon_1+mod_EPS)/2);  %absoprtion inex

R=(abs((nn-1)/(nn+1)))^2; %Reflectivity as defined as a function of the complex refraxctive index

alpha=(2/c)*k*w;


end






