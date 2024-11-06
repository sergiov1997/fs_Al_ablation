function [R,alpha] = reflectivity_metals(Z,lambda,T_e,T_l)
%This function calcultates the reflectivity of various metals as a function of wavelength and e&l temperatures,
% Sergio Vela Liñán, december 2021

if Z==29  %Atomic number of copper

%For copper the paper Optical properties and thermal response of copper films induced by ultrashort-pulsed
%lasers by Yunpeng Ren, J. K. Chen, and Yuwen Zhang. Eqn (4) of the paper is used    
    
    
c=3e8; %speed of light in m/s
w=(2*pi*c/(lambda*(1e-9)));   %angular frecuency of the laser in rad/s

%We will first calculate each term of the Right Hand Side of equation (1)
%of the paper


epsilon_inf=3.686; %first term of eqn 1

w_D=1.34e16; %in rad/s
B_1=0.562; B_2=27.36; B_3=0.242;
Omega_1=3.205e15; Omega_2=3.43e15; Omega_3=7.33e15; %in rad/s
Gamma_1=0.404e15; Gamma_2=0.77e16; Gamma_3=1.12e15; %in rad/s
Phi_1=-8.185; Phi_2=0.226; Phi_3=-0.516; %Parameters as defined by table 1 of the paper

FIRST_RIGHT_SUM=B_1*Omega_1*((exp(1i*Phi_1)./(Omega_1-w-1i*Gamma_1))+(exp(-1i*Phi_1)./(Omega_1+w+1i*Gamma_1)));
SECOND_RIGHT_SUM=B_2*Omega_2*((exp(1i*Phi_2)./(Omega_2-w-1i*Gamma_2))+(exp(-1i*Phi_2)./(Omega_2+w+1i*Gamma_2)));
THIRD_RIGHT_SUM=B_3*Omega_3*((exp(1i*Phi_3)./(Omega_3-w-1i*Gamma_3))+(exp(-1i*Phi_3)./(Omega_3+w+1i*Gamma_3)));

RIGHT_TERM=FIRST_RIGHT_SUM+SECOND_RIGHT_SUM+THIRD_RIGHT_SUM;  %Last term of eqn (1)


%We now need to find the electron relaxation time ---> eqn (2)

A_e=0.12e6;    %This value, in K^(-2)s^(-1), is given in the paper
B_l=3.35e11;   %This value, in K^(-1)*s^(-1), was found from the fact that the electron relaxation time is required to be 10 fs at room temp. (298K)
Tau_e=(B_l*T_l+A_e*T_e.^2).^-1;  %Electron relaxation time, as given in equation (2)

gamma=1./Tau_e; %The damping coefficient from equation (1) equals the reciprocal of the electron relaxation time.

MIDDLE_TERM=w_D.^2./(w.^2+1i*gamma*w); %Central term in eqn(1)

EPSILON=epsilon_inf-MIDDLE_TERM+RIGHT_TERM; %Equation 1

epsilon_1=real(EPSILON);
%epsilon_2=imag(EPSILON);  %Real and imaginary part of epsilon need to be splitted to find the reflectivity

%We now calculate the functions f_1 and f_2, according to equation (5) of
%the paper

mod_EPS=abs(EPSILON);

f_1=sqrt((epsilon_1+mod_EPS)/2);
f_2=sqrt((-epsilon_1+mod_EPS)/2);


%Finally, we find the reflectivity as defined in equation (4) of the paper

R=((f_1-1).^2+f_2.^2)./((f_1+1).^2+f_2.^2);


%We also find the absorption coefficient in m^-1 to be used later, as
%defined in the same equation.

alpha=(2/c)*w*f_2;

elseif Z==79 %atomic number of gold

%For Gold we use the paper Comparison of gold and silver dispersion laws suitable for FDTD
%simulations by A. Vial and T. Laroche. We will apply eqn(2) of the
%paper with the reflectivity and abosorption coefficient as defined as a
%function of the complex electric permitivity. That is, extended Lorentz-Drude
%model


c=3e8; %speed of light in m/s
w=(2*pi*c/(lambda*(1e-9)));   %angular frecuency of the laser in rad/s

%We add a time-dependent relaxation time, got from the paper Picosecond-laser-pulse-induced heat and mass transfer
%by Vadim Kostrykin et al

Ae=1.2e7; %in K^(-2)*s^(-1)
Bl=1.23*10^11; %in K^(-1)*s^(-1), both are defined in the paper

tau=(Ae*T_e.^2+Bl*T_l).^(-1); %Formula for relaxation time
gamma=tau.^(-1);            %The damping coefficient is the inverse of the relaxation tim

%Now we define the required parameters;

epsilon_inf=1.1431; %first term of eqn (2)

w_D=1.3202e16; %in rad/s
A_1=0.26698; A_2=3.0834; 
Omega_1=3.874e15; Omega_2=4.1684e15;  %in rad/s
Gamma_1=4.4642e14; Gamma_2=2.355e15;  %in rad/s
Phi_1=-1.2371; Phi_2=-1.0968; 

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




elseif Z==13 %atomic number of aluminum

%For Aluminium we also use the paper Comparison of gold and silver dispersion laws suitable for FDTD
%simulations by A. Vial and T. Laroche. We will apply eqn(2) of the
%paper with the reflectivity and abosorption coefficient (table 2) as defined as a
%function of the complex electric permitivity. That is, extended Lorentz-Drude
%model


c=3e8; %speed of light in m/s
w=(2*pi*c/(lambda*(1e-9)));   %angular frecuency of the laser in rad/s

%We add a time-dependent relaxation time, got from the paper Picosecond-laser-pulse-induced heat and mass transfer
%by Vadim Kostrykin et al

Ae=2e7; %in K^(-2)*s^(-1)
Bl=7*10^11; %in K^(-1)*s^(-1), both are defined in the paper

tau=(Ae*T_e.^2+Bl*T_l).^(-1); %Formula for relaxation time
gamma=tau.^(-1);            %The damping coefficient is the inverse of the relaxation tim

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

elseif Z==47 %atomic number of silver
    
%%% For silver we combine the papers Short-pulse ablation rates and the two-temperature model
%%% by B.H. Christensen et al. and Fitting the optical constants of gold,
%%% silver, chromium, titanium, and aluminum in the visible bandwidth
%%% by Barchiesi and Grosges. The former is used to obtain the
%%% temperature dependent free ele electron relaxation time while the later
%%% is used for the rest of the parameters of the critical point model. 

c=3e8; %speed of light in m/s
w=(2*pi*c/(lambda*(1e-9)));   %angular frecuency of the laser in rad/s

%We add a time-dependent relaxation time, got from the paper Picosecond-laser-pulse-induced heat and mass transfer
%by Vadim Kostrykin et al

Ae=0.932e7; %in K^(-2)*s^(-1)
Bl=1.02*10^11; %in K^(-1)*s^(-1), both are defined in the paper

tau=(Ae*T_e.^2+Bl*T_l).^(-1); %Formula for relaxation time
gamma=tau.^(-1);            %The damping coefficient is the inverse of the relaxation tim

%Now we define the required parameters;

epsilon_inf=-19.44; %first term of eqn (2)

w_D=1.325e16; %in rad/s
B_1=1727.44; B_2=-742.06; 
Omega_1=3.475e16; Omega_2=4.761e16;  %in rad/s
Gamma_1=5.183e16; Gamma_2=3.232e16;  %in rad/s
Phi_1=1.042; Phi_2=1.754; 

MIDDLE_TERM=w_D.^2./(w.^2+1i*gamma*w);   %Second term of right hand side of eqn (2)

FIRST_RIGHT_SUM=B_1*Omega_1*  ((exp(1i*Phi_1)./(Omega_1-w-1i*Gamma_1))+(exp(-1i*Phi_1)./(Omega_1+w+1i*Gamma_1)));
SECOND_RIGHT_SUM=B_2*Omega_2*((exp(1i*Phi_2)./(Omega_2-w-1i*Gamma_2))+(exp(-1i*Phi_2)./(Omega_2+w+1i*Gamma_2))); %terms as defined in the sum of eqn (2)


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

end