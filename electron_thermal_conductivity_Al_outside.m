function [k_tot,kse,ksi] = electron_thermal_conductivity_Al_outside(T_e,T_l)


function Ce_0f_Te = electron_heat_capacity_interpolation(Te)

 load data_Ce_Al.mat comp_T_Ce_full

 Te_vec=comp_T_Ce_full(1,:);  Ce_vec=comp_T_Ce_full(2,:); 

 cent_Te=100*floor(Te/100);

 vect_subs=Te_vec-cent_Te;
 log_INDEX=find(vect_subs==0);
 INDEX=log_INDEX; modl=mod(floor(Te),100);

 
 Ce_0f_Te=(Ce_vec(INDEX)*(100-modl)+Ce_vec(INDEX+1)*modl)/100;


end

%%%For the electron heat capacity, we combine models given in the Winter
%%%paper with our improved model for the electron heat capacity

 ne=18.1e28; %Atom density of Al, in m^{-3}

h=6.63e-34; hbar=h/(2*pi); me=9.10938e-31; k_B=1.38065e-23;  %phisical constants, in SI units
E_F=(3*pi^2*ne)^(2/3)*hbar^2/(2*me);   %Fermi energy, in J
v_F=sqrt(2*E_F/me);  %Fermi velocity, in m/s
ms=1.05*me; %%effective electron mas for s/p electrons in Al

gamma=90.9705; % In J m^{-3} K^{-2}
beta=-1.467e-8;  % In J m^{-3} K^{-4}

if T_e<46000   %%% Cubic approximation will be used only at its validity range, ie <46000 K
   C_ee=gamma*T_e+beta*T_e^3;
elseif T_e>=46e3 && T_e<=600e3
    C_ee=electron_heat_capacity_interpolation(T_e);
else
   C_ee= 3.613638442183074e+06;
end
% 
% kse=(2/3)*(1/me)*hbar*(E_F^2/(k_B^2*T_e^2))*C_ee;

t=6*k_B*T_e/E_F;

b0=-0.91906; b1=0.79564; b2=0.13456; a0=4e-4;

kse=(1+b0*t^(0.5)+b1*t+b2*t^2)/(a0*t);  %as given in the winter paper

krt=237; %%% in W/(mK), electron-ion contribution at room temp.
Trt=298; %room temperature, in K
Crt=gamma*Trt+beta*Trt^3; %electron heat caacity at room temperature
vsrt=sqrt(v_F^2+3*k_B*Trt/ms); %mean velocity of s/p electrons at room temperature
vs=sqrt(v_F^2+3*k_B*T_e/ms); %mean velocity of s/p electrons at Te

ksi=krt*(C_ee*vs/(Crt*vsrt))*(Trt/T_l);

k_tot=(1/kse+1/ksi)^(-1); %%% Electron heat capacity calculated according to matthisen rule

end

