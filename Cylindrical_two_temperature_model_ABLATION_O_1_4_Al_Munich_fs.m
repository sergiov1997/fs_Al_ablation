function [Tl,Te,SHAPE,SAVED_TIMES]=Cylindrical_two_temperature_model_ABLATION_O_1_4_Al_Munich_fs(dt,t_max,z_min,z_max,dz,MAX_R,dr,INITIAL_Te,INITIAL_Tl,I0,w0,tau,T_ab,rho_omega_vap,lamda)


%%% This function simulates two temperature model with ablation in
%%% cylindrical-coordinates.
%%% INPUT:
%%% dt,t_max---> time step, total simulation time
%%% z_min,z_max,dz---> array corresponding to z-coordinate in cylindrical.
%%% MAX_R,dr---> array corresponding to r-coordinate (goes from -MAX_R to
%%% MAX_R)
%%% INITIAL_Te,INITIAL_Tl---> Initial electron and lattice temperatures
%%% I0---> Peak irradiance
%%% tau---> pulse duration
%%% ke---> electron thermal conductivity
%%% abso--->absorptivity
%%% g---> electron-lattice coupling parameter
%%% w0---> beam radius (defined with (1/e)^2)
%%% (gamma--->electronic heat capacity constant)---> We made it 
%%% dependant on lattice temperature.
%%% Cl---> lattice heat capacity (pero unit volume)
%%% rho_omega_vap---> density x mass specific enthalpy of vaporisation
%%% OUTPUT:
%%% Tl,Te---> Fields caracterizing electron and lattice temperatures as a
%%% function of position and time
%%% SAVED_TIMES---> times corresponding to Te and Tl arrays
%%% SHAPE: (3,0) tensor containing the shape of the piece (first and second
%%% indexes) against time (third index)
%%% Sergio Vela Linan, March 2022

disp('Function is now running')


E_p=1.67206088*w0^2*tau*I0;  %Energy per pulse
Z_atom=13; %Atomic number of Al, material to be used

% 
% NNN=8.445263e28; %Electron number density, in m^{-3}
% k_B=1.38e-23;  %Boltzmann constant, in J/K
% T_F=13.6e4; %Fermi temperature, in K


h0=5e13;     %This coefficient will be used for ablation normal velocity
N_integration_points_Cl=1000;  %number of points to integrate numerically when caluclating Tl
[~,Cl_lim,T_Debye]=lattice_specific_heat(Z_atom,N_integration_points_Cl,500);

t=0:dt:t_max;  nt=length(t);  %time array
z=z_min:dz:z_max; nz=length(z); %z-array
height=max(z);  %height of the piece
r=-MAX_R:dr:MAX_R; nr=length(r);  %r-array
SUM=0; COUNTER=0; Nt=1;
dzz=dz*ones(1,length(r));    %as we will have ablation, the size of the z element needs to depend on r.

%%%For speed and memory saving, we will not save the data at all times, we
%%%will only save at multiples of tau/10 (or tau once the pulse has passed).

MODL=round((tau/10)/dt);
BIG_MODL=round(tau/dt);

[~,ZZZ]=meshgrid(r,z);  %This gives the shape of the piece at the beginning.

%Nt will be the number of steps in which the corresponding data is saved.

for I=1:nt
if t(I)<100*tau
if mod(I,MODL)==0
    Nt=Nt+1;
end
elseif mod(I,BIG_MODL)==0
    Nt=Nt+1;
end
end


function g_0 = electron_phonon_coupling_cte_rho_Al(Te)
%%% for the electron phonon coupling factor of Al we will use the model
%%% given in the supplementary material of the paper 
%%% Winter, Jan & Redka, David & Minar, Jan & Schmidt, Michael & Huber, Heinz. (2023). 
%%% Resolving transient temperature and density during ultrafast laser ablation of aluminum. 
%%% Applied Physics A. 129. 10.1007/s00339-023-06922-5. 


factor=(1.3806e-23)/(1.6022e-19); 
T_eV=factor*Te; % converts Te in K to eV

%%%We now compute g0 in (W/m^3*K)
g_0=3.5e17+1e16*T_eV^3/(1+T_eV^0.5+0.9*T_eV^1.9+0.0224*T_eV^3.7); %function given in the paper for g0
end




% function Ce_analytical = electron_heat_capacity_comp_Ploylog_Al(Te)
% %%% This function computes the full-analytocal solution for electronic 
% %%% heat capacity of Añ we obtained in our paper
% %%% NOTE: Symbolic Math Toolbox needs to be installed to compute the
% %%% polylogarithm 
% %%% For Al, to be used for temperatures higher than 50 kK
% 
% 
% 
%  ne=18.1e28; %Atom density of Al, in m^{-3}
% 
% h=6.63e-34; hbar=h/(2*pi); me=9.10938e-31; k_B=1.38065e-23;  %phisical constants, in SI units
% E_F=(3*pi^2*ne)^(2/3)*hbar^2/(2*me);   %Fermi energy, in J
% 
% T_00=(4/(3*sqrt(pi)))^(2/3)*E_F/k_B;   %%%Scaling temperature T_0, as defined in our paper
% 
% Ce_analytical=-3*pi^1.5*(2*me/h^2)^1.5*k_B^2.5*(2.5*Te^1.5*polylog(5/2,inverting_polylogarithm_3_2(-(T_00/Te)^1.5))-1.5*Te^2.5*(T_00^3/(Te)^4)/polylog(1/2,inverting_polylogarithm_3_2(-(T_00/Te)^1.5)));
%  
% 
% end
% 
% function Li_inc = incomplete_polylog(s,b,z)
% 
% %%% This function calculates the incomplete polylogarithm of order s
% 
% %%% Input ---> s:order of the polylogarithm
% %%%            b: first argument, ie, order limit
% %%%            z:argument of the polylogarithm
% 
% fun = @(x,c) x.^c./(exp(x)/z-1);
% 
% 
% Li_inc  = (1/gamma(1+s-1))*  integral(@(x) fun(x,s-1),b,Inf,'RelTol',1e-12);
% 
% end
% 
% 
% function inv_Li = inverting_polylogarithm_3_2(Y)
% 
% %%% This function calculates the inverse of the polylogarithm of order 3/2
% %%% using Newton-Rapson method 
% 
% myfun = @(x,K0) incomplete_polylog(3/2,0,x)-K0;  % parameterized function
% K0 = Y;                    % parameter
% fun = @(x) myfun(x,K0);    % function of x alone
% 
% %%% We want to find the inverse of the polylogarithm, in other words, for
% %%% each given Y we want to find the solution of the equation Li_{3/2}
% %%% (x)=Y
% 
% if Y>-7.889   %%% We use a table of values to find the startinf points of the iterative methods
%     
%     load aaa.mat aaa
%     
%     [~,N]= min(abs(aaa-Y));
%     
%     START_POINT=-100+(N-1)*0.1;
% else
%     START_POINT=-exp((-3*sqrt(pi)*Y/4)^(2/3));  %%% For Y<< -1, the inverse polilogarithm approaxes this function [1]
% end 
% 
% inv_Li = fzero(fun,START_POINT);  %%This finds the zero of the function Li_{3/2} (x) - Y=0.
% 
% end
% 

function Ce_0f_Te = electron_heat_capacity_interpolation(Te)

 load data_Ce_Al.mat comp_T_Ce_full

 Te_vec=comp_T_Ce_full(1,:);  Ce_vec=comp_T_Ce_full(2,:); 

 cent_Te=100*floor(Te/100);

 vect_subs=Te_vec-cent_Te;
 log_INDEX=find(vect_subs==0);
 INDEX=log_INDEX; modl=mod(floor(Te),100);

 
 Ce_0f_Te=(Ce_vec(INDEX)*(100-modl)+Ce_vec(INDEX+1)*modl)/100;


end




function [Cl,lim,TD]=lattice_specific_heat(Z,N_integration_points,Tl)


%%% This function calculates the lattice-heat capacity according to the paper Two-temperature microscale 
%%% heat transfer model part II: determination of lattice parameters by Ewa
%%% Majchrzak and Jolanta Poteralska. Eqn (19) of the paper is used.

kB=1.38e-23; %Boltzmann constan, in J/K

if Z==29  %atomic number of copper

    na=8.4886e28;  %density of atoms in Cu, in m^{-3}
    TD=340; %Debye temperature of Cu, in K

    x_min=1e-8;
    x_max=TD/Tl;  % Integration vatiable, as it appears in eqn (19)
    dx=(x_max-x_min)/N_integration_points;
    x=x_min:dx:x_max;   

    our_function=x.^4.*exp(x)./(exp(x)-1).^2; %Function to be integrated

    int=integral_regla_del_trapecio(x,our_function);

    Cl=9*na*kB*x_max^(-3)*int;   %Lattice heat capacity, according to eqn (19)

    %disp(['El valor de la integral es ' num2str(int)])
    
    lim=3*na*kB;


elseif Z==47

    na=5.85e28;  %density of atoms in Ag, in m^{-3}
    TD=225; %Debye temperature of Ag, in K

    x_min=1e-8;
    x_max=TD/Tl;   % Integration vatiable, as it appears in eqn (19)
    dx=(x_max-x_min)/N_integration_points;
    x=x_min:dx:x_max;   

    our_function=x.^4.*exp(x)./(exp(x)-1).^2; %Function to be integrated

    int=integral_regla_del_trapecio(x,our_function);

    Cl=9*na*kB*x_max^(-3)*int;   %Lattice heat capacity, according to eqn (19)

    %disp(['El valor de la integral es ' num2str(int)])
    
    lim=3*na*kB;




elseif Z==79 %Atomic number of Gold

    na=5.90e28;  %density of atoms in Ag, in m^{-3}
    TD=165; %Debye temperature of Ag, in K

    x_min=1e-8;
    x_max=TD/Tl;  % Integration vatiable, as it appears in eqn (19)
    dx=(x_max-x_min)/N_integration_points;
    x=x_min:dx:x_max;   

    our_function=x.^4.*exp(x)./(exp(x)-1).^2; %Function to be integrated

    int=integral_regla_del_trapecio(x,our_function);

    Cl=9*na*kB*x_max^(-3)*int;   %Lattice heat capacity, according to eqn (19)

    %disp(['El valor de la integral es ' num2str(int)])
    
        lim=3*na*kB;






elseif Z==26 %Atomic number of Fe

    na=17e28;  %density of atoms in Fe, in m^{-3}
    TD=470; %Debye temperature of Fe, in K

    x_min=1e-8;
    x_max=TD/Tl;  % Integration vatiable, as it appears in eqn (19)
    dx=(x_max-x_min)/N_integration_points;
    x=x_min:dx:x_max;   

    our_function=x.^4.*exp(x)./(exp(x)-1).^2; %Function to be integrated

    int=integral_regla_del_trapecio(x,our_function);

    Cl=9*na*kB*x_max^(-3)*int;   %Lattice heat capacity, according to eqn (19)

    %disp(['El valor de la integral es ' num2str(int)])

    
        lim=3*na*kB;


elseif Z==13 %Atomic number of Al

    na=18.1e28;  %density of atoms in Al, in m^{-3}
    TD=428; %Debye temperature of Al, in K

    x_min=1e-8;
    x_max=TD/Tl;  % Integration vatiable, as it appears in eqn (19)
    dx=(x_max-x_min)/N_integration_points;
    x=x_min:dx:x_max;   

    our_function=x.^4.*exp(x)./(exp(x)-1).^2; %Function to be integrated

    int=integral_regla_del_trapecio(x,our_function);

    Cl=9*na*kB*x_max^(-3)*int;   %Lattice heat capacity, according to eqn (19)

    %disp(['El valor de la integral es ' num2str(int)])
    
    
        lim=3*na*kB;






end
end

%%% To perform the integral, we use the trapecium rule

function I=integral_regla_del_trapecio(x,y)

% Given a function y as a vector of values at points x, returns its 
% definite integral between the first point and each of the other points. 
% Uses a linear interpolation approximation for y(x) between points.
% Input:
%   x    : vector of x points
%   y    : vector of y values at given x points
% Output:
%   I   : value of integral
% Usage example:
%    x = 0:pi/100:2*pi
%    y = sin(x)
%    inty = integral(x,y)


n=length(x);
I=0;
for i=2:n
    I = I + (y(i)+y(i-1))/2 * (x(i)-x(i-1));
end

end 



SAVED_TIMES=zeros(1,Nt);


COUN=1; SAVED_TIMES(1)=t(1);

%%% We now create an array that simply tell us the actual times
%%% corresponding to the data which is saved.

for I=1:nt
if t(I)<100*tau
if mod(I,MODL)==0
    COUN=COUN+1;
    SAVED_TIMES(COUN)=t(I);
end
elseif mod(I,BIG_MODL)==0
    COUN=COUN+1;
    SAVED_TIMES(COUN)=t(I);
end
end

clear COUN;
clear I;

COUN=1;

%%%TTe and TTl are the matrices containing electron and lattice
%%% temperatures during the simulation, at the desired times data will be
%%% transfered from the to Te and Tl, which is an output of the function. 
%%% Same process occurs wih DYN_SHAPE and SHAPE, respectively.

Tl=zeros(nz,nr,Nt); Te=Tl;
TTl=zeros(nz,nr,2); TTe=zeros(nz,nr,2);
SHAPE=zeros(nz,nr,Nt); DYN_SHAPE=zeros(nz,nr,2);
SHAPE(:,:,1)=ZZZ;  DYN_SHAPE(:,:,1)=ZZZ;

%%% TTe and TTl only have two time entries, corresponding to he time at
%%% which we are in the simulation and the previous one

Te(:,:,1)=INITIAL_Te;
Tl(:,:,1)=INITIAL_Tl;

TTe(:,:,1)=INITIAL_Te;
TTl(:,:,1)=INITIAL_Tl; %Initial electron and lattice temperatures

for J=1:nt
for I=3:nr-2
for K=3:nz-2
    
    g=electron_phonon_coupling_cte_rho_Al(TTe(K,I,1));

    if TTl(K,I,1)< 2*T_Debye  %Lattice heat capatity is calculated as a function of lattice temp
        Cl=lattice_specific_heat(Z_atom,N_integration_points_Cl,TTl(K,I,1));
    else
        Cl=Cl_lim;  %At temperatures significantly above Debye temp, it can be assumed to be constant
    end


    if TTe(K,I,1)<46e3
    %%% third order taylor expansion we obtained in our paper
    gamma=90.9705; % In J m^{-3} K^{-2}
    beta=-1.467e-8;  % In J m^{-3} K^{-4}
    Ce_of_Te=gamma*TTe(K,I,1)+beta*TTe(K,I,1).^3; 
    elseif TTe(K,I,1)>=46e3 && TTe(K,I,1)<=600e3
        Ce_of_Te=electron_heat_capacity_interpolation(TTe(K,I,1));
    else
    Ce_of_Te=3.613638442183074e+06;  %%% in J/(m^3K); Limit value for extremely high electron temp.
    end

    [ke,~,~]=electron_thermal_conductivity_Al_outside(TTe(K,I,1),TTl(K,I,1));

    [Rf,abso]=reflectivity_Al_ours_w_RK_freq(lamda,TTe(K,I,1),TTl(K,I,1));

    alpha=ke/Ce_of_Te; %Heat transfer coefficient depends on temperature, so it varies during the simulation
    S=abso*I0*(1- Rf)*exp(-(4*log(2)*((t(J)-3*tau)/tau)^2))*exp(-abso*abs(DYN_SHAPE(K,I,1)-DYN_SHAPE(3,I,1)))*exp(-2*r(I)^2/w0^2);
    %%% S is the source term

    %%% We now calculate the second and third term on the right hand side
    %%% of Te PDE equation.

    SECOND_TERM=-g*dt*(TTe(K,I,1)-TTl(K,I,1))/Ce_of_Te;
    THIRD_TERM=S*dt/Ce_of_Te;

    %%% The following coefficients are used to compute the laplacian of Te
    %%% with cylindrical symmetry. 

    COEFFICIENT_i_m_2_k= -(alpha*dt)/(12*dr^2);
    COEFFICIENT_i_k=1-(2.5*alpha*dt)/(dr^2)-(2.5*alpha*dt)/(dzz(I)^2);
    COEFFICIENT_i_M_2_k=-(alpha*dt)/(12*dr^2);
    COEFFICIENT_i_k_m_2=-(alpha*dt)/(12*dzz(I)^2);
    COEFFICIENT_i_k_m_1=(4*alpha*dt)/(3*dzz(I)^2);
    COEFFICIENT_i_k_M_1=(4*alpha*dt)/(3*dzz(I)^2);
    COEFFICIENT_i_k_M_2=-(alpha*dt)/(12*dzz(I)^2);

    COEFFICIENT_i_m_1_k=(4*alpha*dt)/(3*dr^2)-(alpha*dt)/(2*r(I)*dr);
    COEFFICIENT_i_M_1_k=(4*alpha*dt)/(3*dr^2)+(alpha*dt)/(2*r(I)*dr);

    %%% We now perform the actual iterative proocess, according to eqn...

    SUM=SUM+COEFFICIENT_i_m_2_k*TTe(K,I-2,1)+COEFFICIENT_i_m_1_k*TTe(K,I-1,1)+COEFFICIENT_i_k*TTe(K,I,1);
    SUM=SUM+COEFFICIENT_i_M_1_k*TTe(K,I+1,1)+COEFFICIENT_i_M_2_k*TTe(K,I+2,1)+COEFFICIENT_i_k_m_2*TTe(K-2,I,1);
    SUM=SUM+COEFFICIENT_i_k_m_1*TTe(K-1,I,1)+COEFFICIENT_i_k_M_1*TTe(K+1,I,1)+COEFFICIENT_i_k_M_2*TTe(K+2,I,1);

    TTe(K,I,2)=SUM+SECOND_TERM+THIRD_TERM;
    TTl(K,I,2)=TTl(K,I,1)+g*dt*(TTe(K,I,1)-TTl(K,I,1))/Cl;

    SUM=0; %(The only purpose of SUM is to be able to compute the derivative in 3 steps)
    clear alpha;
    clear S;
    clear g;
    clear Ce_of_Te;
    clear ke;
    clear Rf;
    clear abso;







end




  

end

%%% We now apply the boundary condition. This condition is given by:
%%% \div(\grad(T_{e,l}))*\vec{n}=0, where \vec{n} is the normal vector to
%%% each of the four surfaces.

    TTe(2,:,2)=TTe(3,:,2);
    TTe(1,:,2)=TTe(2,:,2);
    TTe(nz-1,:,2)=TTe(nz-2,:,2);
    TTe(nz,:,2)=TTe(nz-1,:,2);

    TTl(2,:,2)=TTl(3,:,2);
    TTl(1,:,2)=TTl(2,:,2);
    TTl(nz-1,:,2)=TTl(nz-2,:,2);
    TTl(nz,:,2)=TTl(nz-1,:,2);


    TTe(:,2,2)=TTe(:,3,2);
    TTe(:,1,2)=TTe(:,2,2);
    TTe(:,nr-1,2)=TTe(:,nr-2,2);
    TTe(:,nr,2)=TTe(:,nr-1,2);

    TTl(:,2,2)=TTl(:,3,2);
    TTl(:,1,2)=TTl(:,2,2);
    TTl(:,nr-1,2)=TTl(:,nr-2,2);
    TTl(:,nr,2)=TTl(:,nr-1,2);

%     DYN_SHAPE(:,1,2)=DYN_SHAPE(:,1,1);
%     DYN_SHAPE(:,2,2)=DYN_SHAPE(:,2,1);
%     DYN_SHAPE(:,nr-1,2)=DYN_SHAPE(:,nr-1,1);
%     DYN_SHAPE(:,nr,2)=DYN_SHAPE(:,nr,1);

%%% We now simulate the ablation process.

   for I=1:nr

   if  TTl(1,I,2)>T_ab %%% ablation criteria according to paper...

        %disp('Hay ablación')

        %%%We create a normal ablation velocity        

        Vn=h0*(TTl(1,I,2)-T_ab )/rho_omega_vap ;

        %disp(['Vn= ' num2str(Vn)]);

        %%% The following 3 vectors are used to interpolate the data once
        %%% ablation is performed

        z_to_interpolate=DYN_SHAPE(:,I,1);

        TTl_to_interpolate=TTl(:,I,2);

        TTe_to_interpolate=TTe(:,I,2);

        %%% We calculate the new depth reached at each r-point.

        new_depth=DYN_SHAPE(1,I,1)+Vn*dt; 



          

        dzz(I)=(height-new_depth)/nz;          %%% dz changes as z-limits change

        %%% We now calculate the new shape of the piece 

        if length(new_depth:dzz:height)==length(DYN_SHAPE(:,I,1))
            DYN_SHAPE(:,I,2)=transpose(new_depth:dzz(I):height);
        else
            dzz(I)=(height-new_depth)/(nz-1);   %%%This is due to divisibility, z can have an extra element
            DYN_SHAPE(:,I,2)=transpose(new_depth:dzz(I):height);
        end

        n_deleted_elements=length(0:dzz(I):Vn*dt);  %Number of elements deleted due to ablation

        %%% We now delete the entries of vectors to interpolate where
        %%% ablation has taken place.

        for III=1:n_deleted_elements

            z_to_interpolate(1)=[];
            TTl_to_interpolate(1)=[];
            TTe_to_interpolate(1)=[];

        end   

        %%% Finally, we create the new temperature fields by linear
        %%% interpolation

        TTl(:,I,2)=interp1(z_to_interpolate,TTl_to_interpolate,DYN_SHAPE(:,I,2)); 
        TTe(:,I,2)=interp1(z_to_interpolate,TTe_to_interpolate,DYN_SHAPE(:,I,2));
        TTl(1,I,2)=TTl(2,I,2);
        TTe(1,I,2)=TTe(2,I,2);

        clear z_to_interpolate;
        clear TTl_to_interpolate;
        clear TTe_to_interpolate;
        clear n_deleted_elements;
        

    


    else    %%%if ablation criteria is not fulfilled, the shape remains the same
        
        DYN_SHAPE(:,I,2)=DYN_SHAPE(:,I,1);

       
    end
        
    end
    

%%% Finaly, we save the data at the desired times. 

if t(J)<100*tau    
if mod(J,MODL)==0
    COUN=COUN+1;
    Te(:,:,COUN)=TTe(:,:,2);
    Tl(:,:,COUN)=TTl(:,:,2);
    SHAPE(:,:,COUN)=DYN_SHAPE(:,:,2); %Wanted data is saved

    TTe(:,:,1)=[];
    TTl(:,:,1)=[];
    DYN_SHAPE(:,:,1)=[]; %Data corresponding to previous time is deleted
else
    TTe(:,:,1)=[];
    TTl(:,:,1)=[];
    DYN_SHAPE(:,:,1)=[]; %Data corresponding to previous time is deleted


end

elseif mod(J,BIG_MODL)==0
    COUN=COUN+1;
    Te(:,:,COUN)=TTe(:,:,2);
    Tl(:,:,COUN)=TTl(:,:,2);
    SHAPE(:,:,COUN)=DYN_SHAPE(:,:,2); %Wanted data is saved

    TTe(:,:,1)=[];
    TTl(:,:,1)=[];
    DYN_SHAPE(:,:,1)=[]; %Data corresponding to previous time is deleted

else
    TTe(:,:,1)=[];
    TTl(:,:,1)=[];
    DYN_SHAPE(:,:,1)=[]; %Data corresponding to previous time is deleted
end

%%% COUNTER is used to show the percentage of the simulation that has been
%%% performed already
disp(['E_p=' num2str(E_p) ' J Progress is of ' num2str(100*COUNTER/nt) ' %; '])
disp([' max electron temp is ' num2str(max (max (TTe(:,:,1)) )) ' K; max lattice temp is ' num2str(max (max (TTl(:,:,1)) )) ' K' ])



COUNTER=COUNTER+1;




end









end