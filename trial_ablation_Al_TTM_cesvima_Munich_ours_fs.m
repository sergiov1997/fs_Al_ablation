%%%This code will simulate fs ablation of Al after using our previously
%%% defined TTM code
clear variables


AVOID_ZERO_FACTOR=1-pi/400;  %This is to avoin a zero entry in r
%ke=237;
lambda_1=1056 ; %Wavelength, in nm
%g_RT=36e16; %Electron-lattice coupling factor, in J m^{-3} K^{-2}
T_ab=6030; %% We changed the ablation temperature to 0_9 times the critical temperature T_c=6700 K
E_p=33.51*10^(-6); %Energy per pulse in J
tau=320e-15; %Pulse duration
w0=12.162e-6; %Beam radius
abso=1e8; %Absorptivity
Nz=260;
Nr=260;    %Number of z and r elements
MAX_R=4*w0; %Maximum r
dr=AVOID_ZERO_FACTOR*(2*MAX_R)/Nr;
r=-MAX_R:dr:MAX_R;
z_min=0;
z_max=3*w0;
dz=AVOID_ZERO_FACTOR*(z_max-z_min)/Nz; z=z_min:dz:z_max; %z and r arrays
I0=E_p*4*sqrt(log(2))/(pi*w0^2*tau*sqrt(pi)); %Peak irradiance
dt=tau/70; t_max=65*tau; t=0:dt:t_max; %TIme array
rho_omega_vap=4.22e10;  
%gamma=70.361;
%Cl=3.5e6;%density x mass specific enthalpy of vaporisation
[R,~]=meshgrid(r,z);

%%% All is in SI untis


Te_START=298*ones(length(z),length(r),1); Tl_START=Te_START; %Initial electron and lattice temp.


[Tl,Te,SHAPE,SAVED_TIMES]=Cylindrical_two_temperature_model_ABLATION_O_1_4_Al_Munich_fs(dt,t_max,z_min,z_max,dz,MAX_R,dr,Te_START,Tl_START,I0,w0,tau,T_ab,rho_omega_vap,lambda_1);

T_STEPS=length(SAVED_TIMES);

%%% We plot results dynamically

% figure(101)
% for i = 1:T_STEPS
% %     if mod(i,5)==0
%   pcolor(R,SHAPE(:,:,i),Te(:,:,i))
%   colorbar;
%   xlabel('$r[m]$','Interpreter','latex','FontSize',20);
%   ylabel('$z[m]$','Interpreter','latex','FontSize',20);
%   set(gca,'FontSize',20)
%   %zlim([297, max(max(T(:,:,i)))])
%   title(['T_e[K];   t= ' num2str(1e15*SAVED_TIMES(i)) ' fs']);
%   %colormap(autumn)
%   shading interp;
%   pause(0.075)
% %     end
% end
% 
% 
% figure(102)
% for i = 1:T_STEPS
% %     if mod(i,5)==0
%   pcolor(R,SHAPE(:,:,i),Tl(:,:,i))
%   colorbar;
%   xlabel('$r[m]$','Interpreter','latex','FontSize',20);
%   ylabel('$z[m]$','Interpreter','latex','FontSize',20);
%   set(gca,'FontSize',20)
%   %zlim([297, max(max(T(:,:,i)))])
%   title(['T_l[K];   t= ' num2str(1e15*SAVED_TIMES(i)) ' fs']);
%   %colormap(gray)
%   shading interp;
%   pause(0.075)
% %     end
% end
% 
% %%% We plot the final curve of the ablation zone
% 
% figure(103)
% plot(r,SHAPE(1,:,T_STEPS),'ko')
% set(gca,'FontSize',25) 
% xlabel('$r[m]$','Interpreter','latex','FontSize',25);
% ylabel('$z[m]$','Interpreter','latex','FontSize',25);
% title('Final shape of the ablation zone','FontSize',25)
%  
%% Uncomment for saving the Data after simulation
save ("data_after_ablation_Al_12_16_mu_m_320_fs_33_51_mu_J_1056_nm_tmax_65_tau_h0_5e13_T_ab_0_9_Tc_Munich_ours_fs.mat","R","r","Tl","Te","SHAPE","SAVED_TIMES","T_STEPS") 



% 
% 

