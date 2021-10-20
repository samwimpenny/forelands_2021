function FORCE_BALANCE
%%% Model for calculating the force exterted between mountains
%%% and lowlands given input geometric, thermal and density 
%%% parameters based on balancing the vertical lithostatic stress
%%% differences between the mountains and forelands.

% Written by S. Wimpenny, 2017

% Inputs ---------------------------------
%Variables:
Z_lm_vec=100:1:125      ; %[km], lithospheric thickness (mountains)
Z_lf_vec=100:1:125      ; %[km], lithospheric thickness (lowlands)
Z_cm_vec=65:0.5:75      ; %[km], crustal thickness (mountains)
Z_cf_vec=35:0.5:40      ; %[km], crustal thickness (lowlands)
T_moho_f_vec=600:5:700  ; %[degrees], Moho temperature (lowlands)
T_moho_m_vec=700:5:1000 ; %[degrees], Moho temperature (mountains)
T_m_m_vec=1315          ; %[degrees], Mantle temperature (mountains)
T_m_l_vec=1315          ; %[degrees], Mantle temperature (lowlands)
rho_c0_vec=2800         ; %[kg/m^3], reference crustal density
drho=-50                ; %[kg/m^3] chemical depletion of the mantle lithosphere
rho_a_vec=3330          ; %[kg/m^3], reference asthenosphere density
a_crust_vec=3e-5        ; %[1/K], thermal expansivity (crust)
g=9.81                  ; %[m/s^2], gravitational acceleration

%Model parms:
dz=0.01 ; %[km], model discretisation
it=5e5  ; %number of iterations (use > 10,000).

%-----------------------------------------

%initialising matrices
km2m=1000; %kilometers to metres
dz=dz*km2m; %SI units
Z_lf_vec=Z_lf_vec*km2m;
Z_cm_vec=Z_cm_vec*km2m;
Z_cf_vec=Z_cf_vec*km2m;
Z_lm_vec=Z_lm_vec*km2m;

%%% Setting up variables outside the loop
Z_lm=Z_lm_vec(randi([1 numel(Z_lm_vec)],it,1));
Z_lf=Z_lf_vec(randi([1 numel(Z_lf_vec)],it,1));
Z_cm=Z_cm_vec(randi([1 numel(Z_cm_vec)],it,1));
Z_cf=Z_cf_vec(randi([1 numel(Z_cf_vec)],it,1));
T_moho_f=T_moho_f_vec(randi([1 numel(T_moho_f_vec)],it,1));
T_moho_m=T_moho_m_vec(randi([1 numel(T_moho_m_vec)],it,1));
rho_c0=rho_c0_vec(randi([1 numel(rho_c0_vec)],it,1));
rho_a=rho_a_vec(randi([1 numel(rho_a_vec)],it,1));
a_crust=a_crust_vec(randi([1 numel(a_crust_vec)],it,1));
T_m_m=T_m_m_vec(randi([1 numel(T_m_m_vec)],it,1));
T_m_l=T_m_l_vec(randi([1 numel(T_m_l_vec)],it,1));

%%% MONTLE-CARLO TYPE ITERATION %%%
parfor i=1:it
%Calculations
[F(i),h(i)]=iso_num(Z_lm(i),Z_lf(i),Z_cm(i),Z_cf(i),T_moho_m(i),T_moho_f(i),T_m_m(i),T_m_l(i),rho_c0(i),drho,rho_a(i),a_crust(i),g,dz);
% [F_an(i),h_an(i)]=iso1d(rho_c0(i),rho_a(i)+drho,rho_a(i),Z_cm(i),Z_cf(i),Z_lm(i),Z_lf(i),g);
end

% dstore_an=[F_an' h_an'];
dstore=[F' h'];

fid=fopen('FORCE_BALANCE_PUNA.txt','w');
fprintf(fid,'%5.4f %5.4f \n',dstore');
fclose(fid);

save dstore

% Plotting force against elevation

figure(2)
dat=[dstore(:,2)/1e3 dstore(:,1)];
nbinx=30;
nbiny=30;
n=hist3(dat,[nbinx nbiny]);
n1 = n';
n1(size(n,1) + 1, size(n,2) + 1) = 0;
xb = linspace(min(dat(:,1)),max(dat(:,1)),size(n,1)+1);
yb = linspace(min(dat(:,2)),max(dat(:,2)),size(n,1)+1);
h = pcolor(xb,yb,n1);
xlabel('Elevation, km')
ylabel('Force, N/m')
title('#Models - Isostatic Force Balance')

% Plotting a histogram of models with a given force
% hist(dstore(dstore(:,2)>3e3 & dstore(:,2)<3.5e3,1),40)

%%% ------------------------------------------
%Subfuctions: Geotherm calculatuion & Density conversion, isostatic
%balance.
function [F,h]=iso_num(Z_lm,Z_lf,Z_cm,Z_cf,T_moho_m,T_moho_f,T_m_m,T_m_l,rho_c0,drho,rho_a,a_crust,g,dz)
z = 0:dz:Z_lm;       %depth vector
%%% Computing geotherms
[T_mnt]=linear_geotherm(Z_cm,T_moho_m,Z_lm,T_m_m,z);
[T_lwl]=linear_geotherm(Z_cf,T_moho_f,Z_lf,T_m_l,z);
%%% Computing density profiles
[rho_mnt]=temp2dens_nl(rho_c0,drho,Z_cm,z,T_mnt,a_crust,rho_a,Z_lm);
[rho_lwl]=temp2dens_nl(rho_c0,drho,Z_cf,z,T_lwl,a_crust,rho_a,Z_lf);
%%% Isostatic balance at base of the lithosphere.
P_mnt=trapz(z,g*rho_mnt);
P_lwl=trapz(z,g*rho_lwl);
h=-1*(P_mnt-P_lwl)/(rho_a*g); %assume constant asthenosphere density
rho_lwl_is=zeros(size(rho_lwl));
ind_h=find(min(abs(z-h))==(abs(z-h)));
rho_lwl_is(1:ind_h)=0;
rho_lwl_is(ind_h+1:end)=rho_lwl(1:end-ind_h);
F=trapz(z,(cumtrapz(z,g*rho_mnt)-cumtrapz(z,g*rho_lwl_is))); %Computing force

function [T]=linear_geotherm(Z_moho,T_moho,Z_mantle,T_mantle,z)
%Calculating a linear geotherm given the mantle
%temperature and lithosphere potential temperature
ind_moho=find(min(abs(z-Z_moho))==(z-Z_moho));
ind_lith=find(min(abs(z-Z_mantle))==(z-Z_mantle));
Tcrust=(T_moho/Z_moho)*z(1:ind_moho);
Tmantle=((T_moho-T_mantle)/(Z_moho-Z_mantle))*z(ind_moho+1:ind_lith) + ...
    (T_mantle - (Z_mantle*((T_moho-T_mantle)/(Z_moho-Z_mantle))));
Tasth=ones(1,numel(z(ind_lith+1:end)))*T_mantle;
T=[Tcrust Tmantle Tasth];

function [rho]=temp2dens_nl(rho_c0,drho,Z_moho,z,T,a_crust,rho_a,Z_mantle)
%Calculating the density profile based on the following:
%rho = rho_0(1-aT) for crust and lithospheric mantle and
%a linear dependence of expansivity on mantle temperature.
rho=zeros(size(z));
m=9.3e-9 ; c=2.993e-5; %parameters from Bouhifd [1996] assuming linear fit to data
ind_moho=find(min(abs(z-Z_moho))==(z-Z_moho));     %index of moho 
ind_lith=find(min(abs(z-Z_mantle))==(z-Z_mantle)); %index of lithospheric mantle
rho(1:ind_moho)=rho_c0*(1-a_crust*T(1:ind_moho));  %crustal
rho(ind_moho+1:ind_lith)=(rho_a*(1-(m*T(ind_moho+1:ind_lith)+c).*T(ind_moho+1:ind_lith)))+drho; %LM, adding depletion too.
rho(ind_lith+1:end)=(rho_a*(1-(m*T(ind_lith+1:end)+c).*T(ind_lith+1:end)));

function [F,H]=iso1d(rho_c,rho_l,rho_a,Z_cm,Z_cf,Z_lm,Z_lf,g)
%Finding analytical solution to the 1-d isostatic balance
%excluding the effects of thermal buoyancy.
Z_lm=Z_lm-Z_cm;
Z_lf=Z_lf-Z_cf;
H = (1/rho_a) *( rho_c*(Z_cf-Z_cm) + rho_l*(Z_lf-Z_lm) + rho_a*(Z_cm+Z_lm-Z_cf-Z_lf) ); %pressure balance for elevation
%geometric factors
r1=Z_cm-Z_cf-H;
r2=Z_lf+Z_cf+H-Z_cm;
r3=Z_cm+Z_lm-Z_cf-Z_lf-H;
%force calculations
F1=rho_c*g*0.5*(H.^2);
F2=rho_c*g*H*Z_cf;
F3=rho_c*g*H*r1 + (rho_c-rho_l)*0.5*g*(r1.^2);
F4=rho_c*g*r2*(Z_cm-Z_cf) - rho_l*g*r1*r2;
F5=(rho_l-rho_a)*g*0.5*(r3.^2) + g*r3*( rho_l*r2 + rho_c*Z_cm - rho_l*Z_lf - rho_c*Z_cf);
F=F1+F2+F3+F4+F5;
