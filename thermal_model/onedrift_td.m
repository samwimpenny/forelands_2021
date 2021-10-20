function onedrift_td
%%% Script for solving for the temperature distribution 
%%% within the continental lithosphere due to stretching.
%%% This model incorporates the temperature 
%%% and compositional dependence of lithospheric 
%%% materials and can have arbitrary complexity
%%% but assumes the problem can be linearized and 
%%% solved by iteration of the form: A.x^m+1 = b + d

% Inputs:
zuc0  = 20  ;
zlc0  = 30  ;
zl0   = 125 ;
beta  = 1.5 ;
delta = 1.5 ;
trift = 20  ;
Hs_uc = 1.5e-6 ;
Hs_lc = 0.6e-6 ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Model Setup:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parameters prior to deformation:

zc0  = zuc0 + zlc0; %[km], original crustal thickness
zlm0 = zl0 - zc0  ; %[km], original lithospheric mantle
z_nm = 150        ; %[km], Maximum thickness of post-rift lithosphere.

% Model scales:

tend = 200  ; % [Myrs], end of model
dz0  = 1.0  ; % [km], initial spatial step
dt   = 1e5  ; % [yrs], time step

% Material properties: constants

Ts = 0    ; %[deg], surface temperature
Tp = 1335 ; %[deg], mantle potential temperature
g  = 9.81 ; %[m/s^2], acceleration due to gravity

rhosed   = 2200 ;
rho_o_uc = 2800 ; %[kg/m^3], density at 0 degrees
rho_o_lc = 2900 ; 
rho_o_lm = 3300 ; % depleted lithospheric mantle
rho_o_m  = 3300 ; % undepleted asthenospheric mantle

k_uc = 2.5  ; %[W/m/K], Thermal Conductivity
k_lc = 2.5  ;
Cp_uc = 750 ; %J/K/kg, Specific Heat Capacity
Cp_lc = 750 ;
Hs_m  = 0   ; %W/m^3, radiogenic heat production

% Material constants for temperature-dependent parameters

crad = [0.0175,0.000103,2.245e-7,3.407e-11]; % Thermal conductivity - McKenzie '05
clat = [5.3,0.0015,1.753e-2,-1.0365e-4,2.2451e-7,-3.4071e-11];
alpha = [2.832e-5, 0.758e-8]; % Bouhifd '96 thermal expansivity
c = [1580,12230,1694e6]; % Specific heat capacity - Korenaga and Korenaga '16

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Converting to SI units:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

km2m = 1e3;
yr2sec = 60*60*24*365.25;

zuc0 = zuc0*km2m ;
zlc0 = zlc0*km2m ;
zlm0 = zlm0*km2m ;
zc0  = zc0*km2m  ;
zl0  = zl0*km2m  ;
z_nm = z_nm*km2m ;
dz0   = dz0*km2m;
trift = trift*yr2sec*1e6;
dt   = dt*yr2sec;
tend = tend*yr2sec*1e6;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialising Arrays and Indexes:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Base of the model domain

if (zlm0/delta) + (zc0/beta) < z_nm
    dL = zc0*(1-(1/beta)) + zlm0*(1-(1/delta)) ;
    zbase0 = zc0 + zlm0 + dL ;
elseif (zlm0/delta) + (zc0/beta) >= z_nm
    zbase0 = zl0 ;
end

% Model indexing - Lagrangian form.

z = 0:dz0:zbase0;
ind_uc = find(z>=0 & z<=zuc0);
ind_lc = find(z>zuc0 & z<=zc0);
ind_c  = [ind_uc ind_lc];
ind_lm = find(z>zc0 & z<=zl0);
ind_l  = [ind_uc ind_lc ind_lm];
ind_m  = find(z>zl0);
 
% Populating arrays

k   = zeros(size(z)) ; %Thermal conductivity array
rho = zeros(size(z)) ; %Density array
Cp  = zeros(size(z)) ; %Specific Heat capacity array
H   = zeros(size(z)) ; %Radiogenic heating array

% Temperature-independent properties
% in the lithospheric mantle as starting.

Tc     = 900 ; % average mantle temperature
k_lm   = k_ol(Tc,clat,crad);
Cp_lm  = Cp_ol(Tc,c);
rho_lm = rho_ol(Tc,alpha,rho_o_lm);
rho_m  = rho_ol(Tc,alpha,rho_o_m);
alpha_lm = alpha(1) + alpha(2)*Tc ;
gamma = Tp*g*alpha_lm/Cp_lm; % isentrope gradient

% Filling material property arrays

k(ind_uc) = k_uc; % Thermal conductivity
k(ind_lc) = k_lc; 
k(ind_lm) = k_lm;
k(ind_m)  = k_lm;
rho(ind_uc) = rho_o_uc; % Density
rho(ind_lc) = rho_o_lc;
rho(ind_lm) = rho_lm;
rho(ind_m)  = rho_m;
Cp(ind_uc) = Cp_uc; % Specific Heat Capacity
Cp(ind_lc) = Cp_lc;
Cp(ind_lm) = Cp_lm;
Cp(ind_m)  = Cp_lm;
H(ind_uc) = Hs_uc; % Radiogenic heat production
H(ind_lc) = Hs_lc;
H(ind_lm) = Hs_m;
H(ind_m)  = Hs_m;

% Checking the time-step is smaller than 
% the Courant stability condition.

if dt > min(5*rho.*Cp*(dz0^2)./k)
    disp('Possibly unstable time step')
    min(5*rho.*Cp*(dz0^2)./k) / yr2sec;
else
    disp('Stable time step')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setting up initial condition:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculating steady-state profile

[T,k,rho,Cp,H] = steadystate(z,zuc0,zlc0,zlm0,k,rho,Cp,H,Tp,...
    ind_uc,ind_lc,ind_lm,ind_m,dz0,rho_o_m,gamma,rho_o_lm,Ts,...
    clat,crad,alpha,c);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time-stepping:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialising arrays

t = 0:dt:tend; % time vector
nz = numel(z); % number of spatial nodes
nzc  = numel(z(ind_c));
nzlm = numel(z(ind_lm));
nza  = numel(z(ind_m));
nt = numel(t); % number of time nodes

zarr = zeros(nz,nt); % array of z-values
zarr(:,1) = z;
tarr = repmat(t,nz,1); % array of t-values
temptime = zeros(nz,nt); % array of temperature
temptime(:,1) = T;

zcrust = zeros(nz,nt); % array of crustal thickness
zcrust(:,1) = (zarr(:,1)<=zarr(ind_lc(end),1));

karr = repmat(k',1,nt);
rhoarr = repmat(rho',1,nt);
Cparr = repmat(Cp',1,nt);
d = zeros(nz,1);

% Pressure at the base of lithosphere (*1/g)
% zsub = zeros(1,nt)   ; % array of subsidence
% P1 = trapz(zarr(zarr(ind_l,1)<=zl0,1),rhoarr(zarr(:,1)<=zl0,1));
% Pmoho  = trapz(zarr(ind_c,1),g*rhoarr(ind_c,1));

% Lithosphere strain:

vc = ((1/(beta))-1)*(zc0/trift)   ; % crust thinning rate
vl = ((1/(delta))-1)*(zlm0/trift) ; % lithosphere thinning rate

% Time-step loop
epsilon = 0.5; % Convergence parameters
%for i = 2:nt
for i = 2:nt
        
    % Re-sizing z-vector 
    
    if t(i) <= trift
        
        zc   = zc0 + (vc*t(i))  ;  % thickness of crust
        zlm  = zlm0 + (vl*t(i)) ;  % thickness of lithospheric mantle
        dzlm = zlm/nzlm ; % spatial step in lithospheric mantle
        zb   = zbase0 + vc*t(i) + vl*t(i) ; % base of the model domain
        
        zarr(ind_c,i)  = linspace(0,zc,nzc) ;
        zarr(ind_lm,i) = linspace(zc+dzlm,zc+zlm,nzlm)  ;
        zarr(ind_m,i)  = linspace(zc+zlm+dz0,zb,nza) ;
                
    elseif t(i) > trift
        
        zarr(ind_c,i)  = linspace(0,zc,nzc) ;
        zarr(ind_lm,i) = linspace(zc+dzlm,zc+zlm,nzlm)  ;
        zarr(ind_m,i)  = linspace(zc+zlm+dz0,zb,nza) ;
        
    end
    
    zcrust(:,i) = (zarr(:,i) <= zarr(ind_lc(end),i)) ; % tracking crust

    % Defining where the asthenosphere BC to be applied
    
    if zlm+zc > z_nm && zl0 > z_nm
        ind_bc = find(zarr(:,i)>=zlm+zc); % apply at base of lith
    elseif zlm+zc <= z_nm && zl0 > z_nm
        ind_bc = find(zarr(:,i)>=z_nm); % apply at 150 km post-rift
    elseif zlm+zc <= z_nm && zl0 <= z_nm 
        ind_bc = find(zarr(:,i)>=zl0); % apply at original pre-rift zl0
    end
            
    % Iteratively solving for next time-step
    
    dT = 50; % re-set the dT value
    j  = 1 ; % predictor-corrector iterator
    
    while dT>epsilon && j<=5
    
        % Creating coefficient matrices and
        % solving: A-1 * (B.x + d)
        
        if j == 1; % Using old material properties in first iteration
            m = i - 1;
            T = temptime(:,i);
            % Updating material properties
            karr(ind_lm,i)   = k_ol(temptime(ind_lm,m),clat,crad);
            rhoarr(ind_lm,i) = rho_ol(temptime(ind_lm,m),alpha,rho_o_lm);
            Cparr(ind_lm,i)  = Cp_ol(temptime(ind_lm,m),c);
            karr(ind_m,i)    = k_ol(temptime(ind_m,m),clat,crad);
            rhoarr(ind_m,i)  = rho_ol(temptime(ind_m,m),alpha,rho_o_m);
            Cparr(ind_m,i)   = Cp_ol(temptime(ind_m,m),c);            
        elseif j > 1; % using new material properties in next iteration
            m = i;
            % Updating material properties
            karr(ind_lm,i)   = k_ol(temptime(ind_lm,i),clat,crad);
            rhoarr(ind_lm,i) = rho_ol(temptime(ind_lm,i),alpha,rho_o_lm);
            Cparr(ind_lm,i)  = Cp_ol(temptime(ind_lm,i),c);
            karr(ind_m,i)    = k_ol(temptime(ind_m,i),clat,crad);
            rhoarr(ind_m,i)  = rho_ol(temptime(ind_m,i),alpha,rho_o_m);
            Cparr(ind_m,i)   = Cp_ol(temptime(ind_m,i),c);
        end
        
       % dx vectors
       dxnp1_jp1 = zarr(3:end,i) - zarr(2:end-1,i)  ;
       dxnp1_j   = zarr(2:end-1,i) - zarr(1:end-2,i) ;
       dxn_jp1   = zarr(3:end,i-1) - zarr(2:end-1,i-1) ;  
       dxn_j     = zarr(2:end-1,i-1) - zarr(1:end-2,i-1) ;
       Dnp1      = dxnp1_jp1+dxnp1_j ;
       Dn        = dxn_jp1+dxn_j ;
       
       % k vectors at mid-points
       knp1_jph = 0.5*(karr(3:end,m) + karr(2:end-1,m)) ;
       knp1_jmh = 0.5*(karr(2:end-1,m) + karr(1:end-2,m)) ;
       kn_jph   = 0.5*(karr(3:end,i-1) + karr(2:end-1,i-1)) ;
       kn_jmh   = 0.5*(karr(2:end-1,i-1) + karr(1:end-2,i-1)) ;
       
       % Density and specific heat capacity vectors
       rhonp1_j = rhoarr(2:end-1,m);
       Cpnp1_j  = Cparr(2:end-1,m) ;
       rhon_j   = rhoarr(2:end-1,i-1);
       Cpn_j    = Cparr(2:end-1,i-1);
       R = 0.5*(1/dt)*((rhonp1_j.*Cpnp1_j) + (rhon_j.*Cpn_j)) ;
          
       % A Matrix (L.H.S)
       r1 = [ 0; knp1_jph./dxnp1_jp1];
       r2 = [ 0; -1*(Dnp1.*R) - (knp1_jph./dxnp1_jp1) - (knp1_jmh./dxnp1_j) ; 0];
       r3 = [ knp1_jmh./dxnp1_j; 0];
       A  = maketridiag(r1,r2,r3,ind_bc) ;
       
       % B matrix (R.H.S)
       r1 = [ 0; -1*(Dnp1./Dn).*(kn_jph./dxn_jp1)] ;
       r2 = [ 0; -1*(Dnp1.*R) + (Dnp1./Dn).*((kn_jph./dxn_jp1)+(kn_jmh./dxn_j)) ; 0] ;
       r3 = [ -1*(Dnp1./Dn).*(kn_jmh./dxn_j); 0] ;
       B  = maketridiag(r1,r2,r3,ind_bc) ;    
            
       % d Matrix (source term)
       d(1) = 0;
       d(2:end-1) = -1* H(2:end-1)' .* Dnp1 ;
       d(ind_bc) = 0 ;
       
       % Solving for temperature
       temptime(1,i-1) = Ts; %top boundary condition
       temptime(ind_bc,i-1) = Tp + gamma*zarr(ind_bc,i-1); %bottom boundary condition
       temptime(:,i) = A\(B*temptime(:,i-1) + d) ;

       % Checking solution convergence
       dT = (mean(temptime(:,i)-T)^2)^0.5;
       T = temptime(:,i);

       % Iterating predictor-corrector
       j = j + 1;

    end
 
%     % Calculating the sediment-loaded subsidence
%     indl = zarr(:,i) <= zarr(ind_bc(1),i) ;
%     P2 = trapz(zarr(indl,i),rhoarr(indl,i));
%     rho_asth = rho_ol(Tp,alpha,rho_o_m);
%     zsub(i)  = (1/(rhosed-rho_asth)) * ...
%               (P1 - P2 - rho_asth*(zl0 - zarr(ind_bc(1),i)));
%             
%    % P-T-t at the Moho inc. sediment loading
%    if zsub(i) > 0
%        Pmoho = trapz(zarr(ind_c,i),g*rhoarr(ind_c,i)) + (rhosed*g*zsub(i));
%    elseif zsub(i) <= 0
%        Pmoho = trapz(zarr(ind_c,i),g*rhoarr(ind_c,i)) ;
%    end
%    PTmoho(i,:) = [Pmoho/1e9 temptime(ind_c(end),i) t(i)/(yr2sec*1e6)];
   
   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Saving model runs: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Saving temperature-time array as vector:
tt = [tarr(:)/(yr2sec*1e6) zarr(:)/km2m temptime(:)];
fname = ['output_dT/temp_' num2str(zuc0/km2m,'%2d') '_'  num2str(zlc0/km2m,'%2d')  ... 
         '_' num2str(zl0/km2m,'%3d') '_' num2str(beta,'%1.1f') ...
         '_' num2str(trift/(yr2sec*1e6),'%2d') '_' num2str(Hs_uc*1e6,'%2.1f') ...
         '_' num2str(Hs_lc*1e6) '_' num2str(delta,'%1.1f') '.out'];
disp(fname)
dlmwrite(fname,tt,'delimiter','\t')

% % Saving subsidence history:
% subsidence = [t'/(yr2sec*1e6) zsub'/km2m];
% subname = ['output/subs_' num2str(zuc0/km2m,'%2d') '_'  num2str(zlc0/km2m,'%2d')  ... 
%          '_' num2str(zl0/km2m,'%3d') '_' num2str(beta,'%1.1f') ...
%          '_' num2str(trift/(yr2sec*1e6),'%2d') '_' num2str(Hs_uc*1e6,'%2.1f') ...
%          '_' num2str(Hs_lc*1e6) '_' num2str(delta,'%1.1f') '.out'];
% dlmwrite(subname,subsidence,'delimiter','\t')

% % Saving crustal thickness
% zcr = [t(1,:)'/(yr2sec*1e6) max(zarr.*zcrust)'/km2m];
% crname = ['output/crust_' num2str(zuc0/km2m,'%2d') '_'  num2str(zlc0/km2m,'%2d')  ... 
%          '_' num2str(zl0/km2m,'%3d') '_' num2str(beta,'%1.1f') ...
%          '_' num2str(trift/(yr2sec*1e6),'%2d') '_' num2str(Hs_uc*1e6,'%2.1f') ...
%          '_' num2str(Hs_lc*1e6) '_' num2str(delta,'%1.1f') '.out'];
% dlmwrite(crname,zcr,'delimiter','\t')

% % Saving P-T-t path of the Moho
% ptname = ['output/PTmoho_' num2str(zuc0/km2m,'%2d') '_'  num2str(zlc0/km2m,'%2d')  ... 
%          '_' num2str(zl0/km2m,'%3d') '_' num2str(beta,'%1.1f') ...
%          '_' num2str(trift/(yr2sec*1e6),'%2d') '_' num2str(Hs_uc*1e6,'%2.1f') ...
%          '_' num2str(Hs_lc*1e6) '_' num2str(delta,'%1.1f') '.out'];
% dlmwrite(ptname,PTmoho,'delimiter','\t')

% Saving isentropic temperature with depth:
% isentrope = [z'/km2m Tp+gamma*z'];
% dlmwrite('isentrope.out',isentrope,'delimiter','\t')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Subfunctions: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [coeff] = maketridiag(r1,r2,r3,ind_bc)
%%% Makes a tri-diagonal matrix of the form
%%% needed for 1-D temperature dependent heat
%%% flow equation with all values at ind_m 
%%% set to constant temperature along an 
%%% isentrope (i.e. Dirichlet BC).

% Creating tri-diagonal matrix
nj = numel(r2);
coeff = diag(r2) + diag(r1,1) + diag(r3,-1);

% Applying boundary conditions:
coeff(1,1)      = 1;
coeff(ind_bc,:) = 0;
coeff(ind_bc + nj.*(ind_bc-1)) = 1;

coeff = sparse(coeff) ; % Making the matrix sparse

function [T,k,rho,Cp,H] = steadystate(z,z_uc,z_lc,z_lm,k,rho,Cp,H,Tp,...
    ind_uc,ind_lc,ind_lm,ind_m,dz,rho_o_m,gamma,rho_o_lm,Ts,clat,crad,alpha,c)
%%% Calculates a 1-D steady-state solution to 
%%% the heat-flow equation with arbitrary 
%%% conductivity, density and specific heat 
%%% capacity that are temperature-dependent.
%%% Uses a constant parameter approximation 
%%% then iteratively updates until a solution
%%% is found where dT<0.1.

% Variables and array set-up:

epsilon = 0.1; %convergence (mean diff)
dT = 50; %random start seed value

% Finding where to apply the lower BC
% of a mantle isentrope.
ind_bc = find(z>=z_uc+z_lc+z_lm);
ind_l = [ind_uc ind_lc ind_lm];

% Calculating best-guess steady-state
Tconst(ind_l) = z(ind_l) * (Tp+z(ind_lm(end))*gamma)/z(ind_lm(end)) ;
Tconst(ind_bc) = Tp + z(ind_bc)*gamma;
T = Tconst;

% Iteratively solving for temperature-dependent 
% material properties based on constant parameter

i = 1; % iteration counter

while dT > epsilon && i <= 10
   

    % Updating material parameters
    
    k(ind_lm) = k_ol(T(ind_lm),clat,crad);
    rho(ind_lm) = rho_ol(T(ind_lm),alpha,rho_o_lm);
    Cp(ind_lm) = Cp_ol(T(ind_lm),c);
    
    k(ind_m) = k_ol(T(ind_m),clat,crad);
    rho(ind_m) = rho_ol(T(ind_m),alpha,rho_o_m);
    Cp(ind_m) = Cp_ol(T(ind_m),c);

    % Defining coefficient vectors

    r1 = [0 (k(2:end-1) + k(3:end))/(2*dz^2)];
    r2 = [0 -1*(2*k(2:end-1) + k(3:end) + k(1:end-2))/(2*dz^2) 0];
    r3 = [(k(2:end-1) + k(1:end-2))/(2*dz^2) 0];

    % Creating L.H.S coefficient matrix

    coeff = maketridiag(r1,r2,r3,ind_bc);

    % Creating R.H.S vector b

    b         = -1*H ;
    b(1)      = Ts   ; % surface temperature
    b(ind_bc) = Tp + z(ind_bc)*gamma ; % Isentrope

    %Solving for the steady-state temperature

    Tnew = coeff\b' ; % Matrix solution
    Tnew = Tnew'    ; % Transpose into row vector

    % Checking solution convergence

    dT = (mean(Tnew-T)^2)^0.5;
    T = Tnew;
    i = i + 1;

end

function [kout] = k_ol(T,clat,crad)
%%% Calculates the thermal conductivity 
%%% of Olivine based on the temperature 
%%% using the expressions of McKenzie '05
krad = crad(1) - crad(2)*(T-273) + crad(3)*((T-273).^2) - ...
    crad(4)*((T-273).^3);

klat = (clat(1)./(1+clat(2)*T)) + clat(3) + ...
    (clat(4)*((T+273))) + (clat(5)*((T+273).^2)) + ...
    (clat(6)*((T+273).^3));

kout = krad + klat;

function [rhoout] = rho_ol(T,alpha,rho_o)
%%% Calculates the density of olivine given
%%% the temperature using the thermal model
%%% and the expansivity model of Bouhifd '96
%%% but ignores the pressure-dependent 
%%% compressability included by Hoggard et al.
t1 = alpha(1)*(T-273);
t2 = alpha(2)*0.5*(T.^2 - 273^2);
rhoout = rho_o*exp(-1*(t1+t2));

function [Cpout] = Cp_ol(T,c)
%%% Calculates the specific heat capacity of the 
%%% mantle materials based on the formulation 
%%% in Korenaga & Korenaga 2016:
Cpout = c(1) - c(2)*(T).^-0.5 - c(3)*(T).^-3;
