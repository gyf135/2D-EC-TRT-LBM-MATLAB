%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Axis.m: Planar EHD flow with periodic boundary conditions      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lattice Boltzmann sample in Matlab
% Copyright (C) 2017-2018 Yifei Guan
% Address: Mechanical Engineering Building, University of Washington
% E-mail: gyf135@uw.edu
close all;clear all;clc
% GENERAL FLOW CONSTANTS
flag = 10; % If flag = 1, initialize the problem
perturb = 10; %% If perturb = 1, perturb the flow
if flag==1
time = 0;
Nx     = 122;      % number of cells in x-direction
Ny     = 101;      % number of cells in y-direction
LL = 1.22;  % The wavelength of perturbation
Lx = LL;  % Units in m
Ly = 1;  % Units in m
dx = Lx/(Nx);
dy = dx;
xx = linspace(0,Lx-dx,Nx);
yy = linspace(0,Ly,Ny);
[y,x] = meshgrid(yy,xx); % get coordinate of matrix indices
% INITIAL CONDITION: Poiseuille profile at equilibrium
ux = zeros(1,Nx,Ny);
uy = zeros(1,Nx,Ny);
phi = zeros(Nx,Ny);
Ex = zeros(1,Nx,Ny);
Ey = zeros(1,Nx,Ny);
Ex_hat = phi;
charge = zeros(1,Nx,Ny);
CFL = 0.05; %Larger CFL, smaller Omega
dt = dx*CFL; 
% nu = 0.1470588;
% dt = 0.5*dx^2/(3*nu);
% CFL = dt/dx;
cs_square = 1/3/CFL^2;
rho0 = 1600; % Constant density
K = 2.5e-5;
% D2Q9 LATTICE CONSTANTS
t  = [4/9, 1/9,1/9,1/9,1/9, 1/36,1/36,1/36,1/36];
cx = [  0,   1,  0, -1,  0,    1,  -1,  -1,   1]/CFL;
cy = [  0,   0,  1,  0, -1,    1,   1,  -1,  -1]/CFL;
opp = [ 1,   4,  5,  2,  3,    8,   9,   6,   7];
for i=1:9
    cu = (cx(i)*ux+cy(i)*uy)/cs_square;
    fIn(i,:,:) = rho0 .* t(i) .* ...
                   ( 1 + cu + 1/2*(cu.*cu) - 1/2*(ux.^2+uy.^2)/cs_square);
    ceu = (cx(i)*(K.*Ex+ux)+cy(i)*(K.*Ey+uy))/cs_square;
    hIn(i,:,:) = charge .* t(i) .* (1 + ceu + 1/2*(ceu.*ceu)...
        - 1/2*((K.*Ex+ux).^2+(K.*Ey+uy).^2)/cs_square );
end
Kx=(2*pi/Lx)*[0:(Nx/2) (-Nx/2+1):-1];
% Create differential matrix A
e0 = ones(Ny,1);
e1 = [0;e0(2:end)]; e_1 = [e0(1:end-1);0];
e2 = [0;0;e0(3:end)]; e_2 = [e0(1:end-2);0;0];
A = spdiags([-1*e_2,16*e_1,-30*e0,16*e1,-1*e_2],-2:2,Ny-2,Ny-2);
A(1,1) = -15; A(1,2) = -4; A(1,3) = 14; A(1,4) = -6; A(1,5) = 1;
A(Ny-2,Ny-2) = -15; A(Ny-2,Ny-2-1) = -4; A(Ny-2,Ny-2-2) = 14; A(Ny-2,Ny-2-3) = -6; A(Ny-2,Ny-2-4) = 1;
A = A/(12*dx^2);
% Create differential matrix B for solving Ez from phi
B = spdiags([1/12*e_2,-2/3*e_1,0*e0,2/3*e1,-1/12*e2],-2:2,Ny,Ny);
% B(1,1) = -25/12; B(1,2) = 4; B(1,3) = -3; B(1,4) = 4/3; B(1,5) = -1/4;
B(2,1) = -1/4;B(2,2)=-5/6;B(2,3)=3/2;B(2,4)=-1/2;B(2,5)=5/60;
% B(end,end) = 25/12;B(end,end-1) = -4; B(end,end-2) = 3; B(end,end-3) = -4/3;B(end,end-4) = 1/4;
B(end-1,end) = 1/4;B(end-1,end-1) = 5/6; B(end-1,end-2) = -3/2; B(end-1,end-3) = 1/2;B(end-1,end-4) = -5/60;
B=B/dx;
B(1,:) = 0;
B(end,:)=0;

% B = spdiags([-1*e_1,0*e0,1*e1],-1:1,Ny,Ny)/dx*0.5;
% B(1,:) = 0;
% B(end,:)=0;


else
    load('pqfile.mat');
end

CFL =   0.05; %Larger CFL, smaller Omega
dt = dx*CFL; 
nu = 0.1;
% dt = 0.5*dx^2/(3*nu);
% CFL = dt/dx;
cs_square = 1/3/CFL^2;
% uMax   =   1;      % maximum velocity of Poiseuille inflow

% Re     = uMax*Ly/nu;      % Reynolds number
% omega  = 1. / (nu/dt/cs_square+1./2.);      % relaxation parameter\

%omega = 1.0;
% D2Q9 LATTICE CONSTANTS
% t  = [4/9, 1/9,1/9,1/9,1/9, 1/36,1/36,1/36,1/36];
cx = [  0,   1,  0, -1,  0,    1,  -1,  -1,   1]/CFL;
cy = [  0,   0,  1,  0, -1,    1,   1,  -1,  -1]/CFL;
% opp = [ 1,   4,  5,  2,  3,    8,   9,   6,   7];
col = 2:(Ny-1);
in  = 1;   % position of inlet
out = Nx;  % position of outlet
charge0 = 10;
diffu = 6.25e-5;
% omega_charge = 1. / (diffu/cs_square/dt+1./2.);
% omega_c_minus = omega_c_plus;
voltage = 1e4;
eps = 1/voltage;

% obst = zeros(Nx,Ny);
% obst(:,[1,Ny]) = 1;    % Location of top/bottom boundary
% bbRegion = find(obst);
figure
subplot(3,2,1);
title('Velocity field');
hold on;
subplot(3,2,2);
title('Ux');
hold on;
subplot(3,2,3);
title('Pressure');
hold on;
subplot(3,2,4);
title('Uy');
hold on;
subplot(3,2,5);
title('Electric potential');
hold on;
subplot(3,2,6);
title('Charge density');
hold on;

% MAIN LOOP (TIME CYCLES)
maxT=50000;
tPlot = 200;
uw = 0.0;
V=1/12;
omega_plus = 1/(nu/cs_square/dt+1/2)/dt;
omega_minus = 1/(V/(nu/cs_square/dt)+1/2)/dt;
% omega_minus = omega_plus;
% tau = dt/omega;
% disp(['Omega=',num2str(omega)]);
VC=1.0e-6;
omega_c_minus = 1/(diffu/cs_square/dt+1/2)/dt;
omega_c_plus = 1/(VC/(diffu/cs_square/dt)+1/2)/dt;
% obst([1,Nx],:) = 1;    % Location of top/bottom boundary
% bbRegion = find(obst);
% Store values for further analysis
analysis =[];
maxW = [];
M = sqrt(eps/rho0)/K;
T = eps*voltage/K/nu/rho0;
C = charge0 * Ly^2/(voltage*eps);
Fe = K*voltage/diffu;
tic
for cycle = 0:maxT
%     uw = uw + 0.06;
    if mod(cycle,tPlot)==0
        display(cycle);
%         display('Current on the upper plate =');
%         display( K * dx * sum(sum(sum(charge .* Ey))));
    end
    time = time + dt;
    charge = sum(hIn);
%     charge(1,:,1) = charge(1,:,2);%+(charge0-charge(1,:,2));
    charge(1,:,Ny) = charge(1,:,Ny-1);
    % Calculate body force
    forcex = 0 + charge.*Ex;
    forcex(1,:,1:Ny) = forcex(1,:,1:Ny);% + 3000;
    forcey = 0 + charge.*Ey;
    % MACROSCOPIC VARIABLES
    rho = sum(fIn);
%     rho(1,(Nx+1)/2,(Ny+1)/2) = rho0;
ux_old = ux;
uy_old = uy;
    ux  = reshape ( (cx * reshape(fIn,9,Nx*Ny)), 1,Nx,Ny) ./rho + forcex*dt/2./rho;
    uy  = reshape ( (cy * reshape(fIn,9,Nx*Ny)), 1,Nx,Ny) ./rho + forcey*dt/2./rho;
    ux(1,:,1) = -ux_old(1,:,2);
    uy(1,:,1) = -uy_old(1,:,2);
%     ux(1,:,Ny) = uw;%2*uw - ux(1,:,Ny-1);
%     uy(1,:,Ny) = -uy_old(1,:,Ny-1);
%     ux(1,1,:) = 0;
%     uy(1,1,:) = 0;
%     ux(1,Nx,:) = uw;
%     uy(1,Nx,:) = 0;
% ux(1,:,:) = 0;
% uy(1,:,:) = 0;
if (cycle == 0)
if (perturb == 1)
% =========================================================================
% Perturbation for roll patterns
% =========================================================================
uy = reshape((cos(2*pi*y)-1).*cos(2*pi/LL*x),1,Nx,Ny)*0.001;
ux = reshape(LL*sin(2*pi*y).*sin(2*pi/LL*x),1,Nx,Ny)*0.001;
% uy = reshape((sin(2*pi*y)).*cos(2*pi/LL*x),1,Nx,Ny)*0.1;
% ux = reshape(-LL*cos(2*pi*y).*sin(2*pi/LL*x),1,Nx,Ny)*0.1;


% charge = charge + reshape((cos(2*pi*y)-1).*cos(2*pi/LL*x),1,Nx,Ny)*0.1;
end
% =========================================================================
% ux(1,(Nx)/2,(Ny+1)/2) = 0.3;
% ux(1,:,Ny/2) = 0.3;
% uy = reshape((cos(2*pi*y)-1).*cos(2.57*x),1,Nx,Ny);
% ux = reshape(2*pi/2.57*sin(2*pi*y).*sin(2.57*x),1,Nx,Ny);
% uy = reshape((cos(2*pi*y)-1).*cos(2*pi/LL*x),1,Nx,Ny)*0.1;
% ux = reshape(LL*sin(2*pi*y).*sin(2*pi/LL*x),1,Nx,Ny)*0.1;
% charge = charge + reshape((cos(2*pi*y)-1).*cos(2*pi/LL*x),1,Nx,Ny)*0.1;
% for i=1:9
%     cu = (cx(i)*ux+cy(i)*uy)/cs_square;
%     fIn(i,:,:) = rho0 .* t(i) .* ...
%                    ( 1 + cu + 1/2*(cu.*cu) - 1/2*(ux.^2+uy.^2)/cs_square);
%     ceu = (cx(i)*(K.*Ex+ux)+cy(i)*(K.*Ey+uy))/cs_square;
%     hIn(i,:,:) = charge .* t(i) .* (1 + ceu + 1/2*(ceu.*ceu)...
%         - 1/2*((K.*Ex+ux).^2+(K.*Ey+uy).^2)/cs_square );
% end
end
% ux(1,1,(Ny+1)/2) = 0.3;
 % COLLISION STEP
    for i=1:9
       % Compute external forcing term
%        forcepop(i,:,:)=t(i)/cs_square.*(cx(i).*forcex+cy(i).*forcey);
       forcepop(i,:,:)=t(i)/cs_square.*...
       (((cx(i)-ux(1,:,:)) + (cx(i)*ux(1,:,:)+cy(i)*uy(1,:,:))*cx(i)/cs_square).*forcex...
       +((cy(i)-uy(1,:,:)) + (cx(i)*ux(1,:,:)+cy(i)*uy(1,:,:))*cy(i)/cs_square).*forcey);
   % SHOULD BE A PLUS + SIGN
       
%        sourcepop(i,:,:) = (1-dt/(2*tau))*forcepop(i);
       cu = (cx(i)*ux+cy(i)*uy)/cs_square;
       fEq(i,:,:) = rho .* t(i) .* ...
                   ( 1 + cu + 1/2*(cu.*cu) - 1/2*(ux.^2+uy.^2)/cs_square);
               
%        fEq(i,:,:) = rho.*t(i) + t(i)*rho0/cs_square*(cx(i)*ux+cy(i)*uy+...
%            (ux.^2.*(cx(i).^2-cs_square)+ux.*uy.*(cx(i).*cy(i))+...
%            (uy.^2.*(cy(i).^2-cs_square)))/2/cs_square);
%        fOut2(i,:,:) = fIn(i,:,:) - omega .* (fIn(i,:,:)-fEq(i,:,:))+sourcepop(i,:,:)*dt;
       ceu = (cx(i)*(K.*Ex+ux)+cy(i)*(K.*Ey+uy))/cs_square;
       hEq(i,:,:) = charge .* t(i) .* (1 + ceu + 1/2*(ceu.*ceu) ...
           - 1/2*((K.*Ex+ux).^2 + (K.*Ey+uy).^2)/cs_square);                
%        hOut(i,:,:) = hIn(i,:,:) - omega_charge .* (hIn(i,:,:) - hEq(i,:,:));
    end

    for i = 1:9
       forcepop_plus = (forcepop(i,:,:) + forcepop(opp(i),:,:))*0.5;
       forcepop_minus = (forcepop(i,:,:) - forcepop(opp(i),:,:))*0.5;
       sourcepop(i,:,:) = (1-dt/2*omega_plus)*forcepop_plus...
       + (1-dt/2*omega_minus)*forcepop_minus;
       f_plus = 0.5*(fIn(i,:,:)+fIn(opp(i),:,:));
       f_minus = 0.5*(fIn(i,:,:)-fIn(opp(i),:,:));
       f_eq_plus = 0.5*(fEq(i,:,:)+fEq(opp(i),:,:));
       f_eq_minus = 0.5*(fEq(i,:,:)-fEq(opp(i),:,:));
%        fOut(i,:,:) = fIn(i,:,:) - omega .* (fIn(i,:,:)-fEq(i,:,:))+sourcepop(i,:,:)*dt;
       fOut(i,:,:) = fIn(i,:,:) - omega_plus*dt*(f_plus-f_eq_plus)...
           -omega_minus*dt*(f_minus-f_eq_minus)+sourcepop(i,:,:)*dt;
%        if (i==2)
%         test =  omega_plus*dt;
%        display(test);
%        end
       
       
       h_plus = 0.5*(hIn(i,:,:)+hIn(opp(i),:,:));
       h_minus = 0.5*(hIn(i,:,:)-hIn(opp(i),:,:));
       h_eq_plus = 0.5*(hEq(i,:,:)+hEq(opp(i),:,:));
       h_eq_minus = 0.5*(hEq(i,:,:)-hEq(opp(i),:,:));
%        hOut(i,:,:) = hIn(i,:,:) - omega_charge .* (hIn(i,:,:) - hEq(i,:,:));
       hOut(i,:,:) = hIn(i,:,:) - omega_c_plus*dt*(h_plus - h_eq_plus)...
            -omega_c_minus*dt*(h_minus-h_eq_minus);       
    end
%     
%                     disp(charge(1,5,2));

%             disp(fOut(2,5,1));
    
    % OBSTACLE (Full-WAY BOUNCE-BACK)
    for i=1:9
         fOut(i,:,1) = fIn(opp(i),:,1);
         fOut(i,:,Ny) = fIn(opp(i),:,Ny) + 2*t(i)*rho0*(cx(i)*uw)/cs_square;
    end
%     
%     for i=1:9
%          fOut(i,:,Ny) = fIn(opp(i),:,Ny) + 2*t(i)*rho(:,Ny)*(cx(i)*uw)/cs_square;
%     end
    
    % Zero gradient on Ny
    for i = 1:9
        hOut(i,:,Ny) = hOut(i,:,Ny-1);
    end
    % STREAMING STEP
    for i=1:9
       fIn(i,:,:) = circshift(fOut(i,:,:), [0,cx(i),cy(i)]*CFL);
       hIn(i,:,:) = circshift(hOut(i,:,:), [0,cx(i),cy(i)]*CFL);
    end
    % Boundary condition for charges
    for i = 1:9
        hIn(opp(i),:,1) = -hOut(i,:,1)+2*(charge0) * t(i);
    end
    

    %     if i==1

%         disp(Ey(1,5,2));
%     end
%     disp(hIn(1,5,1));
%     display(hIn(1,5,1));
%     display(hIn(2,5,1));
%     disp(fOut(6,4,Ny));
% i = 6;
% display(2*t(i)*rho0*(cx(i)*uw)/cs_square);
    
%     for i = 1:9
%         hIn(opp(i),:,Ny) = -hOut(i,:,Ny)+2*(charge(1,:,Ny)) * t(i);
%     end
    
%     for i = [5,8,9]
%         hIn(i,:,Ny) = hIn(i,:,Ny-1);
%     end
    
    % Fast Poisson solver
%     [phi] = fast_Poisson(phi, Lx, Ly, reshape(-charge/eps, Nx, Ny),voltage); 
%     [phi] = fast_Poisson(A, Lx, Ly, reshape(-charge/eps, Nx, Ny),voltage,Kx); 
% [phi] = slow_Poisson_2(phi, Lx, Ly, reshape(-charge/eps, Nx, Ny),voltage); 
     phi = fast_Poisson_v2(Nx,Ny,Lx, Ly, reshape(-charge/eps, Nx, Ny),voltage);
     phi_hat = fft2(phi);
     for j = 1:Nx
         Ex_hat(j,:) = 1i*Kx(j)*phi_hat(j,:);
         Ey(1,j,:) = -B*reshape(phi(j,:),Ny,1);
     end
     Ex_hat(Nx/2+1,:) = 0;
     Ex(1,:,:) = -reshape(real(ifft2(Ex_hat)),1,Nx,Ny);
%[Ex,Ey] = slow_field(phi,dx);
Ex = reshape(Ex,1,Nx,Ny);
Ey = reshape(Ey,1,Nx,Ny);

%     Ex(1,2:end-1,:) = -(phi(3:end,:)-phi(1:end-2,:))/(2*dx);
%     Ex(1,1,:) = -(phi(2,:)-phi(end,:))/(2*dx);
%     Ex(1,end,:) = -(phi(1,:)-phi(end-1,:))/(2*dx);
%     Ey(1,:,2:end-1) = -(phi(:,3:end)-phi(:,1:end-2))/(2*dx);
   
%     Ey(1,:,1) = Ey(1,:,2);%*2 - Ey(1,:,3);
%     Ey(1,:,end) = Ey(1,:,end-1);%*2 - Ey(1,:,end-2);
        
    if (mod(cycle,tPlot)==1)
        if (mod(cycle,1500) == 1 && cycle ~= 1)
            save('pqfile.mat')
            display('saving..................');
        end
        
        % Store uy for further analysis on modes
        analysis = [analysis reshape(uy(1,:,(Ny+1)/2),Nx,1)];   
        velo = sqrt(uy(1,:,:).^2 + ux(1,:,:).^2); 
        maxW = [maxW;time max(max(velo(1,:,:)))];
        
        
        subplot(3,2,1);
        u = reshape(sqrt(ux.^2+uy.^2),Nx,Ny);
%         u = u(1:10,:);
%         u(bbRegion) = nan;
        imagesc(u');
        axis equal off; drawnow
      %  title('Velocity field');
        subplot(3,2,3);
        u = reshape(rho(1,:,2:Ny-1),Nx,Ny-2);
%         u(bbRegion) = 1;
        imagesc(u');
        axis equal off; drawnow
      %  title('Charge density');
        subplot(3,2,5);
        u = reshape(phi(:,:),Nx,Ny);
%         u = u(1:10,:);
        imagesc(u');
        axis equal off; drawnow
       % title('Electric potential');
       subplot(3,2,2);
        u = reshape(ux(1,:,2:Ny-1),Nx,Ny-2);
%         u(bbRegion) = nan;
        imagesc(u');
        axis equal off; drawnow
       % title('Electric potential');
       subplot(3,2,4);
        u = reshape(uy(1,:,2:Ny-1),Nx,Ny-2);
%         u = u(1:10,:);
%         u(bbRegion) = nan;
        imagesc(u');
        axis equal off; drawnow
       % title('Electric potential');
              subplot(3,2,6);
        u = reshape(charge(1,:,2:Ny-1),Nx,Ny-2);
%         u = u(1:10,:);
%         u(bbRegion) = nan;
        imagesc(u');
        axis equal off; drawnow
       % title('Electric potential');
    end
    
    
end
time1=toc;
%%
profile = reshape(ux(1,out,col),length(col),1);
profile = [profile;0];
profile1 = reshape(rho(1,Nx,col)*cs_square,length(col),1);
profile1 = [profile1;0];
profile2 = reshape(rho(1,1,col)*cs_square,length(col),1);
profile2 = [profile2;0];
disp(max(profile));
disp(max(max(uy)));
disp(max(max(rho)));
% disp('Tau should be close to 1/2');
% disp(['tau=', num2str(dt/omega)]);
save('pqfile.mat')
save('saveA.mat','analysis');
figure
test = reshape(charge(1,:,2:Ny-1),Nx,Ny-2);
surf(test)
title('Charge density');
figure
surf(reshape(ux(1,:,2:Ny-1),Nx,Ny-2));
title('ux');
figure
surf(reshape(uy(1,:,2:Ny-1),Nx,Ny-2));
title('uy');
figure
surf(reshape(K.*Ey(1,:,2:Ny-1),Nx,Ny-2));
title('Drift velocity y');
figure
plot(maxW(:,1),maxW(:,2));
% Current on the upper plate
I = K * dx * sum(sum(charge(1,:,Ny) .* Ey(1,:,Ny)));
display('Current = ');
disp(I);
display(time1)

%%
chargeeee = reshape(charge, Nx,Ny);
chargee = chargeeee(1,:).'/C;
Eyeee = reshape(Ey,Nx,Ny); 
Eyee  = Eyeee(1,:).'*Ly/voltage;
