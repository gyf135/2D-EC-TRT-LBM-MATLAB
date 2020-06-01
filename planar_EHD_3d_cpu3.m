%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Axis.m: 3D EHD flow with periodic boundary conditions      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lattice Boltzmann sample in Matlab
% Copyright (C) 2017-2018 Yifei Guan
% Address: Mechanical Engineering Building, University of Washington
% E-mail: gyf135@uw.edu
close all;clear all;clc
% GENERAL FLOW CONSTANTS
flag = 10;
if flag==1
time = 0;
Nx     = 40;      % number of cells in x-direction
Ny     = 40;      % number of cells in y-direction
Nz     = 41;      % number of cells in z-direciton (non-periodic)
LL = 1; % The wavelength of perturbation
Lx = LL*1;  % Units in m
Ly = LL*1;  % Units in m
Lz = 1;  % Units in m
dx = Lx/(Nx);
dy = dx;
dz = dy;
xx = linspace(0,Lx-dx,Nx);
yy = linspace(0,Ly-dy,Ny);
zz = linspace(0,Lz,Nz);
[y,x,z] = meshgrid(yy,xx,zz); % get coordinate of matrix indices
% INITIAL CONDITION: Poiseuille profile at equilibrium
ux = zeros(1,Nx,Ny,Nz);
uy = zeros(1,Nx,Ny,Nz);
uz = zeros(1,Nx,Ny,Nz);
phi = zeros(Nx,Ny,Nz);
Ex = zeros(1,Nx,Ny,Nz);
Ey = zeros(1,Nx,Ny,Nz);
Ez = zeros(1,Nx,Ny,Nz);
charge = zeros(1,Nx,Ny,Nz);

CFL = 0.01; %Larger CFL, smaller Omega
dt = dx*CFL; 
cs_square = 1/3/CFL^2;
% uMax   =   1;      % maximum velocity of Poiseuille inflow
rho0 = 1600; % Constant density
K = 2.5e-5;
% Re     = uMax*Ly/nu;      % Reynolds number
Kx=(2*pi/Lx)*[0:(Nx/2) (-Nx/2+1):-1];
Ky=(2*pi/Ly)*[0:(Ny/2) (-Ny/2+1):-1];
% D3Q27 Lattice constants
t = [8/27, 2/27, 2/27, 2/27, 2/27, 2/27, 2/27, 1/54, 1/54, 1/54, 1/54,1/54,...
1/54, 1/54, 1/54, 1/54,1/54,1/54, 1/54,...
1/216,1/216,1/216,1/216,1/216,1/216,1/216,1/216];
%     1 2  3 4  5 6  7 8  910 11 12 13 14 15 16 17 181920 21 22 23 24 25 26
cx = [0,1,-1,0, 0,0, 0,1,-1,1,-1,0, 0, 1,-1, 1,-1, 0, 0,1,-1, 1,-1, 1,-1,-1, 1]/CFL;
cy = [0,0, 0,1,-1,0, 0,1,-1,0, 0,1,-1,-1, 1, 0, 0, 1,-1,1,-1, 1,-1,-1, 1, 1,-1]/CFL;
cz = [0,0, 0,0, 0,1,-1,0, 0,1,-1,1,-1, 0, 0,-1, 1,-1, 1,1,-1,-1, 1, 1,-1, 1,-1]/CFL;
opp = [1,3,2,5,4,7,6,9,8,11,10,13,12,15,14,17,16,19,18,21,20,23,22,25,24,27,26];
col = 2:Nz-1;

for i=1:27
    cu = (cx(i)*ux+cy(i)*uy + cz(i)*uz)/cs_square;
    fIn(i,:,:,:) = rho0 .* t(i) .* ...
                   ( 1 + cu + 1/2*(cu.*cu) - 1/2*(ux.^2+uy.^2+uz.^2)/cs_square);
    ceu = (cx(i)*(K.*Ex+ux)+cy(i)*(K.*Ey+uy)+cz(i)*(K.*Ez+uz))/cs_square;
    hIn(i,:,:,:) = charge .* t(i) .* (1 + ceu + 1/2*(ceu.*ceu)...
        - 1/2*((K.*Ex+ux).^2+(K.*Ey+uy).^2+(K.*Ez+uz).^2)/cs_square );
end
% Create differential matrix A
e0 = ones(Nz,1);
e1 = [0;e0(2:end)]; e_1 = [e0(1:end-1);0];
e2 = [0;0;e0(3:end)]; e_2 = [e0(1:end-2);0;0];
A = spdiags([-1*e_2,16*e_1,-30*e0,16*e1,-1*e_2],-2:2,Nz-2,Nz-2);
A(1,1) = -15; A(1,2) = -4; A(1,3) = 14; A(1,4) = -6; A(1,5) = 1; 
A(Nz-2,Nz-2) = -15; A(Nz-2,Nz-2-1) = -4; A(Nz-2,Nz-2-2) = 14; A(Nz-2,Nz-2-3) = -6; A(Nz-2,Nz-2-4) = 1;
A = A/(12*dx^2);
% Create differential matrix B for solving Ez from phi
B = spdiags([1/12*e_2,-2/3*e_1,0*e0,2/3*e1,-1/12*e2],-2:2,Nz,Nz);
% B(1,1) = -25/12; B(1,2) = 4; B(1,3) = -3; B(1,4) = 4/3; B(1,5) = -1/4;
B(2,1) = -1/4;B(2,2)=-5/6;B(2,3)=3/2;B(2,4)=-1/2;B(2,5)=5/60;
% B(end,end) = 25/12;B(end,end-1) = -4; B(end,end-2) = 3; B(end,end-3) = -4/3;B(end,end-4) = 1/4;
B(end-1,end) = 1/4;B(end-1,end-1) = 5/6; B(end-1,end-2) = -3/2; B(end-1,end-3) = 1/2;B(end-1,end-4) = -5/60;
B=B/dx;
B(1,:) = 0;
B(end,:)=0;

else
    load('pqfile_3d.mat');
end
CFL = 0.01; %Larger CFL, smaller Omega
dt = dx*CFL; 
cs_square = 1/3/CFL^2;
% D3Q27 Lattice constants
t = [8/27, 2/27, 2/27, 2/27, 2/27, 2/27, 2/27, 1/54, 1/54, 1/54, 1/54,1/54,...
1/54, 1/54, 1/54, 1/54,1/54,1/54, 1/54,...
1/216,1/216,1/216,1/216,1/216,1/216,1/216,1/216];
%     1 2  3 4  5 6  7 8  910 11 12 13 14 15 16 17 181920 21 22 23 24 25 26
cx = [0,1,-1,0, 0,0, 0,1,-1,1,-1,0, 0, 1,-1, 1,-1, 0, 0,1,-1, 1,-1, 1,-1,-1, 1]/CFL;
cy = [0,0, 0,1,-1,0, 0,1,-1,0, 0,1,-1,-1, 1, 0, 0, 1,-1,1,-1, 1,-1,-1, 1, 1,-1]/CFL;
cz = [0,0, 0,0, 0,1,-1,0, 0,1,-1,1,-1, 0, 0,-1, 1,-1, 1,1,-1,-1, 1, 1,-1, 1,-1]/CFL;
opp = [1,3,2,5,4,7,6,9,8,11,10,13,12,15,14,17,16,19,18,21,20,23,22,25,24,27,26];
col = 2:Nz-1;

charge0 = 10;
diffu = 6.25e-5;
nu = 0.1;
voltage = 1e4;
eps = 1/voltage;
figure
subplot(3,2,1);
title('Velocity field');
hold on;
subplot(3,2,2);
title('Ux');
hold on;
subplot(3,2,3);
title('Uz');
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
maxT=100;
tPlot = 10;
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
tic % measure time
for cycle = 0:maxT
%     uw = uw + 0.06;
    if mod(cycle,1000)==0
        display(cycle);
    end
    time = time + dt;
    charge = sum(hIn);
%     charge(1,:,1) = charge(1,:,2);%+(charge0-charge(1,:,2));
    charge(1,:,:,Nz) = charge(1,:,:,Nz-1);
    % Calculate body force
    forcex = charge.*Ex;
    forcex(1,:,:,2:Nz) = forcex(1,:,:,2:Nz)+0;
    forcey = charge.*Ey;
    forcez = charge.*Ez;
    % MACROSCOPIC VARIABLES
    rho = sum(fIn);
%     rho(1,(Nx+1)/2,(Ny+1)/2) = rho0;
    ux  = reshape ( (cx * reshape(fIn,27,Nx*Ny*Nz)), 1,Nx,Ny,Nz) ./rho + forcex*dt/2./rho;
    uy  = reshape ( (cy * reshape(fIn,27,Nx*Ny*Nz)), 1,Nx,Ny,Nz) ./rho + forcey*dt/2./rho;
    uz  = reshape ( (cz * reshape(fIn,27,Nx*Ny*Nz)), 1,Nx,Ny,Nz) ./rho + forcez*dt/2./rho;
    
    ux(1,:,:,1) = -ux(1,:,:,2);
    uy(1,:,:,1) = -uy(1,:,:,2);
    uz(1,:,:,1) = -uz(1,:,:,2);
%     ux(1,:,:,Nz) = uw;
%     uy(1,:,:,Nz) = 0;
%     uz(1,:,:,Nz) = 0;

    
    
    %Preturbation
if (cycle == 0)
% =========================================================================
% Preturbation for roll patterns
% =========================================================================
% uz = reshape((cos(2*pi*z)-1).*cos(2*pi/LL*x),1,Nx,Ny,Nz)*1;
% ux = reshape(LL*sin(2*pi*z).*sin(2*pi/LL*x),1,Nx,Ny,Nz)*1;
% charge = charge + reshape((cos(2*pi*z)-1).*cos(2*pi/LL*x),1,Nx,Ny,Nz)*0.01;
% =========================================================================
% Perturbation for square patterns
% =========================================================================
% uz = reshape((cos(2*pi*z)-1).*cos(2*pi/LL*x).*cos(2*pi/LL*y),1,Nx,Ny,Nz);
% ux = reshape(2*LL*sin(2*pi*z).*sin(2*pi/LL*x).*cos(2*pi/LL*y),1,Nx,Ny,Nz);
% uy = reshape(2*LL*sin(2*pi*z).*sin(2*pi/LL*y).*cos(2*pi/LL*x),1,Nx,Ny,Nz);
% =========================================================================
% Perturbation for hexagon patterns
% =========================================================================
% L = 3/(4*pi);
% a = 4*pi/3/L;
% uz = (cos(2*pi*z)-1)/3.*(2*cos(2*pi/(sqrt(3)*L)*x).*cos(2*pi/(3*L)*y)+cos(4*pi/(3*L)*y));
% uz = (cos(2*pi*z)-1)/3.*(2*cos(2*pi/1*x*10).*cos(2*pi/2*y*10)+cos(2*pi*3/2*y*10));
% uz = reshape(uz,1,Nx,Ny,Nz);
% ux = 2*pi*sin(2*pi*z)*(4*pi)/(3*sqrt(3)*L*a^2).*sin(2*pi/(sqrt(3)*L)*x).*cos(2*pi/(3*L)*y);
% ux = reshape(ux,1,Nx,Ny,Nz);
% uy = 2*pi*sin(2*pi*z)*(4*pi)/(9*L*a^2).*(cos(2*pi/(sqrt(3)*L)*x)+2*cos(2*pi/(3*L)*y)).*sin(2*pi/(3*L)*y);
% uy = reshape(uy,1,Nx,Ny,Nz);


% Preturbation for square patterns
% ux(1,(Nx)/2,1,(Nz)/2) = 0.3;
% uy(1,1,(Ny)/2,(Nz)/2) = 0.3;
% Preturbation for hexagon patterns
% uy(1,(Nx)/2,1,(Nz)/2) = 0.3;
% ux(1,1,(Ny)/2,(Nz)/2) = 0.3;
% ux = ux+reshape(awgn(ones(Nx*Ny*Nz,1)*10^-1,10,'measured'),1,Nx,Ny,Nz);
% uy = uy+reshape(awgn(ones(Nx*Ny*Nz,1)*10^-1,10,'measured'),1,Nx,Ny,Nz);
% uy(1,:,(Ny+1)/2,(Nz+1)/2) = 0.3;
% uy(1,:,1,(Nz+1)/2) = 0.3;
% uz(1,(Nx+1)/2,(Ny+1)/2,(Nz+1)/2) = 0.03;
end
 % COLLISION STEP
    for i=1:27
       % Compute external forcing term
        cuu = cx(i)*ux(1,:,:,:)+cy(i)*uy(1,:,:,:)+cz(i)*uz(1,:,:,:);
       forcepop(i,:,:,:)=t(i)/cs_square.*(((cx(i)-ux(1,:,:,:))+cuu*cx(i)/cs_square).*forcex...
       +((cy(i)-uy(1,:,:,:)) + cuu*cy(i)/cs_square).*forcey...
       +((cz(i)-uz(1,:,:,:)) + cuu*cz(i)/cs_square).*forcez);
      
       cu = (cx(i)*ux+cy(i)*uy+ cz(i)*uz)/cs_square;
       fEq(i,:,:,:) = rho .* t(i) .* ...
                   ( 1 + cu + 1/2*(cu.*cu) - 1/2*(ux.^2+uy.^2+uz.^2)/cs_square);
               
%        fEq(i,:,:) = rho.*t(i) + t(i)*rho0/cs_square*(cx(i)*ux+cy(i)*uy+...
%            (ux.^2.*(cx(i).^2-cs_square)+ux.*uy.*(cx(i).*cy(i))+...
%            (uy.^2.*(cy(i).^2-cs_square)))/2/cs_square);
%        fOut2(i,:,:) = fIn(i,:,:) - omega .* (fIn(i,:,:)-fEq(i,:,:))+sourcepop(i,:,:)*dt;
       ceu = (cx(i)*(K.*Ex+ux)+cy(i)*(K.*Ey+uy)+cz(i)*(K.*Ez+uz))/cs_square;
       hEq(i,:,:,:) = charge .* t(i) .* (1 + ceu + 1/2*(ceu.*ceu) ...
           - 1/2*((K.*Ex+ux).^2 + (K.*Ey+uy).^2 + (K.*Ez+uz).^2)/cs_square);                
%        hOut(i,:,:) = hIn(i,:,:) - omega_charge .* (hIn(i,:,:) - hEq(i,:,:));
    end  
    
    for i = 1:27
       forcepop_plus = (forcepop(i,:,:,:) + forcepop(opp(i),:,:,:))*0.5;
       forcepop_minus = (forcepop(i,:,:,:) - forcepop(opp(i),:,:,:))*0.5;
       sourcepop(i,:,:,:) = (1-dt/2*omega_plus)*forcepop_plus...
       + (1-dt/2*omega_minus)*forcepop_minus;
       f_plus = 0.5*(fIn(i,:,:,:)+fIn(opp(i),:,:,:));
       f_minus = 0.5*(fIn(i,:,:,:)-fIn(opp(i),:,:,:));
       f_eq_plus = 0.5*(fEq(i,:,:,:)+fEq(opp(i),:,:,:));
       f_eq_minus = 0.5*(fEq(i,:,:,:)-fEq(opp(i),:,:,:));
%        fOut(i,:,:) = fIn(i,:,:) - omega .* (fIn(i,:,:)-fEq(i,:,:))+sourcepop(i,:,:)*dt;
       fOut(i,:,:,:) = fIn(i,:,:,:) - omega_plus*dt*(f_plus-f_eq_plus)...
           -omega_minus*dt*(f_minus-f_eq_minus)+sourcepop(i,:,:,:)*dt;
       h_plus = 0.5*(hIn(i,:,:,:)+hIn(opp(i),:,:,:));
       h_minus = 0.5*(hIn(i,:,:,:)-hIn(opp(i),:,:,:));
       h_eq_plus = 0.5*(hEq(i,:,:,:)+hEq(opp(i),:,:,:));
       h_eq_minus = 0.5*(hEq(i,:,:,:)-hEq(opp(i),:,:,:));
%        hOut(i,:,:) = hIn(i,:,:) - omega_charge .* (hIn(i,:,:) - hEq(i,:,:));
       hOut(i,:,:,:) = hIn(i,:,:,:) - omega_c_plus*dt*(h_plus - h_eq_plus)...
            -omega_c_minus*dt*(h_minus-h_eq_minus);
    end
    
    % ==================================================================
%     display((forcepop(21,1,1,1) - forcepop(opp(21),1,1,1))*0.5);
    % ==================================================================
       
    % OBSTACLE (Full-WAY BOUNCE-BACK)
    for i=1:27
         fOut(i,:,:,1) = fIn(opp(i),:,:,1);
         fOut(i,:,:,Nz) = fIn(opp(i),:,:,Nz) + 2*t(i)*rho0*(cx(i)*uw)/cs_square;
    end
    
    % Zero gradient on Nz
    for i = 1:27
        hOut(i,:,:,Nz) = hOut(i,:,:,Nz-1);
    end
    % STREAMING STEP
    for i=1:27
       fIn(i,:,:,:) = circshift(fOut(i,:,:,:), [0,cx(i),cy(i),cz(i)]*CFL);
       hIn(i,:,:,:) = circshift(hOut(i,:,:,:), [0,cx(i),cy(i),cz(i)]*CFL);
    end
    % Boundary condition for charges
    for i = 1:27
        hIn(opp(i),:,:,1) = -hOut(i,:,:,1)+2*(charge0) * t(i);
    end

    
    % Fast Poisson solver
%     tic
%     [phi] = fast_Poisson_3d(A,Lx, Ly, Lz, reshape(-charge/eps, Nx, Ny, Nz),voltage,Kx,Ky); 
    [phi] = fast_Poisson_3d_v2(Nx,Ny,Nz,Lx,Ly,Lz,reshape(-charge/eps, Nx, Ny, Nz),voltage);
%time2 = toc
%display(time2);
%     phi_hat = fft2(phi);
%     for j = 1:Nx
%         Ex_hat(j,:,:) = 1i*Kx(j)*phi_hat(j,:,:);
%     end
%     Ex_hat(Nx/2+1,:,:) = 0;
%     Ex(1,:,:,:) = -reshape(real(ifft2(Ex_hat)),1,Nx,Ny,Nz);
%     for j = 1:Ny
%         Ey_hat(:,j,:) = 1i*Ky(j)*phi_hat(:,j,:);
%     end
%     Ey_hat(:,Ny/2+1,:) = 0;
%     Ey(1,:,:,:) = -reshape(real(ifft2(Ey_hat)),1,Nx,Ny,Nz);
%     
%     for j = 1:Ny
%         for i = 1:Nx
%             Ez(1,i,j,:) = -B*reshape(phi(i,j,:),Nz,1);
%         end
%     end
    
    
    Ex(1,2:end-1,:,:) = -(phi(3:end,:,:)-phi(1:end-2,:,:))/(2*dx);
    Ex(1,1,:,:) = -(phi(2,:,:)-phi(end,:,:))/(2*dx);
    Ex(1,end,:,:) = -(phi(1,:,:)-phi(end-1,:,:))/(2*dx);
    
    Ey(1,:,2:end-1,:) = -(phi(:,3:end,:)-phi(:,1:end-2,:))/(2*dx);
    Ey(1,:,1,:) = -(phi(:,2,:) - phi(:,end,:))/(2*dx);
    Ey(1,:,end,:) = -(phi(:,1,:) - phi(:,end-1,:))/(2*dx);
    
    Ez(1,:,:,2:end-1) = -(phi(:,:,3:end) - phi(:,:,1:end-2))/(2*dx);   
    Ez(1,:,:,1) = Ez(1,:,:,2);%*2 - Ey(1,:,:,3);
    Ez(1,:,:,end) = Ey(1,:,:,end-1);%*2 - Ey(1,:,:,end-2);
%      
  
    if (mod(cycle,tPlot)==1 && cycle ~=1)
        if (mod(cycle,500)==1) 
            save('pqfile_3d.mat')
            display('saving data..........................');
        end
        % Store uy for further analysis on modes
%         analysis = [analysis reshape(uy(1,:,:),(Nx)*(Ny),1)];   
        subplot(3,2,1);
        u = sqrt(ux.^2+uy.^2+uz.^2);
        u = reshape(u(1,:,:,(Nz+1)/2),Nx,Ny);
%         u = u(1:10,:);
%         u(bbRegion) = nan;
        imagesc(u');
        axis equal off; drawnow
      %  title('Velocity field');
        subplot(3,2,3);
         u = reshape(uz(1,:,:,(Nz+1)/2),Nx,Ny);
%         u(bbRegion) = 1;
        imagesc(u');
        axis equal off; drawnow
      %  title('Charge density');
        subplot(3,2,5);
        u = reshape(phi(:,:,(Nz+1)/2),Nx,Ny);
%         u = u(1:10,:);
        imagesc(u');
        axis equal off; drawnow
       % title('Electric potential');
       subplot(3,2,2);
        u = reshape(ux(1,:,:,(Nz+1)/2),Nx,Ny);
%         u(bbRegion) = nan;
        imagesc(u');
        axis equal off; drawnow
       % title('Electric potential');
       subplot(3,2,4);
        u = reshape(uy(1,:,:,(Nz+1)/2),Nx,Ny);
%         u = u(1:10,:);
%         u(bbRegion) = nan;
        imagesc(u');
        axis equal off; drawnow
       % title('Electric potential');
              subplot(3,2,6);
        u = reshape(charge(1,:,:,(Nz+1)/2),Nx,Ny);
%         u = u(1:10,:);
%         u(bbRegion) = nan;
        imagesc(u');
        axis equal off; drawnow
       % title('Electric potential');
    end
    if (mod(cycle,200) == 1)
        analysis = [analysis reshape(uz(1,:,:,:),Nx*Ny*Nz,1)];
        maxW = [maxW;time max(max(abs(uz(1,:,:,(Nz+1)/2))))];
    end
    
end
time2 = toc;
save('pqfile_3d.mat')
save('saveA.mat','analysis');
%%
%%
figure
test = reshape(charge(1,:,Ny/2,:),Nx,Nz);
surf(test)
title('Charge');
figure
surf(reshape(ux(1,:,Ny/2,:),Nx,Nz));
title('Ux');
figure
surf(reshape(uy(1,:,Ny/2,:),Nx,Nz));
title('Uy');
figure
surf(reshape(uz(1,:,Ny/2,:),Nx,Nz));
title('Uz');
figure
surf(reshape(K.*Ez(1,:,Ny/2,:),Nx,Nz));
title('Drift velocity z');
%%
figure
% h=slice(x,y,z,reshape(ux,Nx,Ny,Nz),[Lx/2],[Ly/2],[Lz/2]);
h=slice(y,x,z,reshape(rho,Nx,Ny,Nz),[Ly/2],[Lx/2],[Lz/2]);
% set(h,'edgecolor','none')
hold on;
colorbar;
ylabel('x');
xlabel('y');
zlabel('z');
title('Pressure');
figure
% h=slice(x,y,z,reshape(ux,Nx,Ny,Nz),[Lx/2],[Ly/2],[Lz/2]);
h=slice(y,x,z,reshape(ux,Nx,Ny,Nz),[Ly/2],[Lx/2],[Lz/2]);
% set(h,'edgecolor','none')
hold on;
colorbar;
ylabel('x');
xlabel('y');
zlabel('z');
title('ux');
figure
% slice(x,y,z,reshape(uy,Nx,Ny,Nz),[Lx/2],[Ly/2],[Lz/2]);
slice(y,x,z,reshape(uy,Nx,Ny,Nz),[Ly/2],[Lx/2],[Lz/2]);
hold on;
colorbar;
ylabel('x');
xlabel('y');
zlabel('z');
title('uy');
figure
% slice(x,y,z,reshape(uz,Nx,Ny,Nz),[Lx/2],[Ly/2],[Lz/2]);
slice(y,x,z,reshape(uz,Nx,Ny,Nz),[Ly/2],[Lx/2],[Lz/2]);
hold on;
colorbar;
ylabel('x');
xlabel('y');
zlabel('z');
title('uz');
figure
% slice(x,y,z,reshape(charge(1,:,:,:),Nx,Ny,Nz),[Lx/2],[Ly/2],[Lz/2]);
slice(y,x,z,reshape(charge(1,:,:,:),Nx,Ny,Nz),[Ly/2],[Lx/2],[Lz/2]);
hold on;
colorbar;
ylabel('x');
xlabel('y');
zlabel('z');
title('Charge density');
figure
plot(maxW(:,1),maxW(:,2));
display(time2);