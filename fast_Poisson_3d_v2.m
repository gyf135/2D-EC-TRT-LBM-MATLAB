function u=fast_Poisson_3d_v2(Nx,Ny,Nz,Lx,Ly,Lz,charge,voltage)
% Nx = 100;
% Ny = 101; % Interior points m = Ny-2
% Lx = 1;
% Ly = 1;
dx = Lx/Nx;
dy = dx;
dz = dx;
% xx = linspace(0,Lx-dx,Nx);
% yy = linspace(dy,Ly-dy,Ny-2);
% [y,x] = meshgrid(yy,xx); % get coordinate of matrix indices
Ne = 2*(Nz-1); % Ne = 2*(m+1)
% ==========================
% tic;
% ==========================
% f  = 4*ones(Ny-2,Nx);
% f  = 40*(sin(x/Lx*2*pi)).';
f = charge(:,:,2:Nz-1);
f(:,:,1) = f(:,:,1) - voltage / dz^2;
g = cat(3,zeros(Nx,Ny,1),f, zeros(Nx,Ny,1),-f(:,:,Nz-2:-1:1)); % odd extension in z
g_hat = fftn(g);
kx = 2*pi/Lx*[0:(Nx/2) (-Nx/2+1):-1];
ky = 2*pi/Lx*[0:(Ny/2) (-Ny/2+1):-1];
kz = 2*pi/(Ne*dy)*[0:(Ne/2) (-Ne/2+1):-1];
[J,I,K]  = meshgrid(ky,kx,kz);
mu = (4/dz^2)*(sin(K*dz/2).^2) + I.^2 + J.^2;
mu(1,1,1) = 1;   % Avoid 0/0 division; vhat(1,1) is known a priori to be 0
v = real(ifftn(-g_hat./mu));
u = v(:,:,1:Nz); % Extract u(i,j) from its odd extension
u(:,:,1) = voltage;
% u = u.';
% ==========================
% time1 = toc;
% ==========================
%  Plot out solution in interior and print out max-norm error
% surf(u);
% xlabel('x')
% ylabel('y')
% zlabel('u')
% title('FD/FFT for 2D Poisson Eqn')
end
