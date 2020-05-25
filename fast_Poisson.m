function u_real=fast_Poisson(A, Lx, Ly, charge,voltage,kx)
Nx = length(charge(:,1));
Ny = length(charge(1,:));
dx = Lx/(Nx);
% x = linspace(0,1-dx,Nx);
% y = linspace(0,2-dx,Ny);
% [X,Y] = meshgrid(x,y);
f = zeros(Nx,Ny);
f(:,1) = voltage;
f(:,Ny) = 0;
%Fourier spectral method
% k = [0:Nx/2-1 0 -Nx/2+1:1];
% kx=(2*pi/Lx)*[0:(Nx/2) (-Nx/2+1):-1];
% e0 = ones(Ny-2,1);
% e1 = [0;e0(2:end)]; e_1 = [e0(1:end-1);0];
% A = spdiags([e_1,-2*e0,e1],-1:1,Ny-2,Ny-2)/dx^2;
K_sq = kx.^2;
F = fft(f);
% B = ones(Nx+1,Ny+1);
% B = (sin(X/Lx*2*pi)).';
% B = charge;
% B = zeros(size(B));
b = fft(charge(1:Nx,2:Ny-1));
% b = b.';
b(:,1) = b(:,1)  - 10/(12*dx^2).*F(:,1);
b(:,2) = b(:,2) + 1/(12*dx^2).*F(:,1);
b(:,Ny-2) = b(:,Ny-2) - 10/(12*dx^2).*F(:,Ny);
b(:,Ny-3) = b(:,Ny-3) + 1/(12*dx^2) .*F(:,Ny);
u = [];
% tic;
for i = 1:Nx
    u = [u (A-K_sq(i)*eye(Ny-2))\b(i,:).'];
end
% time1 = toc;
u_real = ifft(u.');
u_real = real([f(:,1) u_real f(:,Ny)]);
end