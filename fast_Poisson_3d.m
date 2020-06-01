function u_real=fast_Poisson_3d(A,Lx, Ly, Lz, charge,voltage,kx,ky)
Nx = length(charge(:,1,1));
Ny = length(charge(1,:,1));
Nz = length(charge(1,1,:));
dx = Ly/(Ny);

f = zeros(Nx,Ny,Nz);
f(:,:,1) = voltage;
f(:,:,Nz) = 0;
%Fourier spectral method
% k = [0:Nx/2-1 0 -Nx/2+1:1];
% kx=(2*pi/Lx)*[0:(Nx/2) (-Nx/2+1):-1];
% ky=(2*pi/Ly)*[0:(Ny/2) (-Ny/2+1):-1];
% e0 = ones(Nz-2,1);
% e1 = [0;e0(2:end)]; e_1 = [e0(1:end-1);0];
% A = spdiags([e_1,-2*e0,e1],-1:1,Nz-2,Nz-2)/dx^2;
K_xsq = kx.^2;
K_ysq = ky.^2;
F = fft2(f);
% B = ones(Nx+1,Ny+1);
% B = (sin(X/Lx*2*pi)).';
% B = charge;
b = fft2(charge(1:Nx,1:Ny,2:Nz-1));
b(:,:,1) = b(:,:,1) - 10/(12*dx^2).*F(:,:,1);
b(:,:,2) = b(:,:,2) + 1/(12*dx^2).*F(:,:,1);
b(:,:,Nz-2) = b(:,:,Nz-2) - 10/(12*dx^2).*F(:,:,Nz);
b(:,:,Nz-3) = b(:,:,Nz-3) + 1/(12*dx^2) .*F(:,:,Nz);
u = zeros(Nx,Ny,Nz-2);
% tic;
for i = 1:Nx
    for j = 1:Ny
        u(i,j,:) = (A-K_xsq(i)*eye(Nz-2)-K_ysq(j)*eye(Nz-2))\reshape(b(i,j,:),Nz-2,1);        
    end
end
% time1 = toc;
u_real = zeros(Nx,Ny,Nz);
u_real(:,:,2:Nz-1) = real(ifft2(u));
u_real(:,:,1) = voltage;
u_real(:,:,Nz) = 0;
end