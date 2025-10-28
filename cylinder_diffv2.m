function xdot = cylinder_diffv2(t,x,D,R,h,N)
% Computes xdot given t, x, D, h, and N
% Finite difference discretization 
% of diffusion from a cylinder
% x - vector of length NM. 
% N, M are assumed even.
% x(iN + j)) contains concentration
% at radius ri, height zj
% where ri = iR/N and zj = jh/M
% where i = 0,..., N - 1
% and j = 0,..., M - 1
% Concentration at R and h is 0.
c = reshape(x, N+1, []); % c(1, 1) -> concentration at r=0, h=0
M = length(c(:,1)) - 1;
delr = R/N;
delz = h/M;
cdot=zeros(N + 1,M + 1); % cdot(1, 1) -> time derivative of concentration at r=0,h=0
Dr = D/(delr^2);
Dz = D/(delz^2);
cdot(1, 1) = 4*Dr*(c(2, 1) - c(1, 1)) + 2*Dz*(c(1, 2) - c(1, 1));
for i=1:N-1
    cdot(i+1, 1) = (Dr/(2*i))*((2*i + 1)*c(i + 2, 1) + (2*i - 1)*c(i, 1) - 4*i*c(i + 1, 1))...
            + 2*Dz*(c(i + 1, 2) - c(i + 1, 1));
    for j=1:M-1
    cdot(1, j + 1) = 4*Dr*(c(2, j + 1) - c(1, j + 1)) + Dz*(c(1, j + 2) + c(1, j) - 2*c(1, j + 1));
    cdot(i + 1,j + 1) = (Dr/(2*i))*((2*i + 1)*c(i + 2, j + 1) + (2*i - 1)*c(i, j + 1) - 4*i*c(i + 1, j + 1))...
            + Dz*(c(i + 1, j + 2) + c(i + 1, j) - 2*c(i + 1, j + 1));
    end
end

xdot = [reshape(cdot, [], 1)];
end

