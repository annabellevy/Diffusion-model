% Run Cylinder Diffusion.
R = 1;  % Radius
h = 1;
N = 50;
M = 50;
D = 1; % Diffusivity. 
T = 1;  % Length of time interval. 
tarr = [0:T/100:T];  % Array of time points in the interval 0 to T.
                     % The default is equally spaced 101 points including 0 and T.

% Initial concentration is constant inside. 
S = 0;
n = N/2;
m=M/2;
for i=1:n
    S = S + (4 * (2*i - 1) + 2 * (2*i)) * (1 + 4 * m + 2 * (m - 1));
end
S = 4*pi*R^2*h/(9*N^2*M) * S;
C0=1/S;
c0=C0*ones(N,M);
c0 = vertcat(c0, zeros(1, N));
c0 = horzcat(c0, zeros(M+1, 1));
c0 = [reshape(c0, (N+1)*(M+1), 1);0];

% Solve the ODE 
%%[tt,xt] = ode23s(@cylinder_diff,tarr,c0,[],D,R,h,N);
odeset('JPattern',S);
[tt,xt] = ode23s(@cylinder_diff,tarr,c0,[],D,R,h,N);
tt=tt';xt=xt';
% Qt is the array of fraction released. 
Qt=xt(end,:);
plot(tarr, Qt);
% To plot, use plot(tarr,Qt)