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
c0 = [reshape(c0, (N+1)*(M+1), 1)];

% Solve the ODE 
%%[tt,xt] = ode23s(@cylinder_diff,tarr,c0,[],D,R,h,N);
odeset('JPattern',S);
[tt,xt] = ode23s(@cylinder_diffv2,tarr,c0,[],D,R,h,N);
xt=xt';

% xt(iN + j, t) contains concentration
% at radius ri, height zj, time t
% where ri = iR/N and zj = jh/M
% where i = 0,..., N - 1
% and j = 0,..., M - 1
xt = reshape(xt, (N+1), (M+1), length(tarr));

% Fraction of quantity released.
% Qdot = -2 \int_0^h \int_0^R c(r, z, t) 2 pi r dr dz
% Use Simpson's rule.Assumes N, M are even.
% Let n=N/2, m=M/2.
% For each t, Q is the fraction remaining in the cylinder
% so Qt(t) should be 1 - Q
Qt = zeros(length(tarr) - 1, 1); % Array of fraction released
Q = 0;
n = N/2;
m=M/2;
for t=1:length(tarr)
   Q = 0;
    for i=1:n
        Q = Q + 4 * (2*i - 1) * xt(2*i, 1, t) + 2 * (2*i) * xt(2*i + 1, 1, t);
        for j=1:m
            Q = Q + 4 * (2*i - 1) * (4 * xt(2*i, 2*j, t) + 2 * xt(2*i, 2*j + 1, t))...
                + 2 * (2*i) * (4 * xt(2*i + 1, 2*j, t) + 2 * xt(2*i + 1, 2*j + 1, t));
        end
        Qt(t, 1) = 1 - (4*pi*R^2*h/(9*N^2*M) * Q);
    end
end

plot(tarr, Qt);
