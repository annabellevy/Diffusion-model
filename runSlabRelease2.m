% Run Slab Diffusion.
h = 1/2;  % h is half thickness of the slab. 
N = 100;  % Number of grid points. 500 is good. N must be even.
D = 1; % Diffusivity.
T = 1; % Length of time interval.
tarr = [0:T/100:T];  % Array of time points in the interval 0 to T.
                     % The default is equally spaced 101 points including 0 and T.


% Weights for Simpson's rule
weights = zeros(N, 1);
weights(1, 1) = 2;
for i=1:N/2
    weights(2*i, 1) = 2;
    weights(2*i + 1) = 4;
end
weights(N, 1) = 1;

% Initial concentration is constant inside. 
% Vector c0 = c0(1),c0(2),..,c0(N) 
% Initial concentration at locations 
% 0/N, h/N, 2h/N,..,h(N-1)/N. 
% We want Q0=\int_0^h c0(r)dr=1. 
% Using Simpsons rule.
S = 0;
for i=1:N
    S = S + weights(i, 1) * 1;
end
C0 = 2*S*h/(N*3);
c0=C0*ones(N,1);

% Solve the ODE
[tt,xt] = ode23s(@slab_diff,tarr,c0,[],D,h);
tt=tt'; xt=xt';

Qt = zeros(length(tarr), 1);
for t=1:length(tarr)
    Q = 0;
    for i=1:N
        Q = Q + weights(i) * xt(i, t);
    end
    Q = 2*Q*h/(3*N);
    Qt(t) = 1 - Q;
end

Qrt = repmat(Qt', 11, 1);
hold on;
Qrt = Qrt';
rarr = ones(11, 1)';

surf(rarr, tarr, Qrt, 'EdgeColor', 'red', 'LineWidth', 2)