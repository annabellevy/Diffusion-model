R = 1;  % Radius
h = 1/2; % Height
N_0 = 3; % Initial N
M_0 = 3; % Initial M
h_0 = 1 / (N_0 * M_0); % Initial mesh size
D = 1; % Diffusivity. 
T = 1;  % Length of time interval. 
tarr = [0:T/100:T];  % Array of time points in the interval 0 to T.
                     % The default is equally spaced 101 points including 0 and T.
Qt_acc = []; % Qt accumulator
             % Qt_acc(i, j) = Fractional release at time t_j for k = i
             % i.e. each row is a full run of the model with mesh size h_i
num_sizes = 12; % Number of different mesh sizes
alpha = 2;
harr = zeros(1, num_sizes); % This will be an array of mesh sizes h_k raised to alpha

% Should we have h_k = (1/k)*h_0 or h_k = (1/N*M) = (1/k^2)*h_0 ??
% Below uses the latter

for k=1:num_sizes

    % Update mesh size
    N = k * N_0;
    M = k * M_0;
    harr(1, k) = (1 / (N*M))^alpha;

    % Weights for Simpson's rule
    weights = zeros(M, N);
    n = N/2;
    m = M/2;
    for i=1:n
        weights(2*i, 1) = 4 * (2*i - 1);
        weights(2*i + 1, 1) = 2 * (2*i);
        for j=1:m
            weights(2*i, 2*j) = 4 * 4 * (2*i - 1);
            weights(2*i, 2*j + 1) = 4 * 2 * (2*i - 1);
            weights(2*i + 1, 2*j) = 2 * 4 * (2*i);
            weights(2*i + 1, 2*j + 1) = 2 * 2 * (2*i);
        end
    end

    % Initial concentration is constant inside. 
    S = 0;
    
    for i=1:N
        for j=1:M
            S = S + weights(i, j) * 1;
        end
    end
    S = 4*pi*R^2*h/(9*N^2*M) * S;
    C0=1/S;
    
    c0=C0*ones(N,M);
    c0 = vertcat(c0, zeros(1, N));
    c0 = horzcat(c0, zeros(M+1, 1));
    c0 = [reshape(c0, (N+1)*(M+1), 1)];
    
    % Solve the ODE 
    [tt,xt] = ode23s(@cylinder_diffv2,tarr,c0,[],D,R,h,N);
    xt=xt';
    
    % xt(iN + j, t) contains concentration
    % at radius ri, height zj, time t
    % where ri = iR/N and zj = jh/M
    % where i = 0,..., N - 1
    % and j = 0,..., M - 1
    xt = reshape(xt, (N+1), (M+1), length(tarr));
    
    % Fraction of quantity released.
    % Q = -2 \int_0^h \int_0^R c(r, z, t) 2 pi r dr dz
    % Use Simpson's rule.Assumes N, M are even.
    % Let n=N/2, m=M/2.
    % For each t, Q is the fraction remaining in the cylinder
    % so Qt(t) should be 1 - Q
    n = N/2;
    m = M/2;
    
    Qt = zeros(length(tarr), 1);
    for t=1:length(tarr)
        Q = 0;
        for i=1:N
            for j=1:M
                Q = Q + weights(i, j) * xt(i, j, t);
            end
        end
        Q = 4*pi*R^2*h/(9*N^2*M) * Q;
        Qt(t, 1) = 1 - Q;
    end
    Qt_acc = [Qt_acc, Qt(:)];
    
end

% Plot Q_t(h_k) against h_k^alpha for each t
for t = 1:length(tarr)
    hold on;
    plot(harr, Qt_acc(t, :));
end