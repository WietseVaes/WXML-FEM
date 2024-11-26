
%% Compute your solution u by first building you M,P,S and b matrices and then solving the system.
% [M,P,S,b, Mtilde, Ptilde, Stilde, hvec] = Build(x,y,elmat,elmatbd,f,g, vx, vy, Id, h, t);

% Identify interior nodes (nodes not in Id)
interior_nodes = setdiff(1:length(x), Id);

% Initialize uu for interior nodes and time steps
uu = zeros(length(interior_nodes), nt);

% Initial Dirichlet boundary condition
uu(:, 1) = h(interior_nodes, 1);

[M, P, S, b, Mtilde, Ptilde, Stilde] = build_tind(x, y, elmat, elmatbd, f(:, (1)), g(:, (1)), vx(:, (1)), vy(:, (1)), Id);

i70 = 1;
u = zeros(length(x), nt);
u(Id, i70) = h(Id, i70);               % Dirichlet b  oundary condition at the current step
u(interior_nodes, i70) = uu(:, i70); % Interior solution for current time step

% Time-stepping loop
for i70 = 1:nt-1
    
    % Build the matrices for the current time step
    [Pn, bn, Ptilden] = Build(x, y, elmat, elmatbd, f(:, (i70)), g(:, (i70)), vx(:, (i70)), vy(:, (i70)), Id);
 

    % Dirichlet boundary conditions for previous and current time steps
    hn = h(Id, i70);
    hnp1 = h(Id, i70+1);
    % % Compute k1
    k1 = bn + (M - (Pn - D * S)) * uu(:, i70) + (Mtilde - (Ptilden - D * Stilde)) * hn;

    % Compute k2
    [P_half, b_half, Ptilde_half] = Build(x, y, elmat, elmatbd, ...
        f2(:,i70), ...
        g2(:, i70), ...
        vx2(:, i70), ...
        vy2(:, i70), ...
        Id);

    k2 = b_half + (M - (P_half - D * S)) * (uu(:, i70) + 0.5 * Dt * k1) + ...
         (Mtilde - (Ptilde_half - D * Stilde)) * hn;

    % Compute k3
    k3 = b_half + (M - (P_half - D * S)) * (uu(:, i70) + 0.5 * Dt * k2) + ...
         (Mtilde - (Ptilde_half - D * Stilde)) * hn;

    % Compute k4
    [Pnp1, bnp1, Ptildenp1] = Build(x, y, elmat, elmatbd, ...
        f(:, i70+1), g(:, i70+1), vx(:, i70+1), vy(:, i70+1), Id);
    k4 = bnp1 + (M - (Pnp1 - D * S)) * (uu(:, i70) + Dt * k3) + ...
         (Mtilde - (Ptildenp1 - D * Stilde)) * hn;


    % [P, b, Ptilde] = Build(x, y, elmat, elmatbd, f(:, i70), g(:, i70), vx(:, i70), vy(:, i70), Id);
    % hn = h(Id, i70); % Dirichlet boundary condition at current time step
    % k1 = (M - (P - D * S)) \ ...
    %      (b + M * uu(:, i70) + Mtilde * hn - (Mtilde - (Ptilde - D * Stilde)) * hn);
    % 
    % % Stage 2: Compute k2
    % [P_half, b_half, Ptilde_half] = Build(x, y, elmat, elmatbd, ...
    %     (f(:, i70) + f(:, i70+1)) / 2, ...
    %     (g(:, i70) + g(:, i70+1)) / 2, ...
    %     (vx(:, i70) + vx(:, i70+1)) / 2, ...
    %     (vy(:, i70) + vy(:, i70+1)) / 2, Id);
    % h_half = (h(Id, i70) + h(Id, i70+1)) / 2; % Boundary condition at t + Dt/2
    % k2 = (M - (P_half - D * S)) \ ...
    %      (b_half + M * (uu(:, i70) + 0.5 * Dt * k1) + Mtilde * h_half - ...
    %       (Mtilde - (Ptilde_half - D * Stilde)) * h_half);
    % 
    % % Stage 3: Compute k3
    % k3 = (M - (P_half - D * S)) \ ...
    %      (b_half + M * (uu(:, i70) + 0.5 * Dt * k2) + Mtilde * h_half - ...
    %       (Mtilde - (Ptilde_half - D * Stilde)) * h_half);
    % 
    % % Stage 4: Compute k4
    % [P_next, b_next, Ptilde_next] = Build(x, y, elmat, elmatbd, f(:, i70+1), g(:, i70+1), vx(:, i70+1), vy(:, i70+1), Id);
    % h_next = h(Id, i70+1); % Dirichlet boundary condition at t + Dt
    % k4 = (M - (P_next - D * S)) \ ...
    %      (b_next + M * (uu(:, i70) + Dt * k3) + Mtilde * h_next - ...
    %       (Mtilde - (Ptilde_next - D * Stilde)) * h_next);
    
    % Update solution using RK4 formula
    uu(:, i70+1) = uu(:, i70) + M \ (Mtilde*(hn-hnp1)) + M \ ((Dt / 6) * (k1 + 2 * k2 + 2 * k3 + k4));
    u(Id, i70+1) = hnp1;               % Dirichlet b  oundary condition at the current step
    u(interior_nodes, i70+1) = uu(:, i70+1); % Interior solution for current time step
    disp([num2str(i70/nt*100,3), '% done.'])
end
