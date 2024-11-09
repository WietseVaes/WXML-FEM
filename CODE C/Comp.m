
%% Compute your solution u by first building you M,P,S and b matrices and then solving the system.
% [M,P,S,b, Mtilde, Ptilde, Stilde, hvec] = Build(x,y,elmat,elmatbd,f,g, vx, vy, Id, h, t);

u = zeros(length(x), nt);

% Identify interior nodes (nodes not in Id)
interior_nodes = setdiff(1:length(x), Id);

% Initialize uu for interior nodes and time steps
uu = zeros(length(interior_nodes), nt);

% Initial Dirichlet boundary condition
u(Id, 1) = h(x(Id), y(Id), t(1));

% Time-stepping loop
for i70 = 2:nt
    % Current and previous time values
    current_t = t(i70);
    previous_t = t(i70 - 1);
    
    % Build the matrices for the current time step
    [M, P, S, b, Mtilde, Ptilde, Stilde] = Build(x, y, elmat, elmatbd, ...
                                                    @(x, y) f(x, y, current_t), ...
                                                    @(x, y) g(x, y, current_t), ...
                                                    @(x, y) vx(x, y, current_t), ...
                                                    @(x, y) vy(x, y, current_t), ...
                                                    Id);

    % Dirichlet boundary conditions for previous and current time steps
    hn = h(x(Id), y(Id), previous_t);
    hnp1 = h(x(Id), y(Id), current_t);
    
    % Update solution for interior points at the current time step
    uu(:, i70) = (M - Dt * (P - D * S)) \ (Dt * b + M * uu(:, i70 - 1) + Mtilde * hn - (Mtilde - Dt * (Ptilde - D * Stilde)) * hnp1);
    
    % Update the full solution vector with boundary conditions
    u(Id, i70) = hnp1;               % Dirichlet b  oundary condition at the current step
    u(setdiff(1:length(x), Id), i70) = uu(:, i70); % Interior solution for current time step

end