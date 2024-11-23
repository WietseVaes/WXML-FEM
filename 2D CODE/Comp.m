
%% Compute your solution u by first building you M,P,S and b matrices and then solving the system.
% [M,P,S,b, Mtilde, Ptilde, Stilde, hvec] = Build(x,y,elmat,elmatbd,f,g, vx, vy, Id, h, t);



% Identify interior nodes (nodes not in Id)
interior_nodes = setdiff(1:length(x), Id);

% Initialize uu for interior nodes and time steps
uu = zeros(length(interior_nodes), nt);

% Initial Dirichlet boundary condition
uu(:, 1) = h(interior_nodes, 1);

[M, P, S, b, Mtilde, Ptilde, Stilde] = build_tind(x, y, elmat, elmatbd, f(:, (1)), g(:, (1)), vx(:, (1)), vy(:, (1)), Id);


% Time-stepping loop
for i70 = 2:nt
    
    % Build the matrices for the current time step
    [P, b, Ptilde] = Build(x, y, elmat, elmatbd, f(:, (i70)), g(:, (i70)), vx(:, (i70)), vy(:, (i70)), Id);
 

    % Dirichlet boundary conditions for previous and current time steps
    hn = h(Id, i70-1);
    hnp1 = h(Id, i70);
    
    % Update solution for interior points at the current time step
    uu(:, i70) = (M - Dt * (P - D * S)) \ (Dt * b + M * uu(:, i70 - 1) + Mtilde * hn - ((Mtilde - Dt * (Ptilde - D * Stilde)) * hnp1));
end

u = zeros(length(x), nt);

for i70=1:nt
    u(Id, i70) = h(Id, i70);               % Dirichlet b  oundary condition at the current step
    u(interior_nodes, i70) = uu(:, i70); % Interior solution for current time step
end