
%% Compute your solution u by first building you M,P,S and b matrices and then solving the system.
% [M,P,S,b, Mtilde, Ptilde, Stilde, hvec] = Build(x,y,elmat,elmatbd,f,g, vx, vy, Id, h, t);



% Identify interior nodes (nodes not in Id)
interior_nodes = setdiff(1:length(x), Id);

% Initialize uu for interior nodes and time steps
uu = zeros(length(interior_nodes), nt);

% Initial Dirichlet boundary condition
uu(:, 1) = h(interior_nodes, 1);

[M, S, Mtilde,  Stilde] = build_tind(x, y, elmat, elmatbd, f(:, (1)), g(:, (1)), vx(:, (1)), vy(:, (1)), Id);

u = zeros(length(x), nt);
u(Id, 1) = h(Id, 1);               % Dirichlet b  oundary condition at the current step
u(interior_nodes, 1) = uu(:, 1); % Interior solution for current time step

% Time-stepping loop
if Time_Meth == "RK4"
    for i70 = 1:nt-1
        % Build the matrices for the current time step
        [Pn, bn, Ptilden] = Build(x, y, elmat, elmatbd, f(:, i70), g(:, i70), vx(:, i70), vy(:, i70), Id);
        [P2, b2, Ptilde2] = Build(x, y, elmat, elmatbd, f2(:, i70), g2(:, i70), vx2(:, i70), vy2(:, i70), Id);
        [Pnp1, bnp1, Ptildenp1] = Build(x, y, elmat, elmatbd, f(:, i70+1), g(:, i70+1), vx(:, i70+1), vy(:, i70+1), Id);
        % Dirichlet boundary conditions for previous and current time steps
        hn = h(Id, i70);
        hn2 = h2(Id, i70);
        hnp1 = h(Id, i70+1);
        Dthn = Dth(Id, i70);
        Dthn2 = Dth2(Id, i70);
        Dthnp1 = Dth(Id, i70+1);

        k1 = M \ (bn + ((Pn - D * S) * uu(:, i70) + (Ptilden - D * Stilde) * hn - Mtilde * Dthn));

        k2 = M \ (b2 + (P2 - D * S) * (uu(:, i70) + Dt/2 * k1)  + (Ptilde2 - D * Stilde) * hn2 - Mtilde * Dthn2);

        k3 = M \ (b2 + (P2 - D * S) * (uu(:, i70) + Dt/2 * k2)  + (Ptilde2 - D * Stilde) * hn2 - Mtilde * Dthn2);

        k4 = M \ (bnp1  + (Pnp1 - D * S) * (uu(:, i70) + Dt * k3) + (Ptildenp1 - D * Stilde) * hnp1- Mtilde * Dthnp1);

        % Update solution for interior points at the current time step
        uu(:, i70+1) = uu(:, i70) + (Dt / 6) * (k1 + 2 * k2 + 2 * k3 + k4);

        u(Id, i70+1) = h(Id, i70+1);               % Dirichlet b  oundary condition at the current step
        u(interior_nodes, i70+1) = uu(:, i70+1);   % Interior solution for current time step
    end
elseif Time_Meth == "CN"
    for i70 = 1:nt-1
        % Build the matrices for the current time step
        [Pn, bn, Ptilden] = Build(x, y, elmat, elmatbd, f(:, (i70)), g(:, (i70)), vx(:, (i70)), vy(:, (i70)), Id);
        [Pnp1, bnp1, Ptildenp1] = Build(x, y, elmat, elmatbd, f(:, (i70+1)), g(:, (i70+1)), vx(:, (i70+1)), vy(:, (i70+1)), Id);


        % Dirichlet boundary conditions for previous and current time steps
        hn = h(Id, i70);
        hnp1 = h(Id, i70+1);
        Dthn = Dth(Id, i70);
        Dthnp1 = Dth(Id, i70+1);

        % Update solution for interior points at the current time step
        RHS = M*uu(:, i70) + Dt/2*(bn+bnp1) + Dt/2*((Ptildenp1-D*Stilde)*hnp1-Mtilde*Dthnp1 + (Ptilden-D*Stilde)*hn-Mtilde*Dthn + (Pn-D*S)*uu(:,i70));
        uu(:, i70+1) = (M - Dt/2 * (Pnp1 - D * S)) \ RHS;


        u(Id, i70+1) = h(Id, i70+1);               % Dirichlet b  oundary condition at the current step
        u(interior_nodes, i70+1) = uu(:, i70+1); % Interior solution for current time step
    end
elseif Time_Meth == "BDF2"
    uu(:, 2) = h(interior_nodes, 2);
    u(Id, 2) = h(Id, 2);               % Dirichlet b  oundary condition at the current step
    u(interior_nodes, 2) = uu(:, 2);
    for i70 = 3:nt
        % Build the matrices for the current time step
        [P, b, Ptilde] = Build(x, y, elmat, elmatbd, f(:, (i70)), g(:, (i70)), vx(:, (i70)), vy(:, (i70)), Id);


        % Dirichlet boundary conditions for previous and current time steps
      
        hn = h(Id, i70);
        Dthn = Dth(Id, i70);

        % Update solution for interior points at the current time step
        uu(:, i70) = (M - 2/3*Dt * (P - D * S)) \ (2 / 3 * Dt * b + 4/3*M * uu(:, i70 - 1) - 1/3 * M * uu(:, i70 - 2) +  2 / 3*Dt * ((Ptilde - D * Stilde) * hn - Mtilde * Dthn));


        u(Id, i70) = h(Id, i70);               % Dirichlet b  oundary condition at the current step
        u(interior_nodes, i70) = uu(:, i70); % Interior solution for current time step
    end
else
    for i70 = 2:nt
        % Build the matrices for the current time step
        [P, b, Ptilde] = Build(x, y, elmat, elmatbd, f(:, (i70)), g(:, (i70)), vx(:, (i70)), vy(:, (i70)), Id);


        % Dirichlet boundary conditions for previous and current time steps
        hn = h(Id, i70-1);
        hnp1 = h(Id, i70);

        % Update solution for interior points at the current time step
        uu(:, i70) = (M - Dt * (P - D * S)) \ (Dt * b + M * uu(:, i70 - 1) + Mtilde * hn - ((Mtilde - Dt * (Ptilde - D * Stilde)) * hnp1));


        u(Id, i70) = h(Id, i70);               % Dirichlet b  oundary condition at the current step
        u(interior_nodes, i70) = uu(:, i70); % Interior solution for current time step
    end
end