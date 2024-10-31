
%% Compute your solution u by first building you M,P,S and b matrices and then solving the system.
[M,P,S,b, Mtilde, Ptilde, Stilde, hvec] = Build(x,y,elmat,elmatbd,f,g, vx, vy, Id,h);

uu = (M- P + D*S)\(b-((Mtilde - Ptilde + D*Stilde)*hvec));

u = zeros(length(x),1);
u(Id) = hvec; %Puts the dirichlet boundary condition in the solution matrix
u(setdiff(1:length(x), Id)) = uu; %puts approximations in the solution matrix