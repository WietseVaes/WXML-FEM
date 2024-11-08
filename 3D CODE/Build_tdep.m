function [P, b, Ptilde] = Build_tdep(x, y, z, elmat, elmatbd, Bmat, f, g, vx, vy, vz, Id)

%took out f and g from input...
n       = length(x);

P       = sparse(n,n); %Velocity matrix

b 		= zeros(n,1); %Right hand side


for i1 = 1:length(elmat(:,1)) % for all internal elements
    Pelem = GenerateElementMatrix_tdep(elmat(i1,:), x, y, z, Bmat{i1}, vx, vy, vz);

    P(elmat(i1,:), elmat(i1,:)) = P(elmat(i1,:), elmat(i1,:)) + Pelem;

    belem = GenerateElementVector(elmat(i1,:),x,y,z,f);

    b(elmat(i1,:)) = b(elmat(i1,:)) + belem;

end


for i1 = 1:length(elmatbd(:,1)) % for all boundary elements
    belem = GenerateBoundaryElementVector(elmatbd(i1,:), x, y, z, g); %calculating b vector over current neumann boundary element (belem)

    b(elmatbd(i1,:)) = b(elmatbd(i1,:)) + belem;

end

Ptilde = P(:,Id);
Ptilde(Id,:) = [];
P(Id,:) = [];
P(:,Id) = [];

b(Id) = [];