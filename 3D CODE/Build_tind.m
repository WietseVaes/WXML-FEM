function [M, S, Mtilde, Stilde, Bmat] = Build_tind(x, y, z, elmat, Id)

%took out f and g from input...
n = length(x);

M       = sparse(n,n); %Mass matrix

S 		= sparse(n,n); %Stiffness matrix


for i1 = 1:length(elmat(:,1)) % for all internal elements

    [Melem, Selem, X] = GenerateElementMatrix_tind(elmat(i1,:), x, y, z); %Generate small matrices (Melem, Pelem, Selem)
    
    Bmat{i1} = X;

    M(elmat(i1,:), elmat(i1,:)) = M(elmat(i1,:), elmat(i1,:)) + Melem;
    S(elmat(i1,:), elmat(i1,:)) = S(elmat(i1,:), elmat(i1,:)) + Selem;

end

Mtilde = M(:,Id);
Mtilde(Id,:) = [];
M(Id,:) = [];
M(:,Id) = [];

Stilde = S(:,Id);
Stilde(Id,:) = [];
S(Id,:) = [];
S(:,Id) = [];
