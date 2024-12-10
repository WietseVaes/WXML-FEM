function [M,  S, Mtilde, Stilde] = build_tind(x,y,elmat, Id)

n = length(x);

M       = sparse(n,n); %Mass matrix

S 		= sparse(n,n); %Stiffness matrix



for i1 = 1:length(elmat(:,1)) % for all internal elements
    [Melem, Selem] = GenerateElementMatrix_tind(elmat(i1,:),x,y);
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