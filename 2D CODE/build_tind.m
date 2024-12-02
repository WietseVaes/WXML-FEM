function [M,  S, Mtilde, Stilde] = build_tind(x,y,elmat,elmatbd,f,g,vx,vy, Id)

%took out f and g from input...
n = length(x);

M       = sparse(n,n); %Mass matrix

P       = sparse(n,n); %Velocity matrix

S 		= sparse(n,n); %Stiffness matrix

b 		= zeros(n,1); %Right hand side


for i1 = 1:length(elmat(:,1)) % for all internal elements
    [Melem, ~, Selem] = GenerateElementMatrix(elmat(i1,:),x,y,vx,vy); %Generate small matrices (Melem, Pelem, Selem)


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