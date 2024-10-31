function [M,P,S,b, Mtilde, Ptilde, Stilde, hvec] = Build(x,y,elmat,elmatbd,f,g,vx,vy, Id,h)

%took out f and g from input...
n = length(x);

M       = sparse(n,n); %Mass matrix

P       = sparse(n,n); %Velocity matrix

S 		= sparse(n,n); %Stiffness matrix

b 		= zeros(n,1); %Right hand side


for i1 = 1:length(elmat(:,1)) % for all internal elements
    [Melem, Pelem, Selem] = GenerateElementMatrix(elmat(i1,:),x,y,vx,vy); %Generate small matrices (Melem, Pelem, Selem)

    for i = 1:length(elmat(i1,:))
        for j = 1:length(elmat(i1,:))
            M(elmat(i1,i),elmat(i1,j)) = M(elmat(i1,i),elmat(i1,j)) + Melem(i,j);
            S(elmat(i1,i),elmat(i1,j)) = S(elmat(i1,i),elmat(i1,j)) + Selem(i,j);
            P(elmat(i1,i),elmat(i1,j)) = P(elmat(i1,i),elmat(i1,j)) + Pelem(i,j);

        end
    end

    belem = GenerateElementVector(elmat(i1,:),x,y,f); %calculating b vector over current element (belem)
    
    for i = 1:length(elmat(i1,:))
        b(elmat(i1,i)) = b(elmat(i1,i)) + belem(i);
           
    end
    
    % add belem to big b


end


for i1 = 1:length(elmatbd(:,1)) % for all boundary elements 
    belem = GenerateBoundaryElementVector(elmatbd(i1,:),x,y,g); %calculating b vector over current neumann boundary element (belem)
    
    for i = 1:length(elmatbd(i1,:))
        b(elmatbd(i1,i)) = b(elmatbd(i1,i)) + belem(i);
           
    end
    % add belem to big b

end

Mtilde = M(:,Id);
Mtilde(Id,:) = [];
M(Id,:) = [];
M(:,Id) = [];

Stilde = S(:,Id);
Stilde(Id,:) = [];
S(Id,:) = [];
S(:,Id) = [];

Ptilde = P(:,Id);
Ptilde(Id,:) = [];
P(Id,:) = [];
P(:,Id) = [];

b(Id) = [];


hvec = h(x(Id),y(Id));
size(hvec)