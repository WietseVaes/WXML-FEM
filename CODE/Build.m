function [M,P,S,b] = Build(x,y,elmat,f,elmatbd,g)

n = size(x);

M 		= sparse(n,n); %Mass matrix

P       = sparse(n,n); %Velocity matrix

S 		= sparse(n,n); %Stiffness matrix

b 		= zeros(n,1); %Right hand side


for i1 = 1:length(elmat(:,1)) % for all internal elements
    GenerateElementMatrix; %Generate small matrices (Melem, Pelem, Selem)

    %Add small matrices to big matrices M,P and S

    GenerateElementVector; %calculating b vector over current element (belem)
    
    % add belem to big b


end


for i1 = 1:length(elmatbd(:,1)) % for all boundary elements 
    GenerateBoundaryElementVector; %calculating b vector over current neumann boundary element (belem)
    
    % add belem to big b

end

