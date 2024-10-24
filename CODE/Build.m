function [M,P,S,b] = Build(x,y,elmat,elmatbd)

%took out f and g from input...
n = length(x);

M       = sparse(n,n); %Mass matrix

P       = sparse(n,n); %Velocity matrix

S 		= sparse(n,n); %Stiffness matrix

b 		= zeros(n,1); %Right hand side


for i1 = 1:length(elmat(:,1)) % for all internal elements
    [Melem, Pelem, Selem] = GenerateElementMatrix(elmat(i1,:),x,y); %Generate small matrices (Melem, Pelem, Selem)
    elmat(i1,:)

    
    for i = 1:length(elmat(i1,:))
        for j = 1:length(elmat(i1,:))
            M(elmat(i1,i),elmat(i1,j)) = M(elmat(i1,i),elmat(i1,j)) + Melem(i,j);
            S(elmat(i1,i),elmat(i1,j)) = S(elmat(i1,i),elmat(i1,j)) + Selem(i,j);
        end
    end
    %M(i11,i11) = M(i11,i11) + MelemDiag;
    %M(i12,i12) = M(i12,i12) + MelemDiag;
    %M(i13,i13) = M(i13,i13) + MelemDiag;
  
    %Add small matrices to big matrices M,P and S

    GenerateElementVector; %calculating b vector over current element (belem)
    
    % add belem to big b


end


for i1 = 1:length(elmatbd(:,1)) % for all boundary elements 
    GenerateBoundaryElementVector; %calculating b vector over current neumann boundary element (belem)
    
    % add belem to big b

end

