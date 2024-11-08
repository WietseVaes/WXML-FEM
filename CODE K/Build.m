function [P, b, Ptilde] = Build(x,y,elmat,elmatbd,f,g,vx,vy, Id)

%took out f and g from input...
n = length(x);

P       = sparse(n,n); %Velocity matrix


b 		= zeros(n,1); %Right hand side


for i1 = 1:length(elmat(:,1)) % for all internal elements
    [Melem, Pelem, Selem] = GenerateElementMatrix(elmat(i1,:),x,y,vx,vy); %Generate small matrices (Melem, Pelem, Selem)

    P(elmat(i1,:), elmat(i1,:)) = P(elmat(i1,:), elmat(i1,:)) + Pelem;

    belem = GenerateElementVector(elmat(i1,:),x,y,f); %calculating b vector over current element (belem)


    b(elmat(i1,:)) = b(elmat(i1,:)) + belem;

    % add belem to big b


end


for i1 = 1:length(elmatbd(:,1)) % for all boundary elements
    belem = GenerateBoundaryElementVector(elmatbd(i1,:),x,y,g); %calculating b vector over current neumann boundary element (belem)

    b(elmatbd(i1,:)) = b(elmatbd(i1,:)) + belem;

    % add belem to big b

end

Ptilde = P(:,Id);
Ptilde(Id,:) = [];
P(Id,:) = [];
P(:,Id) = [];

b(Id) = [];


%hvec = h(x(Id),y(Id));
%size(hvec)