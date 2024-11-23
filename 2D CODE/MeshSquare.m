%% Function

function [x,y,elmat,elmatbd, Id, In] = getMesh(dom_range,n)
% Computes the mesh of the FEM
%
% Input: dom_range; range of domain, for example: {[-1,1],[-2,2]}, i.e., x
%                   has range -1 to 1 and y has range -2 to 2
%        n; discretization in 1 dimension.
%
% Output: x; nx1 vector of x components of the mesh
%         y; nx1 vector of y components of the mesh  
%         elmat; element matrix that contains the indices of x and y
%           components that compose a triangle in the rows.
%         elmatbd; element matrix that contains the indices of x and y
%           components that compose the boundary in the rows.
%         Id: contains indices i1 of x and y for which (x(i1),y(i1)) lies
%           on the Dirichlet boundary.
%         In: contains indices i1 of x and y for which (x(i1),y(i1)) does
%           not lie on the Dirichlet boundary.

%USE DELAUNEY TO DEFINE ELMAT
    x = linspace(dom_range{1}(1), dom_range{1}(2), n); %x coordinates
    y = linspace(dom_range{2}(1), dom_range{2}(2), n); %y coordinates
    [X,Y] = meshgrid(x,y);
    X = X(:);
    Y = Y(:);
    elmat = delaunay(X,Y);

elmatbd = [];
for col = 1:(n-1)
    elmatbd = [elmatbd; (col - 1)*n + 1, (col)*n + 1]; %bottom
    elmatbd = [elmatbd; col*n, (col + 1)*n]; %top
end

for row = 1:(n-1)
    elmatbd = [elmatbd; row, row + 1]; %left
    elmatbd = [elmatbd; n*(n-1) + row, n*(n-1) + row + 1]; %right
end

Id = [];

for i = 1:n
    Id = [Id; i, 1]; %bottom
    Id = [Id; i, n]; %top
    Id = [Id; 1, i]; %left
    Id = [Id; n, i]; %right
end

Id = unique(Id,'rows');

In = [];

for row = 2:(n-1)
    for col = 2:(n-1)
        In = [In; row,col];
    end
end


figure;
trimesh(elmat,X,Y, 'Color', [1,0, 1]);
hold on;
plot(x(Id(:,1)), y(Id(:,2)), 'ro', 'MarkerFaceColor', 'r');
plot(x(In(:,1)), y(In(:,2)), 'ro', 'MarkerFaceColor', 'g');
hold off;



end
dom_range = {[-2,2], [-2,2]};
n = 5;
[x,y,elmat,elmatbd, Id, In] = getMesh(dom_range,n)

%% Testing

n = 10;
dom_range = {[-2,2], [-2,2]};

x = linspace(dom_range{1}(1), dom_range{1}(2), n); %x coordinates
y = linspace(dom_range{2}(1), dom_range{2}(2), n); %y coordinates
[X,Y] = meshgrid(x,y);
X = X(:)
Y = Y(:)
elmat = zeros((n-1)*(n-1)*2, 3);
index = 1;

for col = 1:(n-1)
    for row = 1:(n-1)
        bl = (col-1)*n + row;
        tl = (col-1)*n + row + 1;
        br = col*n + row;
        tr = col*n + row + 1;


        bottomTri = [bl, tl, br];
        topTri = [br, tr, tl];

        elmat(index,:) = bottomTri;
        elmat(index + 1, :) = topTri;
        

        index = index + 2;
        end
    end

elmatbd = [];

for col = 1:(n-1)
    elmatbd = [elmatbd; (col - 1)*n + 1, (col)*n + 1]; %bottm
    elmatbd = [elmatbd; col*n, (col + 1)*n]; %top
end

for row = 1:(n-1)
    elmatbd = [elmatbd; row, row + 1]; %left
    elmatbd = [elmatbd; n*(n-1) + row, n*(n-1) + row + 1]; %right
end


Id = [];

for i = 1:n
    Id = [Id; i, 1]; %bottom
    Id = [Id; i, n]; %top
    Id = [Id; 1, i]; %left
    Id = [Id; n, i]; %right
end

Id = unique(Id,'rows')

In = [];

for row = 2:(n-1)
    for col = 2:(n-1)
        In = [In; row,col];
    end
end

In

%gridCords = []

%for i = 1:length(X)
    %gridCords = [gridCords; X(i), Y(i)];
%end


figure;
trimesh(elmat,X,Y, 'Color', [1,0, 1]);
hold on;
plot(x(Id(:,1)), y(Id(:,2)), 'ro', 'MarkerFaceColor', 'r');
plot(x(In(:,1)), y(In(:,2)), 'ro', 'MarkerFaceColor', 'g');
hold off;