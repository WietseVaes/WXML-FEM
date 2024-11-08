function [X, Y, Z, elmat, elmatbd, Id] = Mesh_rectangle(dom_range, n)


x = linspace(dom_range{1}(1),dom_range{1}(2),n);
y = linspace(dom_range{2}(1),dom_range{2}(2),n);
z = linspace(dom_range{3}(1),dom_range{3}(2),n);

X = zeros(n^3,1);
Y = zeros(n^3,1);
Z = zeros(n^3,1);

j = 1;
for i1 = 1:length(x)
    for i2 = 1:length(y)
        for i3 = 1:length(z)
            X(j) = x(i1);
            Y(j) = y(i2);
            Z(j) = z(i3);
            j = j + 1;
        end
    end
end

DT = delaunayTriangulation(X,Y,Z);

% tetramesh(DT,'FaceAlpha',0.3);
elmat = DT.ConnectivityList;
[elmatbd, ~] = convexHull(DT);

Ibd = unique(elmatbd)';

%% Manual definition Dirichlet Boundary
Id = [];
% Case: x = dom_range{1}(1), or
%       y = dom_range{2}(2) & dom_range{3}(1)+(dom_range{3}(2)-dom_range{3}(1))/3 < z < dom_range{3}(1)+2*(dom_range{3}(2)-dom_range{3}(1))/3
for i1 = Ibd
    if abs(X(i1) - dom_range{1}(1)) <= eps 
        Id = [Id,i1];
    end
    if abs(Y(i1) - dom_range{2}(2)) <= eps && dom_range{3}(1)+(dom_range{3}(2)-dom_range{3}(1))/3 < Z(i1) && Z(i1) < dom_range{3}(1)+2*(dom_range{3}(2)-dom_range{3}(1))/3
        Id = [Id,i1];
    end
end






