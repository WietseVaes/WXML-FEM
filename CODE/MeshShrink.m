function [x_total,y_total,elmat,elmatbd, Id, In] = MeshShrink(bnd_type, dom_range,n, Dir_int,f,g,h)
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


%                   2D: {'Kite'},

%                       {'Star',a,w}: a = indent (0<a<1); w = #petals
%                       {'Star',.3,7}

%                       {'Circles',c,R}: c = centers in array; R =
%                       radii in array;
%                       {'Circles',[0 + 0i],[2]}

%                       {'C Curve',a,b,c,d}, a = sharpennes (bigger =
%                       sharper); b = cavity controller (smaller = less
%                       cavity); c = fatness (larger = fatter all around);
%                       d = direction (1 = cavity left; exp(1i*(pi-theta)) = rotation theta)  
%                       {'C Curve',3,2.5,1,1} 
%                       {'C Curve',3,2.5,1,exp(-1i*3*pi/4)}

%DUPLICATE DATA POINTS IN ELMAT WHEN USING C CURVE....
%bnd_type = {'C Curve',3,2.5,1,1} ;
shrinker = -0.001;
[s, x, y, keep, P] = setup_bd(bnd_type,n, dom_range, shrinker);


%BD EQUIDISTRIBUTED
bd = P;
bd_x = real(bd).';
bd_y = imag(bd).';

x_total = [x; bd_x];
y_total = [y; bd_y];

DB_tf = DirichletBD(bd_x, bd_y, s,Dir_int);
Id = find(DB_tf == 1) + size(x,1);
In = find(DB_tf == 0) + size(x,1);

[found, indices] = ismember(bd_x, x_total);
elmatbd = indices(found);
elmatbd = [elmatbd, circshift(elmatbd, -1)];

%EVEN PLOT
figure; 

plot(real(s.original), imag(s.original), 'r-', 'LineWidth',2); hold on
plot(real(s.Z), imag(s.Z), 'b-', 'LineWidth',2')
plot(P, 'ok', 'MarkerSize',8 )
scatter(x_total(Id), y_total(Id), 'o', 'MarkerFaceColor', [0, 1, 0], ...
    'MarkerEdgeColor', [0, 1, 0],'DisplayName', 'Dirichlet Boundary'); hold on%  red color

scatter(x_total(In), y_total(In), 'o', 'MarkerFaceColor', [1, 0, 0], ...
    'MarkerEdgeColor', [1, 0, 0], 'DisplayName', 'Neumann Boundary'); hold on %  green 

elmat = delaunay(x_total,y_total);
elmat = correct_elmat(elmat, x_total, y_total,bd_x,bd_y);
triplot(elmat,x_total,y_total); hold on
scatter(x_total,y_total,'*')
title('elamat with equidistant boundary pts');

end


%ELIMINATING BASED ON VERTICES

function elmat = correct_elmat(elmat, x_total, y_total, bd_x, bd_y)
keep_triangle = true(size(elmat, 1), 1);
    
 for ind = 1:size(elmat, 1)
    Index = elmat(ind, :);  
        
    triangle_vertices_x = x_total(Index);
    triangle_vertices_y = y_total(Index);
        
    % Check if all vertices of the triangle are boundary vertices
    if all(ismember(triangle_vertices_x, bd_x) & ismember(triangle_vertices_y, bd_y))
        keep_triangle(ind) = false; 
    end
end
    
    elmat = elmat(keep_triangle, :);
end