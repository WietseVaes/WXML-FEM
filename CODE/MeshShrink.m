function [x,y,elmat,elmatbd, Id, In] = Mesh(dom_range,n)
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
bnd_type = {'C Curve',3,2.5,1,1} ;
fattener = -0.2;
[s, x, y, keep] = setup_bd(bnd_type,n, dom_range, fattener);

%ADD BOUNDARY POINTS BEFORE DRAWING DELAUNAY
%BD normal
line_bd_pts = s.Zorigional(s.t);
line_bd_x = real(line_bd_pts).';
line_bd_y = imag(line_bd_pts).';

%BD skinny
bd_pts_skinny = s.Z(s.t);
bd_x_skinny = real(bd_pts_skinny).';
bd_y_skinny = imag(bd_pts_skinny).';

elmatbd = [line_bd_x,line_bd_y];
x_total = [x; line_bd_x];
y_total = [y; line_bd_y];

DB_tf = DirichletBD(line_bd_x, line_bd_y, s);

dirichlet_x = line_bd_x(DB_tf == 1);
dirichlet_y = line_bd_y(DB_tf == 1);

neuman_x = line_bd_x(DB_tf == 0);
neuman_y = line_bd_y(DB_tf == 0);



figure; 
plot(line_bd_x, line_bd_y, 'b-', 'LineWidth',2); hold on
scatter(dirichlet_x, dirichlet_y, 'o', 'MarkerFaceColor', [1, 0, 0], ...
    'MarkerEdgeColor', [1, 0, 0],'DisplayName', 'Dirichlet Boundary'); hold on%  red color
scatter(neuman_x, neuman_y, 'o', 'MarkerFaceColor', [0, 1, 0], ...
    'MarkerEdgeColor', [0, 1, 0], 'DisplayName', 'Neumann Boundary'); hold on %  green 
elmat_total = delaunay(x_total,y_total);
elmat_total = correct_elmat_vert(elmat_total, x_total, y_total,line_bd_x,line_bd_y);%Apply both filters....?
triplot(elmat_total,x_total,y_total); hold on
scatter(x_total,y_total,'*')
title('elamat with boundary pts');


%ORIGIONAL WITH OUT PTS ON BOUNDARY 

figure;
plot(bd_x_skinny, bd_y_skinny, 'b-', 'LineWidth',2); hold on
elmat = delaunay(x,y);
elmat = correct_elmat(elmat, x, y);
triplot(elmat,x,y); hold on
scatter(x,y,'*')
title('Origional elmat');

end

%ELIMINATING BASED ON AREA --> DOESNT GET ALL CASES 

function elmat = correct_elmat(elmat, x, y)
Area = zeros(size(elmat,1),1);

for i1 = 1:size(elmat,1)
    Index = elmat(i1,:);

    xc = x(Index);
    yc = y(Index);

    Area(i1) = abs(det([ones(length(xc),1) xc yc]));
     
end
elmat(Area > mean(Area),:) = [];

end

%ELIMINATING BASED ON VERTICES --> DOESNT GET ALL CASES

function elmat = correct_elmat_vert(elmat, x_total, y_total, bd_x, bd_y)
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