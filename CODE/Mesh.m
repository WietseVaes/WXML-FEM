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
bnd_type = {'Circles',[0 + 0i],[2]};
fattener = 0;
[s, x, y, keep] = setup_bd(bnd_type,n, dom_range, fattener);
elmat = delaunay(x,y);

elmat = correct_elmat(elmat, x, y);

triplot(elmat,x,y); hold on
%triplot(elmatbd,x,y,'k');
scatter(x,y,'*')


end

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
