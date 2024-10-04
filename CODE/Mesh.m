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



end