function belem = GenerateBoundaryElementVector(indices, x, y, g)
% Generate the matrices of b over a considered boundary element i1 (use g). Use the
% quadrature rule we discussed. you should have the result be belem, a 2x1
% vector.
%
% Input:
%   x, y: Coordinates of the points
%   elmat: Element connectivity matrix for boundary elements (2 vertices per boundary element)
%   g: Function handle for boundary data (e.g., g(x, y))
%
% Output:
%   belem: Boundary element vector (2x1 for each boundary element)
belem = zeros(2, 1);

% Get the coordinates of the boundary element vertices
xe = x(indices);
ye = y(indices);
Length = sqrt((xe(2) - xe(1))^2 + (ye(2) - ye(1))^2);

% Evaluate the function g at the midpoint of the boundary element
% g_mid = g((xe(1) + xe(2)) / 2, (ye(1) + ye(2)) / 2);

% Using quadrature rule, compute belem
belem = (Length / 2) * g(xe,ye);

end
