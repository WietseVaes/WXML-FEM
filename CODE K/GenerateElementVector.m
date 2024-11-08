function belem = GenerateElementVector(indices,x, y, f)
% Generate the interior element vector belem using the quadrature rule.
% 
% Input:
%   x, y: Coordinates of the points
%   elmat: Element connectivity matrix for interior elements (3 vertices per triangular element)
%   f: Function handle for interior source term (e.g., f(x, y))
%
% Output:
%   belem: Interior element vector (3x1 for each triangular element)

belem = zeros(3, 1);
xe = x(indices);
ye = y(indices);
Area = abs(det([ones(3,1) xe ye])/2);

% Evaluate the function f at each vertex
f_values = f(indices);

for i = 1:3
    belem(i) = (Area / 3) * f_values(i);  % Only f(x_i) contributes for phi_i
end

end
