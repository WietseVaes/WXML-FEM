function belem = GenerateBoundaryElementVector(indices, x, y, z, g)

Md = 4;

% Get the coordinates of the boundary element vertices
xc = x(indices);
yc = y(indices);
zc = z(indices);

if all(xc(1) == xc)
    AREA = abs(det([ones(Md-1, 1), yc, zc]));
elseif all(yc(1) == yc)
    AREA = abs(det([ones(Md-1, 1), xc, zc]));
elseif all(zc(1) == zc)
    AREA = abs(det([ones(Md-1, 1), xc, yc]));
end

% Using quadrature rule, compute belem
belem = AREA / factorial(Md-1) * g(indices);

end
