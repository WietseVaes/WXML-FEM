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
else
    v1 = [xc(2) - xc(1), yc(2) - yc(1), zc(2) - zc(1)];
    v2 = [xc(3) - xc(1), yc(3) - yc(1), zc(3) - zc(1)];

    AREA = 0.5*norm(cross(v1,v2));
end

% Using quadrature rule, compute belem
belem = AREA / factorial(Md-1) * g(indices);

end
