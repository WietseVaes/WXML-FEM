function belem = GenerateElementVector(indices, x, y, z, f)

Md = 4;

xc = x(indices);
yc = y(indices);
zc = z(indices);

Area = abs(det([ones(Md,1) xc yc zc]));

belem = (Area / factorial(Md)) * f(indices);  % Only f(x_i) contributes for phi_i

end
