function Pelem = GenerateElementMatrix_tdep(indices,x,y,z,Bmat,vx,vy,vz)

Md = 4;

xc = x(indices);
yc = y(indices);
zc = z(indices);

%% Melem
elemArea = abs(det([ones(Md, 1), xc, yc, zc]));

Pelem = zeros(Md);

for i=1:length(xc)
    for j=1:length(yc)
        Pelem(i,j) = (elemArea / factorial(Md))*((Bmat(2,i) * vx(indices(j)))+(Bmat(3,i) * vy(indices(j)))+(Bmat(4,i) * vz(indices(j))));
    end
end

end