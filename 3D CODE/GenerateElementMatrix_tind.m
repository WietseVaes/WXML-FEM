function [Melem, Selem, X] = GenerateElementMatrix_tind(indices, x, y, z)

Md = 4;

xc = x(indices);
yc = y(indices);
zc = z(indices);

%% Melem
elemArea = abs(det([ones(Md, 1), xc, yc, zc]));
Melem = ((elemArea)/factorial(Md+1))*(1 + eye(Md));


A = [ones(Md, 1), xc, yc, zc];

X = A \ eye(Md);

%% Selem

Selem = zeros(length(xc));

for i=1:length(xc)
    for j=1:length(yc)
        Selem(i,j) = (elemArea/factorial(Md-1))*((X(2,j) * X(2,i))+(X(3,j) * X(3,i))+(X(4,j) * X(4,i)));
    end
end


end

