% Generate P over a considered element i1 through
% theory and numerical methods. you need to output Pelem
% which is a 3x3 matrix in this case

function [Pelem] = GenerateElementMatrix(indices,x_total,y_total,vx,vy)

elemArea = get_area(indices,x_total,y_total);

xc = x_total(indices);
yc = y_total(indices);

A = zeros(length(xc));

for i=1:length(xc)
    A(i,:) = [1,xc(i),yc(i)];
end

X = A\eye(size(A));

Pelem = zeros(length(xc));

for i=1:length(xc)
    for j=1:length(yc)
        Pelem(i,j) = (elemArea / 3)*((X(2,i) * vx(indices(j)) )+(X(3,i) * vy(indices(j)) ));
    end
end


function area = get_area(indices, x_total, y_total)
xc = x_total(indices);
yc = y_total(indices);
area  = 0.5*abs(det([ones(3,1), xc, yc]));
end

end

