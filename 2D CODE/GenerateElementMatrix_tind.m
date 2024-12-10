% Generate the matrices of M and S over a considered element i1 through
% theory and numerical methods. you need to output Melem, Pelem, Selem
% which are 3x3 matrices in this case

function [Melem, Selem] = GenerateElementMatrix_tind(indices,x_total,y_total)

%Melem
Md = 3;
elemArea = get_area(indices,x_total,y_total);
Melem = ((2*elemArea)/factorial(Md+1))*(1 + eye(3));

%Selem

%Solve for alpha, beta, gamma

xc = x_total(indices);
yc = y_total(indices);

%phi1
A = zeros(length(xc));

for i=1:length(xc)
    A(i,:) = [1,xc(i),yc(i)];
end

X = A\eye(size(A));

%x arrays store alpha, beta, then gamma

Selem = zeros(length(xc));

for i=1:length(xc)
    for j=1:length(yc)
        Selem(i,j) = (elemArea)*((X(2,j) * X(2,i) )+(X(3,j) * X(3,i) ));
    end
end



function area = get_area(indices, x_total, y_total)
xc = x_total(indices);
yc = y_total(indices);
area  = 0.5*abs(det([ones(3,1), xc, yc]));
end


end

