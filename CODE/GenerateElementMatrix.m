% Generate the matrices of M, P and S over a considered element i1 through
% theory and numerical methods. you need to output Melem, Pelem, Selem
% which are three by three matrices in this case

function [Melem, Pelem, Selem] = GenerateElementMatrix(indices,x_total,y_total)

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
b1 = [1,0,0]';

for i=1:length(xc)
    A(i,:) = [1,xc(i),yc(i)];
end

x1 = A\b1; %solve Ax = b

%phi2

b2 = [0,1,0]';
x2 = A\b2;

%phi3

b3 = [0,0,1]';
x3 = A\b3;

X = [x1,x2,x3];
%x arrays store alpha, beta, then gamma

Selem = zeros(length(xc));

for i=1:length(xc)
    for j=1:length(yc)
        Selem(i,j) = (elemArea / 3)*((X(2,j) * X(2,i) )+(X(3,j) * X(3,i) ));
    end
end


%Pelem


function area = get_area(indices, x_total, y_total)
xc = x_total(indices);
yc = y_total(indices);
area  = 0.5*abs(det([ones(3,1), xc, yc]));
end

Pelem = 0;
end

