function [X, Y, Z, Xorig, Yorig, Zorig, keep, elmat, elmatbd, In] = Mesh_room(dom_range, n)

%add keeep and all original 
x = linspace(dom_range{1}(1),dom_range{1}(2),n);
y = linspace(dom_range{2}(1),dom_range{2}(2),n);
z = linspace(dom_range{3}(1),dom_range{3}(2),n);
dx = (dom_range{1}(2) - dom_range{1}(1))/n;
dz = (dom_range{3}(2) - dom_range{3}(1))/n;

X = zeros(n^3,1);
Y = zeros(n^3,1);
Z = zeros(n^3,1);

j = 1;
for i1 = 1:length(x)
    for i2 = 1:length(y)
        for i3 = 1:length(z)
            X(j) = x(i1);
            Y(j) = y(i2);
            Z(j) = z(i3);
            j = j + 1;
        end
    end
end

%get rid of points outside LEW 230(corners by traingles) first 
Xorig = X;
Yorig = Y;
Zorig = Z;

%filter = ~(X < 0.5 & Y >5) & ~(X > 3.5 & Y > 5);
filter = ~(X < 0.5 & (2*X + 5) < Y) & ~(X > 3.5 & (-2*X + 13) < Y);

X = X(filter);
Y = Y(filter);
Z = Z(filter);

keep = filter;

nx = round(0.5/dx);
nz = round(8/dz);

x_leftwall = linspace(0, 0.5, n);
z_leftwall = linspace(0, 8, n);
x_rightwall = linspace(3.5, 4, n); 
z_rightwall = linspace(0, 8, n);

[XLWall, ZLWall] = meshgrid(x_leftwall, z_leftwall);
[XRWall, ZRWall] = meshgrid(x_rightwall, z_rightwall);
YLWall = 2*XLWall + 5;
XLWall = XLWall(:);
YLWall = YLWall(:);
ZLWall = ZLWall(:);

YRWall = -2* XRWall +13;
XRWall = XRWall(:);
YRWall = YRWall(:);
ZRWall = ZRWall(:);

size(X)
size(XLWall)
size(ZLWall)

X = [X;XLWall; XRWall];
Y = [Y;YLWall; YRWall];
Z = [Z;ZLWall; ZRWall];

size(X)

Xorig = [Xorig;XLWall; XRWall];
Yorig = [Yorig;YLWall; YRWall];
Zorig = [Zorig;ZLWall; ZRWall];

keep = [keep; ones(size(XLWall)); ones(size(XRWall))];

DT = delaunayTriangulation(X,Y,Z);

elmat = DT.ConnectivityList;
[elmatbd, ~] = convexHull(DT);
Ibd = unique(elmatbd(:))';

In = Ibd;








