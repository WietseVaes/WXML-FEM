function [X, Y, Z, Xorig, Yorig, Zorig, keep, elmat, elmatbd, In] = Mesh_room(dom_range, n)

%add keeep and all original 
x = linspace(dom_range{1}(1),dom_range{1}(2),n);
y = linspace(dom_range{2}(1),dom_range{2}(2),n);
z = linspace(dom_range{3}(1),dom_range{3}(2),n);

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

DT = delaunayTriangulation(X,Y,Z);

elmat = DT.ConnectivityList;
[elmatbd, ~] = convexHull(DT);

%filter = ~(X < 0.5 & Y >5) & ~(X > 3.5 & Y > 5);
filter = ~(X < 1 & (X + 5) < Y) & ~(X > 3 & (-X + 9) < Y);

filtind = find(~filter);

% X = X(filter);
% Y = Y(filter);
% Z = Z(filter);

keep = filter;

for i1 = size(elmatbd,1):-1:1
    if any(ismember(elmatbd(i1,:), filtind))
        elmatbd(i1,:) = [];
    end
end

for i1 = size(elmat,1):-1:1
    if sum(~keep(elmat(i1,:)))>1
        elmat(i1,:) = [];
    elseif sum(~keep(elmat(i1,:)))==1
        elmatbd = [elmatbd; elmat(i1,keep(elmat(i1,:)))];
    end
end

% counter = 0;
% for i1 = flip(filtind .')
%     elmat(elmat >= i1) = elmat(elmat >= i1) - 1;
%     elmatbd(elmatbd >= i1) = elmatbd(elmatbd >= i1) - 1;
%     counter = counter + 1;
% end


 trisurf(elmatbd,X,Y,Z, 'FaceColor', 'cyan'); hold on

In = unique(elmatbd(:)).';








