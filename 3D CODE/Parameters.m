%% Parameters function
% Choose a curve s (options in Mesh), amount spacing you want in your x,
% domain you want, define your f and define you g.

% Define the domain range (xmin, xmax for x, ymin, ymax for y)
dom_range = {[-10, 10], [-10, 10], [-10, 10]};  % Domain in x and y
dom_range = {[0, 4], [0, 6], [0, 8]};
% Set the number of points for discretization
if ~exist('n')
    n = 10;  % You can adjust this for more refined spacing
end

T = 1;
Dt = 0.01;
t = 0:Dt:T;
nt = length(t);

D = 1;

%[x, y, z, elmat, elmatbd, In] = Mesh_rectangle(dom_range, n);

[x, y, z, Xorig, Yorig, Zorig, keep, elmat, elmatbd, In] = Mesh_room(dom_range,n);
% interior forcing 
f = (x+y+z+t).^0 - 1;  % Example function

% Dirichlet boundary condition
h = (x+y+z+t).^0 - 1;

% velocity field
vx = (x+y+z+t).^0 ;
vy = (x+y+z+t).^0 ;
vz = (x+y+z+t).^0 ;

% Flux, manual definition
g = zeros(length(x),nt); % Zeros everywhere, but...

Ibd = unique(elmatbd)';
%In = setdiff(Ibd, Id);
Id = setdiff(Ibd, In);

%Insulation everywhere that is not in a circle centered at x = .25 y = -.3 with radius
%.75 when z = top. (heated walls at Dirichlet bd). Incomming flux here over
%all time
for i1 = In
    %lets the flux be on the 2d circle not the 3d cylinder. 
    xc = x(i1); yc = y(i1); zc = z(i1);

    if (abs(yc - dom_range{2}(1)) <= eps)&& (1 <= xc) && (xc <= 2.5) && (zc < 4)
        g(i1,:) = 1; % -1 is incomming, since we take -g in pde.
    end

    if (yc == dom_range{2}(2)) && (1 <= xc) && (xc <= 3) && (zc < 5) && (zc >2)
        g(i1,:) = 1; % -1 is incomming, since we take -g in pde.
    end

    if (abs(yc - (2 * xc + 5)) <= eps) && (zc -5 < eps) && (zc >2)&& (0.1 - xc <= eps) && (xc-0.4 <= eps)
        g(i1,:) = 1; % -1 is incomming, since we take -g in pde.
    end

    if (abs(yc - (-2 * xc + 13)) <= eps) && (zc -5 < eps) && (zc-2  > eps) && (3.6 - xc <= eps) && (xc-3.9 <= eps)
        g(i1,:) = 1; % -1 is incomming, since we take -g in pde.
    end

end




scatter3(x(Ibd(g(Ibd,1) == 0)),y(Ibd(g(Ibd,1) == 0)),z(Ibd(g(Ibd,1) == 0)),".b"); hold on;
scatter3(x(g(:,1) ~= 0),y(g(:,1) ~= 0),z(g(:,1) ~= 0),"pentagram");
scatter3(x(Id),y(Id),z(Id),"*r"); 
xlabel("x"); ylabel("y"); zlabel("z"); legend("no Flux", "Incoming flux", "Dirichlet")



