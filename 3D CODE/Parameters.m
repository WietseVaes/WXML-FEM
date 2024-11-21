%% Parameters function
% Choose a curve s (options in Mesh), amount spacing you want in your x,
% domain you want, define your f and define you g.

% Define the domain range (xmin, xmax for x, ymin, ymax for y)
dom_range = {[0, 4], [0, 6], [0, 5]};

% Set the number of points for discretization
if ~exist('n')
    n = 40;  % You can adjust this for more refined spacing
end

T = 1;
Dt = 0.01;
t = 0:Dt:T;
nt = length(t);

D = 10;

%[x, y, z, elmat, elmatbd, In] = Mesh_rectangle(dom_range, n);

[x, y, z, Xorig, Yorig, Zorig, keep, elmat, elmatbd, In] = Mesh_room(dom_range,n);
% interior forcing 
f = (x+y+z+t).^0 - 1;  % Example function

% Dirichlet boundary condition
h = (x+y+z+t).^0 - 1;

% velocity field
vx = (x+y+z+t).^0 - 1 ;
vy = -10*(x+y+z+t).^0 ;
vz = -(x+y+z+t).^0 ;

% Flux, manual definition
g = zeros(length(x),nt); % Zeros everywhere, but...

Id = [];

%Insulation everywhere that is not in a circle centered at x = .25 y = -.3 with radius
%.75 when z = top. (heated walls at Dirichlet bd). Incomming flux here over
%all time

for i1 = In
    xc = x(i1); yc = y(i1); zc = z(i1);

    if (abs(yc - dom_range{2}(1)) <= eps)&& (xc <= 2*(dom_range{1}(2)-dom_range{1}(1))/5) && (xc >= (dom_range{1}(2)-dom_range{1}(1))/5) && (zc-3*dom_range{3}(2)/5 <= eps) %DOOR
        g(i1,:) = -1; % -1 is incomming, since we take -g in pde.
    elseif (zc <= 2*(dom_range{3}(2)-dom_range{3}(1))/3) && (zc >= 1*(dom_range{3}(2)-dom_range{3}(1))/3)
        if (abs(yc - dom_range{2}(2)) <= eps) && (xc <= 2*(dom_range{1}(2)-dom_range{1}(1))/3) && (xc >= (dom_range{1}(2)-dom_range{1}(1))/3)           
            g(i1,:) = 1;
        elseif (yc - 5.25 > eps) && (yc - 5.75 < eps) 
            if (xc - 0.25 > eps) && (xc - 0.75 < eps) 
                g(i1,:) = 1; % -1 is incomming, since we take -g in pde.
            elseif (xc-3.75 <= eps) && (3.25 - xc <= eps)
                g(i1,:) = 1;
            end
        end
    end

end


scatter3(x(In(g(In,1) == 0)),y(In(g(In,1) == 0)),z(In(g(In,1) == 0)),".k"); hold on;
scatter3(x(In(g(In,1) ~= 0)),y(In(g(In,1) ~= 0)),z(In(g(In,1) ~= 0)),"pentagram");
scatter3(x(Id),y(Id),z(Id),"*r"); 
xlabel("x"); ylabel("y"); zlabel("z"); legend("no Flux", "Incoming flux", "Dirichlet")



