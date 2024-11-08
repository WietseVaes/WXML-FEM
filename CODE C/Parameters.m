%% Parameters function
% Choose a curve s (options in Mesh), amount spacing you want in your x,
% domain you want, define your f and define you g.

%                   2D: {'Kite'},

%                       {'Star',a,w}: a = indent (0<a<1); w = #petals
%                       {'Star',.3,7}

%                       {'Circles',c,R}: c = centers in array; R =
%                       radii in array;
%                       {'Circles',[0 + 0i],[2]}

%                       {'C Curve',a,b,c,d}, a = sharpennes (bigger =
%                       sharper); b = cavity controller (smaller = less
%                       cavity); c = fatness (larger = fatter all around);
%                       d = direction (1 = cavity left; exp(1i*(pi-theta)) = rotation theta)
%                       {'C Curve',3,2.5,1,1}
%                       {'C Curve',3,2.5,1,exp(-1i*3*pi/4)}

bnd_type = {'Star',.3,7};

% Define the domain range (xmin, xmax for x, ymin, ymax for y)
dom_range = {[-4, 4], [-4, 4]};  % Domain in x and y

% t interval of the parametrization for Dirichlet boundary condition full
% boudnary is [0,2*pi]
Dir_int = [0,2*pi];

% Set the number of points for discretization
if ~exist('n')
    n = 400;  % You can adjust this for more refined spacing
end

dx = (dom_range{1}(2)- dom_range{2}(1))/n;
T = 1;
Dt = 0.001;
%Dt = dx^2;
t = 0:Dt:T;
nt = length(t);


D = 1;

% Define the forcing term function f (interior source)
f = @(x, y) (x+y) .^0;  % Example function

% Define the boundary term function g ( Neumann boundary conditions)
g = @(x, y) (x+y) .^0 - 1;

% Define the source term function h (interior source)
h = @(x, y) (x+y) .^0 - 1;









%% Test case for accuracy testing

s = bd_curve(bnd_type);

usol = @(x,y,t) x .* cos(pi*y) .* sin(pi*t);

Dxu = @(x,y,t) cos(pi*y) .* sin(pi*t) +(x.^0 - 1);
Dxxu = @(x,y,t) (x+y+t).^0 - 1;

Dyu = @(x,y,t) -pi * x .* sin(pi*y) .* sin(pi*t);
Dyyu = @(x,y,t) - pi^2 * x .* cos(pi*y).* sin(t) .* sin(pi*t);

Dtu = @(x,y,t) pi * x .* cos(pi*y) .* cos(pi*t);

vx = @(x,y,t) x-y.^2+t;
Dxvx = @(x,y,t) (x+y+t).^0;

vy = @(x,y,t) cos(pi*x) .* y - t;
Dyvy = @(x,y,t) cos(pi*x) + (y+t).^0 - 1;


% forcing here
f = @(x,y,t) Dtu(x,y,t) + Dxvx(x,y,t) .* usol(x,y,t) + vx(x,y,t) .* Dxu(x,y,t) + Dyvy(x,y,t) .* usol(x,y,t) + vy(x,y,t) .* Dyu(x,y,t) - D*(Dxxu(x,y,t) + Dyyu(x,y,t));

g1 = @(x,y,t) s.nx(inv_gam(x,y,s)) .* (usol(x,y,t) .* vx(x,y,t) - D * Dxu(x,y,t));
g2 = @(x,y,t) s.ny(inv_gam(x,y,s)) .* (usol(x,y,t) .* vy(x,y,t) - D * Dyu(x,y,t));

g = @(x,y,t) -(g1(x,y,t) + g2(x,y,t));

h = @(x,y,t) usol(x,y,t);


function s = bd_curve(bnd_type)
    switch bnd_type{1}
        case 'Kite'
            a = 2;
            gam = chebfun(@(t) cos(t) + .65*cos(a*t)-.65 + 1i*(1.5*sin(t)), [0, 2*pi]);
            dgam = diff(gam);
            dgam2 = diff(dgam);
        case 'Star'
            a = bnd_type{2}; w = bnd_type{3};
            gam = chebfun(@(t) (2 + 2*a*cos(w*t)).*exp(1i*t), [0, 2*pi]);
        case 'Circles'
            c = bnd_type{2}; % centers
            R = bnd_type{3}; %radii
            gam = chebfun(@(t) R.*exp(1i*t)-c, [0,2*pi]);
        case 'C Curve'
            a = bnd_type{2}; %controls how "sharp" the corners are, bigger = sharper
            b = bnd_type{3}; %2.5 (smaller = less extreme cavity/less depth into the C)
            c = bnd_type{4}; % .1 (larger = "fatter all around")
            d = bnd_type{5}; % cavity opens to the left; choose as exp(1i*(0--2*pi)) to rotate
            t = chebfun(@(x) x, [0, 2*pi]);
            r = 3+c*tanh(a*cos(t));
            th = b*sin(t);
            gam = d*exp(1i*th).*r;
    end
    dgam = diff(gam);
    dgam2 = diff(dgam);
    s.Z = gam;
    s.Zp = dgam;
    
    s.nx = @(t) real(-1i .* s.Zp(t) ./ abs(s.Zp(t)));
    s.ny = @(t) imag(-1i .* s.Zp(t) ./ abs(s.Zp(t)));
end


function t = inv_gam(x,y,s)

t = zeros(length(x),1);
for i1 = 1:length(x)
    rts = roots(s.Z-(x(i1)+1i*y(i1)));
    if length(rts) == 2
        t(i1) = rts(1);
    else
        t(i1) = rts;
    end
end
end