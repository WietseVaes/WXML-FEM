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

bnd_type = {'Kite'};

% Define the domain range (xmin, xmax for x, ymin, ymax for y)
dom_range = {[-2, 2], [-2, 2]};  % Domain in x and y

% t interval of the parametrization for Dirichlet boundary condition full
% boudnary is [0,2*pi]
Dir_int = [0,pi];

% Set the number of points for discretization
n = 50;  % You can adjust this for more refined spacing

D = 1;

% Define the forcing term function f (interior source)
f = @(x, y) (x+y) .^0;  % Example function

% Define the boundary term function g ( Neumann boundary conditions)
g = @(x, y) (x+y) .^0 - 1;

% Define the source term function h (interior source)
h = @(x, y) (x+y) .^0 - 1;









%% Test case for accuracy testing

s = bd_curve(bnd_type);

usol = @(x,y) x .* cos(pi*y);

Dxu = @(x,y) cos(pi*y) +(x.^0 - 1);
Dxxu = @(x,y) (x+y).^0 - 1;

Dyu = @(x,y) -pi * x .* sin(pi*y);
Dyyu = @(x,y) - pi^2 * x .* cos(pi*y);

vx = @(x,y) x-y.^2;
Dxvx = @(x,y) (x+y).^0;

vy = @(x,y) cos(pi*x) .* y;
Dyvy = @(x,y) cos(pi*x) + y.^0 - 1;


% forcing here
f = @(x,y) usol(x,y) + Dxvx(x,y) .* usol(x,y) + vx(x,y) .* Dxu(x,y) + Dyvy(x,y) .* usol(x,y) + vy(x,y) .* Dyu(x,y) - D*(Dxxu(x,y) + Dyyu(x,y));

g1 = @(x,y) s.nx(inv_gam(x,y,s)) .* (usol(x,y) .* vx(x,y) - D * Dxu(x,y));
g2 = @(x,y) s.ny(inv_gam(x,y,s)) .* (usol(x,y) .* vy(x,y) - D * Dyu(x,y));

g = @(x,y) -(g1(x,y) + g2(x,y));

h = @(x,y) usol(x,y);


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
    t(i1) = roots(s.Z-(x(i1)+1i*y(i1)));
end

end