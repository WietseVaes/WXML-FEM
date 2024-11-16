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

bnd_type = {"Annulus", 3, 1};

% Define the domain range (xmin, xmax for x, ymin, ymax for y)
dom_range = {[-5, 5],[-5,5]};  % Domain in x and y

% t interval of the parametrization for Dirichlet boundary condition full
% boudnary is [0,2*pi]
Dir_int = [0.01,0.02];

% Set the number of points for discretization
if ~exist('n')
    n = 60;  % You can adjust this for more refined spacing
end

dx = (dom_range{1}(2)- dom_range{2}(1))/n;
T = 1;
Dt = 0.001;
%Dt = dx^2;
t = 0:Dt:T;
nt = length(t);




[x,y,elmat,elmatbd, Id, In] = MeshShrink(bnd_type, dom_range, n, Dir_int);


D = 10;

% Define the forcing term function f (interior source)
f =  (y + t) .^0 - 1;  % Example function

% Define the boundary term function g ( Neumann boundary conditions)
g = (x + t) .^0 - 1;

% Define the source term function h (interior source)
h = (x  + t) .^0 .^0 - 1;

vx = -5*y  + (t.^0-1);
vy = 5*x + (t.^0-1);

for i1 = unique(elmatbd).'
    if (sqrt(x(i1).^2 + y(i1).^2)) > bnd_type{3}+.1 && abs(angle(x(i1) + y(i1)*1i) - pi/2) < .1
        g(i1,:) = 15;
    end
    if (sqrt(x(i1).^2 + y(i1).^2)) < bnd_type{2}-.1 && abs(angle(x(i1) + y(i1)*1i) - pi/2) < .1*bnd_type{2}/bnd_type{3}
        g(i1,:) = -15;
    end
end




% % %% Test case for accuracy testing
% 
% s = bd_curve(bnd_type);
% 
% usol =  x .* cos(pi*y) .* sin(pi*t);
% 
% Dxu =  cos(pi*y) .* sin(pi*t) +(x.^0 - 1);
% Dxxu =  (x+y+t).^0 - 1;
% 
% Dyu = -pi * x .* sin(pi*y) .* sin(pi*t);
% Dyyu =  - pi^2 * x .* cos(pi*y).* sin(t) .* sin(pi*t);
% 
% Dtu = pi * x .* cos(pi*y) .* cos(pi*t);
% 
% vx = x-y.^2+t;
% Dxvx =  (x+y+t).^0;
% 
% vy =  cos(pi*x) .* y - t;
% Dyvy =  cos(pi*x) + (y+t).^0 - 1;
% 
% % 
% % s = bd_curve(bnd_type);
% % 
% % usol =  x .* cos(pi*y) .* sin(pi*t);
% % 
% % Dxu =  cos(pi*y) .* sin(pi*t) +(10);
% % Dxxu =  (x+y+t).^0 -1;
% % 
% % Dyu = -pi * x .* sin(pi*y) .* sin(pi*t);
% % Dyyu =  - pi^2 * x .* cos(pi*y).* sin(t) .* sin(pi*t);
% % 
% % Dtu = pi * x .* cos(pi*y) .* cos(pi*t);
% % 
% % vx = x-y.^2+t;
% % Dxvx =  (x+y+t).^0;
% % 
% % vy =  cos(pi*x) .* y - t;
% % Dyvy =  cos(pi*x) + (y+t).^0 - 1;
% % 
% 
% % forcing here
% f =  Dtu + Dxvx .* usol + vx .* Dxu + Dyvy .* usol + vy .* Dyu - D*(Dxxu + Dyyu);
% 
% 
% g1 = zeros(length(x),nt);
% g2 = zeros(length(x),nt);
% 
% for i1 = In
%     for i2 = 1:nt
%         g1(i1,i2) =  s.nx(inv_gam(x(i1),y(i1),s)) .* (usol(i1,i2) .* vx(i1,i2) - D * Dxu(i1,i2));
%         g2(i1,i2) =  s.ny(inv_gam(x(i1),y(i1),s)) .* (usol(i1,i2) .* vy(i1,i2) - D * Dyu(i1,i2));
%     end
% end
% 
% g =  -(g1 + g2);
% 
% h = usol;


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
        % case 'Circles'    
        %     R = bnd_type{3}; % Radii for circles
        %     c = bnd_type{2}; % Centers for circles
        % 
        %     % Check for multiple centers
        %     if numel(c) > 1
        %         gam = chebfun(@(t) (R(1).* exp(1i * t) - c) + (R(2).* exp(1i * t) - c), [0, 2*pi]);
        %     else
        %         gam = chebfun(@(t) R(1).* exp(1i * t) - c(1), [0, 2*pi]);
        %     end
        % 
        case 'C Curve'
            a = bnd_type{2}; %controls how "sharp" the corners are, bigger = sharper
            b = bnd_type{3}; %2.5 (smaller = less extreme cavity/less depth into the C)
            c = bnd_type{4}; % .1 (larger = "fatter all around")
            d = bnd_type{5}; % cavity opens to the left; choose as exp(1i*(0--2*pi)) to rotate
            t = chebfun(@(x) x, [0, 2*pi]);
            r = 3+c*tanh(a*cos(t));
            th = b*sin(t);
            gam = d*exp(1i*th).*r;
        case 'Annulus'
            R = bnd_type{2};
            r = bnd_type{3};%  radius
        %NN = length (R);
        %a = 2*pi / NN;
        %t = chebfun('t', [0,2*pi], 'splitting',' on');
        %gam = R(1)* exp(1i*t)-c(1);
        %for i1 = 2:length(c)
            %gam = join(gam,R(i1)*exp(i1*t) - c(i1));
        %end
            t = chebfun('t',[0,pi],'splitting','on');
            gam = join(r*exp(2i*t), R*exp(2i*t));

        %gam = chebfun(@(t) (R.* exp(1i * t) - c) + (r.* exp(1i * t) - c), [0, 2*pi]);  


    end
    dgam = diff(gam);
    dgam2 = diff(dgam);
    s.Z = gam;
    s.Zp = dgam;
    s.Zpp = dgam2;
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