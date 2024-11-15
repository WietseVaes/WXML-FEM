function [s, xx, yy, keep, P] = setup_bd(bnd_type,n, rnge, shrinker)
% Boundary domain type: 
%                   2D: {'Kite'},

%                       {'Star',a,w}: a = indent (0<a<1); w = #petals
%                       {'Star',.3,7}

%                       {'Ci ircles',c,R}: c = centers in array; R =
%                       radii in array;
%                       {'Circles',[0 + 0i],[2]}

%                       {'C Curve',a,b,c,d}, a = sharpennes (bigger =
%                       sharper); b = cavity controller (smaller = less
%                       cavity); c = fatness (larger = fatter all around);
%                       d = direction (1 = cavity left; exp(1i*(pi-theta)) = rotation theta)  
%                       {'C Curve',3,2.5,1,1} 
%                       {'C Curve',3,2.5,1,exp(-1i*3*pi/4)}

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
      %       % Check for multiple centers
      %     if numel(c) > 1
      %         gam = chebfun(@(t) (R(1).* exp(1i * t) - c) + (R(2).* exp(1i * t) - c), [0, 2*pi]);
      %     else
      %         gam = chebfun(@(t) R(1).* exp(1i * t) - c(1), [0, 2*pi]);
      %     end
    case 'C Curve'
        a = bnd_type{2}; %controls how "sharp" the corners are, bigger = sharper
        b = bnd_type{3}; %2.5 (smaller = less extreme cavity/less depth into the C)
        c = bnd_type{4}; % .1 (larger = "fatter all around")
        d = bnd_type{5}; % cavity opens to the left; choose as exp(1i*(0--2*pi)) to rotate
        t = chebfun(@(x) x, [0, 2*pi]);
        r = 3+c*tanh(a*cos(t));
        th = b*sin(t);
        gam = d*exp(1i*th).*r;
    case 'Torus'
        R = bnd_type{2};
        r = bnd_type{3};%  radius
        %NN = length (R);
        %a = 2*pi / NN;
        %t = chebfun('t', [0,2*pi], 'splitting',' on');
        %gam = R(1)* exp(1i*t)-c(1);
        %for i1 = 2:length(c)
            %gam = join(gam,R(i1)*exp(i1*t) - c(i1));
        %end
        %Anulus definition (next 2 lines)
        t = chebfun('t',[0,pi],'splitting','on');
        gam = join(r*exp(2i*t),R*exp(2i*t));

end

dgam = diff(gam);
dgam2 = diff(dgam);
s.Z = gam;
s.original = s.Z;
s.Zp = dgam;
s.Zpp = dgam2;

s.xp = @(t) s.Zp(t);
s.sp = @(t) abs(s.xp(t) );
s.tang = @(t) s.xp(t)./s.sp(t);
s.nx = @(t) -1i*s.tang(t);  

%SHRINKER HAS ISSUE W TORUS


if ~strcmp(bnd_type{1}, 'Torus')
    s.Z = chebfun(@(t) gam(t) + shrinker*s.nx(t), [0, 2*pi]);
end

s.t = linspace(0,2*pi,100).';

[xx, yy, keep] = autosample_domain(n, rnge, s.Z,bnd_type);

xrnge = rnge{1};
h = (xrnge(2) - xrnge(1)) / n;
L = arcLength(gam);
N = floor(L/h);

%EQUIDISTRIBUTING POINTS ALONG A CONTOUR
if strcmp(bnd_type{1}, 'Torus')
    %I think you can just consider gam as one it doesnt seem to make a
    %differnce 
    % lenout = cumsum(abs(diff(R*exp(2i*t))));
    % lenin = cumsum(abs(diff(r*exp(2i*t))));
    % len = lenout + lenin;
    len = cumsum(abs(diff(gam))); 
    
    %This makes it worse! 
    g = inv(len);
    

    P = gam(g((0:N-1)*h));

else
    len = cumsum(abs(diff(gam)));
    g = inv(len);
    P = gam(g((0:N-1)*h));

end

 
end

function [xx, yy, keep] = autosample_domain(n, rnge, gam,bnd_type)
xrng = rnge{1}; yrng = rnge{2};
x = linspace(xrng(1), xrng(2), n).';
y = linspace(yrng(1), yrng(2), n);
[xx,yy] = meshgrid(x.', y);
zz = xx + 1i*yy;
sz = size(zz);
zz = zz(:);
%now keep only the points exterior to the object

keep = incontour(gam, {zz},bnd_type);
xx = xx(keep);
yy = yy(keep);


end

function out = incontour(gamma, xs,bnd_type)
% determine whether x is inside gamma (y = 1), or outside(y=0).
% here gamma is a closed contour.

x = xs{1};
t = linspace(0, 2*pi, 501);
t = t(1:end-1).';
%N = (length(t)+1)/2;
%sz = size(x);
% to avoid sampling too close to the boundary, we "push out" the 
% boundary in the normal direction. There will be issues here
% if the boundary has corners or cusps, or potentially related to
% nonconvexity. For now, it should work.


n = normal(gamma);
z = gamma(t);
zz = z;

%out = inpolygonc(x,zz);
if strcmp(bnd_type{1}, 'Torus')
    out = inpolygonc(x,gamma,bnd_type);
else
    out = inpolygonc(x,zz, bnd_type);
end
end

function T = inpolygonc(z,w,bnd_type) %complex version of inpolygon

    if strcmp(bnd_type{1}, 'Torus')

        outer_radius = bnd_type{2};
        inner_radius = bnd_type{3};
        %IMPLEMENT SHRINKER HERE
        in_outer_circle = ( (real(z)).^2 + (imag(z)).^2 ) < outer_radius^2;

        in_inner_circle = ( (real(z)).^2 + (imag(z)).^2 ) < inner_radius^2;

        % Points in the annulus are inside the outer circle but not the inner circle
        T = in_outer_circle & ~in_inner_circle;

        % Visualize the result

        x = real(z);
        y = imag(z);
        figure;
    else 
        [T, ON] = inpolygon(real(z),imag(z),real(w),imag(w));
        %correct so that points on the edge are not included:
        T(ON) = false;
    end



end