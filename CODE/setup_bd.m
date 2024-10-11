function [s, xx, yy, keep] = setup_bd(bnd_type,n, rnge)
% Boundary domain type: 
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
        gam = @(t) chebfun(@(t) R.*exp(1i*t)-c);
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
s.Zpp = dgam2;
[xx, yy, keep] = autosample_domain(n, rnge, gam);

s.t = linspace(0,2*pi,1000);
s.xp = s.Zp(s.t);
s.sp = abs(s.xp);
s.tang = s.xp./s.sp;
s.nx = chebfun( -1i*s.tang, [0, 2*pi]);


end

function [xx, yy, keep] = autosample_domain(n, rnge, gam)
xrng = rnge{1}; yrng = rnge{2};
x = linspace(xrng(1), xrng(2), n).';
y = linspace(yrng(1), yrng(2), n);
[xx,yy] = meshgrid(x.', y);
zz = xx + 1i*yy;
sz = size(zz);
zz = zz(:);
%now keep only the points exterior to the object

keep = incontour(gam, {zz});
xx = xx(keep);
yy = yy(keep);
end

function out = incontour(gamma, xs)
% determine whether x is inside gamma (y = 1), or outside(y=0).
% here gamma is a closed contour.

if isscalar(xs)
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
    z = gamma(t); nn = n(t); nn = nn./sqrt(nn(:,1).^2 + nn(:,2).^2);
    nn = nn(:,1) + 1i*nn(:,2);
    zz = z + .05*nn;
    out = inpolygonc(x,zz);
else
    xx = xs{1}; yy = xs{2}; zz = xs{3};
    NN = 100;
    theta = linspace(0, 2*pi, NN);
    phi = linspace(0, 2*pi, NN);

    [Theta, Phi] = meshgrid(theta, phi);

    X = Theta*0; Y = Theta*0; Z = Theta*0;

    for i1 = 1:NN
        for i2 = 1:NN
            G = gamma(Theta(i1,i2),Phi(i1,i2));
            X(i1,i2) = G(1); Y(i1,i2) = G(2); Z(i1,i2) = G(3);
        end
    end

    fv = surf2patch(X, Y, Z, 'triangles');
    vertices = fv.vertices;
    faces = fv.faces;

    out = inpolyhedron(faces, vertices, [xx',yy',zz']);
end
%dgamma = diff(gamma);
%y = gamma(t);
%out = sum(pi/N*(1./(y-x.').*dgamma(t)));
%out = out(:);
%out(abs(out)>1e-5) = 1; %in contour
%out(abs(out)<=1e-5) = 0; %outsided contour
%out = logical(out);
end

function T = inpolygonc(z,w) %complex version of inpolygon
[T, ON] = inpolygon(real(z),imag(z),real(w),imag(w));
%correct so that points on the edge are not included:
T(ON) = false;
end