function S = scattering(M,u0,n0,n,d,a,N)
%DESCRIPTION: Solving scattering problem: Iu = u0 + Au, or u_j = u0_j - 4pi\sum_1^M {G(x_j,x_m)h(x_m)a^(2-kappa)u_m}
%SYNTAX     : scattering(M,u0,n0,n,a,d)
%INPUT      : M  : Number of equations and unknows (particles), or size of A
%             u0 : Initial field
%             n0 : Original refraction coefficient
%             n  : Desired refraction coefficient
%             d  : Distance between two particles
%             a  : Radius of one particle
%             N  : Continuous distribution function of particles
%             e.g. scattering -> no input values, it will take default values defined
%             below
%OUTPUT     : S  : The solution to the scattering problem
%AUTHOR     : NHAN TRAN

% INITIALIZING SOME CONSTS:
% Speed of light in optics
c = 3*10^10;
% Frequency in optics
f = 10^14;
% Wave number k = 2pi/lambda
k = 2*pi*f/c;
% Power const with respect to the radius of particles: kappa in [0,1]
kappa = 0.9;
% Volume of the domain Q that contains all particles
VolQ = 1;

% CHECKING INPUT VALUES:
if (nargin < 6)
    % Radius of one particle
    a = 10^(-2);
end
if (nargin < 1)
    % Number of particles: M = O(1/a)
    M = round(1/a);
end
if(M <= 0)
    error('Size of a matrix must be positive!');
end
if (nargin < 7)
    % Continuous distribution function of particles
    N = ones(1,M);
end
if (nargin < 5)
    % Distance between two particles: d = O(a^(1/3))
    d = ((a^(2-kappa))/VolQ)^(1/3);
end
if (nargin < 4)
    % Desired refraction coefficient
    n = ones(1,M).*sqrt(0.2);
end
if (nargin < 3)
    % Original refraction coefficient
    n0 = ones(1,M);
end
if (nargin < 2)
    % Initial field satisfies Helmholtz equation in R^3
    u0 = ones(1,M);
end

% INITIALIZING SOME FUNCTIONS: 
% Note: 
% If p and h are discrete: all indices in p and h are integer (the order of particles)
% If p and h are continuous: all indices in p and h are real (the values of the position of particles)

% p = p1 +i*p2
p = (k^2).*(n0.*n0-n.*n);
p1 = real(p);
p2 = imag(p);
% h = h1 +i*h2 is a continuous function, Im(h) <= 0
for t=1:M
    h1(t) = p1(t)/(4*pi*N(t));
    h2(t) = p2(t)/(4*pi*N(t));
end
h = h1 + i.*h2;

%printInputs(c,f,k,kappa,VolQ,N,a,d,n,n0,u0,M,p,h);

%Position of particles
[x,y,z] = position(M,d);

% Matrix of the scattering system
A = formA(M,x,y,z,k,h,a,kappa);
% Identity matrix
I = eye(M); 

% SOLVING THE SYSTEM OF EQUATION: Iu = u0 + Au , or u_j = u0_j - 4pi\sum_1^M {G(x_j,x_m)h(x_m)a^(2-kappa)u_m}
S = gmres(I-A,u0');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function printInputs(c,f,k,kappa,VolQ,N,a,d,n,n0,u0,M,p,h)

display(' Speed of light in optics:');
display(c);
display(' Frequency in optics:');
display(f);
display(' Wave number k = 2pi/lambda:');
display(k);
display(' Power const with respect to the radius of particles: kappa is in [0,1]:');
display(kappa);
display(' Volume of the domain Q that contains all particles:');
display(VolQ);
display(' Continuous distribution function of particles:');
display(N);
display(' Radius of one particle:');
display(a);
display(' Distance between two particles: d = O(a^(1/3)):');
display(d);
display(' Desired refraction coefficient:');
display(n);
display(' Original refraction coefficient:');
display(n0);
display(' Initial field satisfies Helmholtz equation in R^3:')
display(u0);
display(' Number of particles, M = O(1/a):');
display(M);
display(' Function p(x):');
display(p);
display(' Function h(x):');
display(h);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x,y,z] = position(M,d)
% Set the position for each particle (uniformly distributed)

% Number of particles on a side of a cube of size 1
b = ceil(M^(1/3));

x0 = -0.5-d;
y0 = -0.5-d;
z0 = -0.5-d;

% The first particle [x1,y1,z1] is at the left end bottom corner of the
% cube and is called particle number 1.

for s=1:b
    x(s) = x0 + d*s;
    y(s) = y0 + d*s;
    z(s) = z0 + d*s;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function A = formA(M,x,y,z,k,h,a,kappa)
% Create a matrix with Green-function values

A = zeros(M);
a2k = a^(2-kappa);
for s=1:M
    for t=1:M
        if (s~=t)
            A(s,t) = -4*pi*green(s,t,x,y,z,k,M)*h(t)*a2k;
        end
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xs,ys,zs] = particle2position(s,x,y,z,M)
% Return the position in the 3D cube of particle s

% The first particle [x1,y1,z1] is at the left end bottom corner of the
% cube and is called particle number 1. The next one will be on the same
% row, go to the right. When finishing the first line, go to the second line
% and start at the first column again. When finishing the first plane, move
% up.

% Number of particles on a side of a cube of size 1
b = ceil(M^(1/3));
% Find the plane where the particle s is on
plane = floor((s-1)/(b^2));
% [x1,x2,x3] is an array index
x3 = plane + 1;
x2 = mod((s-1), b) + 1;
t = mod(s-1,b^2) + 1;
x1 = floor((t-1)/b) + 1;
% Find the position of [x1,x2,x3] in Cartesian coordinates
xs = x(x1);
ys = y(x2);
zs = z(x3);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function G = green(s,t,x,y,z,k,M)
% Create a Green function in 3D

[xs,ys,zs] = particle2position(s,x,y,z,M);
[xt,yt,zt] = particle2position(t,x,y,z,M);

% Distance from particle s to particle t in 3D
xy = sqrt((xs-xt)^2 + (ys-yt)^2 + (zs-zt)^2);

G = exp(i*k*xy)/(4*pi*xy);

end

end