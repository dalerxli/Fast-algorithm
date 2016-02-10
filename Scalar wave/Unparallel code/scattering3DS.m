function S = scattering3DS(a,kappa,alpha,d,M,N,n,n0,u0,draw)
%DESCRIPTION: Solving scattering problem: Iu = u0 + Au, or u_j = u0_j - 4pi\sum_1^M {G(x_j,x_m)h(x_m)a^(2-kappa)u_m}
%SYNTAX     : scattering(M,u0,n0,n,a,d)
%INPUT      : M   : Number of equations and unknows (particles), or size of A
%             u0  : Initial field
%             n0  : Original refraction coefficient
%             n   : Desired refraction coefficient
%             d   : Distance between two particles
%             a   : Radius of one particle a<<d
%             N   : Continuous distribution function of particles
%             kappa: Power const with respect to the radius of particles: kappa in [0,1]
%             alpha: a unit vector that indicates the direction of plane wave
%             draw: Draw the cube Q if draw = 1
%             e.g. scattering3DS -> no input values, it will take default values defined
%             below
%OUTPUT     : S   : The solution to the scattering problem
%AUTHOR     : NHAN TRAN

% INITIALIZING SOME CONSTS:
% Speed of light in optics
c = 3*10^10;
% Frequency in optics
f = 10^14;
% Wave number k = 2pi/lambda
k = 2*pi*f/c;
% Volume of the domain Q that contains all particles
VolQ = 1;

% CHECKING INPUT VALUES:
if (nargin < 1)
    % Radius of one particle
    a = 10^(-2);
end
if (nargin < 2)
    % Power const with respect to the radius of particles: kappa in [0,1]
    kappa = 0.9;
end
if (nargin < 3)
    % alpha is a unit vector that indicates the direction of plane wave
    alpha = [1,0,0];
end
if (nargin < 4)
    % Distance between two particles: d = O(a^(1/3))
    d = ((a^(2-kappa))/VolQ)^(1/3);
end
if (nargin < 5)
    % Number of particles: M = O(1/a)
    M = round(1/d)^3;
end
if(M <= 0)
    error('Size of a matrix must be positive!');
end
if (nargin < 6)
    % Continuous distribution function of particles
    N = ones(1,M);
end

% GLOBAL VARIABLES:
% Number of particles on a side of a cube of size 1
b = ceil(M^(1/3));

%Position of particles
[x,y,z] = positionParticles(d);

if (nargin < 7)
    % Desired refraction coefficient
    n = ones(1,M).*sqrt(0.2);
end
if (nargin < 8)
    % Original refraction coefficient
    n0 = ones(1,M);
end
if (nargin < 9)
    % Initial field satisfies Helmholtz equation in R^3
    u0 = initField(alpha, k);
end
if (nargin < 10)
    % Draw the cube Q if draw = 1
    draw = 0;
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
h1 = zeros(1,M);
h2 = zeros(1,M);
for t1=1:M
    h1(t1) = p1(t1)/(4*pi*N(t1));
    h2(t1) = p2(t1)/(4*pi*N(t1));
end
h = h1 + 1i.*h2;

if (draw == 1)
    drawCube(M);
end

% Matrix of the scattering system
A = formA(M,k,h,a,kappa);
% Identity matrix
I = eye(M);

% SOLVING THE SYSTEM OF EQUATION: Iu = u0 + Au , or u_j = u0_j - 4pi\sum_1^M {G(x_j,x_m)h(x_m)a^(2-kappa)u_m}
S = gmres(I-A,u0');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function [x,y,z] = positionParticles(d)
        % Set the position for each particle (uniformly distributed)
        
        % % Number of particles on a side of a cube of size 1
        % b = ceil(M^(1/3));
        
        x0 = -0.5-d;
        y0 = -0.5-d;
        z0 = -0.5-d;
        
        % The first particle [x1,y1,z1] is at the left end bottom corner of the
        % cube and is called particle number 1.
        x = zeros(1,b);
        y = zeros(1,b);
        z = zeros(1,b);
        for s=1:b
            x(s) = x0 + d*s;
            y(s) = y0 + d*s;
            z(s) = z0 + d*s;
        end
        
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function drawCube(M)
        
        for s=1:M
            [xs,ys,zs] = particle2position(s);
            plot3(xs,ys,zs,'r.','MarkerSize',20);
            hold on;
        end
        box on;
        grid on;
        axis equal;
        
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function A = formA(M,k,h,a,kappa)
        % Create a matrix with Green-function values
        
        A = zeros(M);
        a2k = a^(2-kappa);
        for s=1:M
            for t=1:M
                if (s~=t)
                    A(s,t) = -4*pi*green(s,t,k)*h(t)*a2k;
                end
            end
        end
        
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function [xs,ys,zs] = particle2position(s)
        % Return the position in the 3D cube of particle s
        
        % The first particle [x1,y1,z1] is at the left end bottom corner of the
        % cube and is called particle number 1. The next one will be on the same
        % row, go to the right. When finishing the first line, go to the second line
        % and start at the first column again. When finishing the first plane, move
        % up.
        
        % % Number of particles on a side of a cube of size 1
        % b = ceil(M^(1/3));
        
        % Find the plane where the particle s is on
        plane = floor((s-1)/(b^2));
        % [x1,x2,x3] is an array index
        x3 = plane + 1;
        x2 = mod((s-1), b) + 1;
        t = mod(s-1,b^2);
        x1 = floor(t/b) + 1;
        % Find the position of [x1,x2,x3] in Cartesian coordinates
        xs = x(x1);
        ys = y(x2);
        zs = z(x3);
        
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function G = green(s,t,k)
        % Create a Green function in 3D
        
        [xs,ys,zs] = particle2position(s);
        [xt,yt,zt] = particle2position(t);
        
        % Distance from particle s to particle t in 3D
        xy = sqrt((xs-xt)^2 + (ys-yt)^2 + (zs-zt)^2);
        
        G = exp(1i*k*xy)/(4*pi*xy);
        
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function u0 = initField(alpha, k)
        % Create an inittial field u0 satisfying Helmholtz equation in R^3
        
        u0 = zeros(1,M);
        for s=1:M
            [xs,ys,zs] = particle2position(s);
            u0(s) = exp(1i*k*(alpha*[xs,ys,zs]'));
        end
        
    end


end