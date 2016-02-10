function S = scattering3DP(a,kappa,alpha,d,M,P,N,n,n0,draw,k,VolQ)
%DESCRIPTION: Solving scattering problem in a domain Q: Iu = u0 + Au, or u_j = u0_j - 4pi\sum_1^P {G(x_j,x_m)h(x_m)N(t)|small cube|u_m}
%SYNTAX     : scattering(M,u0,n0,n,a,d)
%INPUT      : M   : Number of equations and unknows (particles), or size of A
%             n0  : Original refraction coefficient
%             n   : Desired refraction coefficient
%             d   : Distance between two particles
%             a   : Radius of one particle a<<d
%             N   : Continuous distribution function of particles
%             P   : The number of small cubes when partitioning the domain Q of M particles
%             kappa: Power const with respect to the radius of particles: kappa in [0,1]
%             alpha: a unit vector that indicates the direction of plane wave
%             draw: Draw the cube Q if draw = 1
%             e.g. scattering3DP -> no input values, it will take default values defined
%             below
%OUTPUT     : S   : The solution to the scattering problem
%AUTHOR     : NHAN TRAN - nhantran@math.ksu.edu

% INITIALIZING SOME CONSTS:
% Speed of light in optics
% c = 3*10^10;
% % Frequency in optics
% f = 10^14;
% % Wave number k = 2pi/lambda
% k = 2*pi*f/c;
% % Volume of the domain Q that contains all particles
% VolQ = 1;
PI4 = pi*4;

% CHECKING INPUT VALUES:
if (nargin < 1)
    % Radius of one particle
    a = 10^(-3);
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
    % Number of small cubes after partitioning the big cube Q
    P = round((M^(1/3))/7)^3;
end
if (nargin < 7)
    % Continuous distribution function of particles
    N = 1;
end
if (nargin < 8)
    % Desired refraction coefficient
    n = sqrt(0.2);
end
if (nargin < 9)
    % Original refraction coefficient
    n0 = 1;
end
if (nargin < 10)
    % Draw the cube Q if draw = 1
    draw = 0;
end

% GLOBAL VARIABLES:
% Number of particles on a side of a cube of size 1
b = ceil(M^(1/3));
% Number of small cubes on a side of a cube of size 1
nC = ceil(P^(1/3));
% Size of the big cube Q
sQ = (b-1)*d;
% Size of one small cube
sC = sQ/nC;
%Position of small cubes in Cartesian coordinates
[x,y,z] = positionCubes();

% INITIALIZING SOME FUNCTIONS:
% Note:
% If p and h are discrete: all indices in p and h are integer (the order of particles)
% If p and h are continuous: all indices in p and h are real (the values of the position of particles)

% p = p1 +i*p2
p = (k^2)*(n0*n0-n*n);
p1 = real(p);
p2 = imag(p);
% h = h1 +i*h2 is a continuous function, Im(h) <= 0
h1 = p1/(PI4*N);
h2 = p2/(PI4*N);
h = h1 + 1i*h2;

% Initial field satisfies Helmholtz equation in R^3
u0P = initField();

% Matrix of the scattering system
A = formA();

% SOLVING THE SYSTEM OF EQUATION: Iu = u0 + Bu (A = I - B), or u_j = u0_j - 4pi\sum_1^P {G(x_j,x_m)h(x_m)N(t)|small cube|u_m}
S = gmres(A,u0P,[],1e-3,size(A,1));

%VISUALIZE
if (draw == 1)
    drawCube(P);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function [x,y,z] = positionCubes()
        % Set the position for each small cube (uniformly distributed)
        
        % % Number of particles on a side of a cube of size 1
        % b = ceil(M^(1/3));
        % % Number of small cubes on a side of a cube of size 1
        % nC = ceil(P^(1/3));
        % % Size of the big cube Q
        % sQ = (b-1)*d;
        % % Size of one small cube = distance between 2 cubes
        % sC = sQ/nC;
        
        x0 = -0.5-sC/2;
        y0 = -0.5-sC/2;
        z0 = -0.5-sC/2;
        
        % The first small cube [x1,y1,z1] is at the left end bottom corner of the
        % big cube Q and is called cube number 1.
        x = zeros(1,nC);
        y = zeros(1,nC);
        z = zeros(1,nC);
        for s=1:nC
            x(s) = x0 + sC*s;
            y(s) = y0 + sC*s;
            z(s) = z0 + sC*s;
        end
        
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function drawCube(P)
        
        % % Number of particles on a side of a cube of size 1
        % b = ceil(M^(1/3));
        % % Number of small cubes on a side of a cube of size 1
        % nC = ceil(P^(1/3));
        % % Size of the big cube Q
        % sQ = (b-1)*d;
        % % Size of one small cube = distance between 2 cubes
        % sC = sQ/nC;
        
        for s=1:P
            [xs,ys,zs] = cube2position(s);
            %plot3(xs,ys,zs,'s','MarkerSize',sC);
            drawACube(xs,ys,zs,sC);
            hold on;
        end
        box on;
        grid on;
        axis equal;
        AZ=-20;         % azimuth
        EL=25;          % elevation
        view(AZ,EL);    % orientation of the axes
        
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function drawACube(xc,yc,zc,size)
        
        x1=([0 1 1 0 0 0;1 1 0 0 1 1;1 1 0 0 1 1;0 1 1 0 0 0]-0.5)*size+xc;
        y1=([0 0 1 1 0 0;0 1 1 0 0 0;0 1 1 0 1 1;0 0 1 1 1 1]-0.5)*size+yc;
        z1=([0 0 0 0 0 1;0 0 0 0 0 1;1 1 1 1 0 1;1 1 1 1 0 1]-0.5)*size+zc;
        for i=1:6
            pa=patch(x1(:,i),y1(:,i),z1(:,i),'w');
            set(pa,'edgecolor','k')
        end
        
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function A = formA()
        % Create a matrix with Green-function values
        
        A = zeros(P);
        % % Number of particles on a side of a cube of size 1
        % b = ceil(M^(1/3));
        % % Number of small cubes on a side of a cube of size 1
        % nC = ceil(P^(1/3));
        % % Size of the big cube Q
        % sQ = (b-1)*d;
        % % Size of one small cube
        % sC = sQ/nC;
        % Volume of one small cube
        VolCube = sC^3;
        NhV = PI4*h*N*VolCube;        
        
        for s=1:P
            for t=1:P
                if (s~=t)
                    A(s,t) = green(s,t)*NhV;
                else
                    A(s,s) = 1;
                end
            end
        end
        
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function [xs,ys,zs] = cube2position(s)
        % Return the position in the 3D cube of the small cube s
        
        % The first small cube [x1,y1,z1] is at the left end bottom corner of the
        % big cube Q and is called cube number 1. The next one will be on the same
        % row, go to the right. When finishing the first line, go to the second line
        % and start at the first column again. When finishing the first plane, move
        % up.
        
        % % Number of small cubes on a side of a cube of size 1
        % nC = ceil(P^(1/3));
        
        % Find the plane where the particle s is on
        plane = floor((s-1)/(nC^2));
        % [x1,x2,x3] is an array index
        x3 = plane + 1;
        x2 = mod((s-1), nC) + 1;
        t = mod(s-1,nC^2);
        x1 = floor(t/nC) + 1;
        % Find the position of [x1,x2,x3] in Cartesian coordinates
        xs = x(x1);
        ys = y(x2);
        zs = z(x3);
        
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function G = green(s,t)
        % Create a Green function in 3D
        
        [xs,ys,zs] = cube2position(s);
        [xt,yt,zt] = cube2position(t);
        
        % Distance from cube s to cube t in 3D
        xy = sqrt((xs-xt)^2 + (ys-yt)^2 + (zs-zt)^2);
        G = exp(1i*k*xy)/(PI4*xy);
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function u0P = initField()
        % Create an inittial field u0 satisfying Helmholtz equation in R^3
        
        u0P = zeros(P,1);
        for s=1:P
            [xs,ys,zs] = cube2position(s);
            u0P(s) = exp(1i*k*(alpha*[xs,ys,zs]'));
        end
        
    end

end