function E = ScatteringCompare3D_SP(a,kappa,alpha,d,M,P,N,n,n0,draw,k,VolQ)
%DESCRIPTION: Solving scattering problem in a domain Q using LAS S (original system) and P (reduced system)
%SYNTAX     : ScatteringCompare3D_SP(M,u0,n0,n,d,a,N,P)
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
%             e.g. ScatteringCompare3D_SP -> no input values, it will take default values defined
%             below
%OUTPUT     : E   : The distance between 2 solutions of systems S and P
%AUTHOR     : NHAN TRAN - nhantran@math.ksu.edu

% INITIALIZING SOME CONSTS:

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
    % Number of small cubes after partitioning the big cube Q
    P = round((M^(1/3))/7)^3;
end
if (nargin < 7)
    % Continuous distribution function of particles
    N = ones(1,M);
end
if (nargin < 8)
    % Desired refraction coefficient
    n = ones(1,M).*sqrt(0.2);
end
if (nargin < 9)
    % Original refraction coefficient
    n0 = ones(1,M);
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
% Set position of particles based on N and d
[x,y,z] = positionParticles(d);

% SOLVING 2 SYSTEM S AND P:
str = sprintf('\nComputing...');
disp(str);
str = sprintf('\nSolving the reduced system P:');
disp(str);
tic
u2 = scattering3DP(a,kappa,alpha,d,M,P,N,n,n0,draw,k,VolQ);
toc
str = sprintf('\n2-norm of the solution of the reduced system P: %e',norm(u2,2));
disp(str);

str = sprintf('\nSolving the original system S:');
disp(str);
tic
u1 = scattering3DS(a,kappa,alpha,d,M,N,n,n0,draw,k,VolQ);
toc
str = sprintf('\n2-norm of the solution of the original system S: %e',norm(u1,2));
disp(str);

% FIND THE DISTANCE BETWEEN 2 SOLUTIONS u1 AND u2
orphanCount = 0;
dist = zeros(1,P);
nParticleCube = zeros(1,P);
for s1=1:M
    t = FindCube(s1);  % t <= P
    if(t > P)
        orphanCount = orphanCount + 1;
        continue;
    end
    diff = abs(u1(s1) - u2(t));
    dist(t) = dist(t) + diff;
    nParticleCube(t) = nParticleCube(t) + 1;
end
% Get average
for s1=1:P
    dist(s1) = dist(s1)/nParticleCube(s1);
end

if (orphanCount>0)
    str = sprintf('%d particles are not in any cube!', orphanCount);
    disp(str);
end

% The error/distance between 2 solutions u1 and u2
E = sum(dist);

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
        t1 = mod(s-1,b^2);
        x1 = floor(t1/b) + 1;
        % Find the position of [x1,x2,x3] in Cartesian coordinates
        xs = x(x1);
        ys = y(x2);
        zs = z(x3);
        
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function CN = FindCube(s)
        % Find cube number that contains particle s
        
        % % Number of particles on a side of a cube of size 1
        % b = ceil(M^(1/3));
        % % Number of small cubes on a side of a cube of size 1
        % nC = ceil(P^(1/3));
        % % Size of the big cube Q
        % sQ = (b-1)*d;
        % % Size of one small cube
        % sC = sQ/nC;
        
        [xs,ys,zs] = particle2position(s);
        xs = xs + 0.5;
        ys = ys + 0.5;
        zs = zs + 0.5;
        % Find the matrix order of the cube that contains particle s having position [xs,ys,zs]
        xc = floor(xs/sC) + 1;
        yc = floor(ys/sC) + 1;
        zc = floor(zs/sC) + 1;
        if (xc > nC)
            xc = nC;
        end
        if (yc > nC)
            yc = nC;
        end
        if (zc > nC)
            zc = nC;
        end
        
        % Return the cube number
        CN = yc + (xc - 1)*nC + (zc-1)*nC^2;
        
    end

end