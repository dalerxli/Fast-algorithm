str = sprintf('WAVE SCATTERING PROBLEM: COMPARE SOLUTIONS OF TWO SYSTEMS S (ORIGINAL) AND P (REDUCED)\n');
disp(str);

c = 3*10^10; % Speed of light in optics
f = 10^14; % Frequency in optics
k = 2*pi*f/c; % Wave number k = 2pi/lambda
kappa = 0.9; % Power const with respect to the radius of particles: kappa in [0,1]
alpha = [1,0,0]; % alpha is a unit vector that indicates the direction of plane wave
VolQ = 1; % Volume of the domain Q that contains all particles
a = 10^(-4); % Radius of one particle
d = ((a^(2-kappa))/VolQ)^(1/3); % Distance between two particles: d = O(a^(1/3))
M = round(1/d)^3; % Number of particles: M = O(1/a)
P = round((M^(1/3))/7)^3;  % Number of small cubes after partitioning the big cube Q
N = ones(1,M); % Continuous distribution function of particles
n = ones(1,M).*sqrt(0.2); % Desired refraction coefficient
n0 = ones(1,M); % Original refraction coefficient
draw = 0;

str = sprintf('\nINPUT');
disp(str);
printInputs(c,f,k,kappa,VolQ,a,d,M,P,alpha);

E = ScatteringCompare3D_SP(a,kappa,alpha,d,M,P,N,n,n0,draw);
str = sprintf('\nOUTPUT\n\nThe distance between the solutions of the S system (orignal) and the P system (reduced) is: %e', E);
disp(str);