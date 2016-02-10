str = sprintf('WAVE SCATTERING PROBLEM: COMPARE SOLUTIONS OF TWO SYSTEMS S (ORIGINAL) AND P (REDUCED)\n');
disp(str);

% For optics:
% c = 3.0e+10; % Speed of light in optics
% f = 1.0e+14; % Frequency in optics
% k = 2*pi*f/c; % Wave number k = 2pi/lambda

% For acoustic waves:
c = 34400; % Speed of light in optics
f = 1000; % Frequency in optics
k = 2*pi*f/c; % Wave number k = 2pi/lambda

kappa = 0.9; % Power const with respect to the radius of particles: kappa in [0,1]
alpha = [1,0,0]; % alpha is a unit vector that indicates the direction of plane wave
VolQ = 1; % Volume of the domain Q that contains all particles
a = 10^(-3); % Radius of one particle
d = ((a^(2-kappa))/VolQ)^(1/3); % Distance between two particles: d = O(a^(1/3))
M = round(1/d)^3; % Number of particles: M = O(1/d^3)
P = round((M^(1/3))/7)^3;  % Number of small cubes after partitioning the big cube Q
N = 1; % Continuous distribution function of particles
n = sqrt(0.2); % Desired refraction coefficient
n0 = 1; % Original refraction coefficient
draw = 0;

str = sprintf('\nINPUT');
disp(str);
printInputs(c,f,k,kappa,VolQ,a,d,M,P,alpha);

tic
scattering3DS(a,kappa,alpha,d,M,N,n,n0,draw,k,VolQ);
E = ScatteringCompare3D_SP(a,kappa,alpha,d,M,P,N,n,n0,draw,k,VolQ);
str = sprintf('\nOUTPUT\n\nThe distance between the solutions of the S system (orignal) and the P system (reduced) is: %e', E);
disp(str);
toc

setpref('Internet','SMTP_Server','mail.math.ksu.edu');
setpref('Internet','E_mail','nhantran@math.ksu.edu');
sendmail('nhantran@math.ksu.edu', 'MATLAB 3D scattering', 'Scattering testing done!');