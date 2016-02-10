function S = Green3D(a,M)
%DESCRIPTION: Solving scalar wave scattering problem in 3D using 2 methods:
%             1. Standard matrix-vector multiplication
%             2. Convolution and FFT
%SYNTAX     : Green3D(a,d,M)
%INPUT      : a   : The radius of particles 
%             M   : Number of particles (Number of equations and unknows)
%OUTPUT     : SolDiff: The solution difference between the two methods
%AUTHOR     : NHAN TRAN - nhantran@math.ksu.edu

disp('SOLVING SCALAR WAVE SCATTERING PROBLEM BY MANY SMALL PARITCLES IN A UNIT CUBE:');
PI4 = 4*pi;
% For acoustic waves:
c = 34400; % Speed of wave
f = 1000; %Frequency
k = 2*pi*f/c; % Wave number k = 2pi/lambda
kappa = 0.99; % Power const with respect to the radius of particles: kappa in [0,1]
alpha = [1,0,0]; % alpha is a unit vector that indicates the direction of plane wave
%a = 10^(-6) % Radius of one particle
b = ceil(M^(1/3)); %Number of particles on a side of a unit cube
d = 1/(b-1); %Distance between neighboring particles
n = -1+0.001*1i; % Desired refraction coefficient n^2
n0 = 1; % Original refraction coefficient n_0^2
VolQ = 1; %Volume of the domain that contains all the particles
N = M*a^(2-kappa)/VolQ; % Continuous distribution function of particles

% p = p1 +i*p2
p = (k^2)*(n0*n0-n*n);
p1 = real(p);
p2 = imag(p);
% h = h1 +i*h2 is a continuous function, Im(h) <= 0
h1 = p1/(PI4*N);
h2 = p2/(PI4*N);
h = h1 + 1i*h2;

a2k = (a^(2-kappa))*PI4;
ha2k = h*a2k;

disp('Setup vector right-hand-side:');
u1 = InitField(b); %this is a cube
u2 = Convert3DTo1D(u1,b); %this is a vector

% Standard matrix-vector multiplication:
disp('Setup Green matrix:');
GM = GreenMatrix(M,d);
tic
disp('Standard matrix*vector:');
Sol1 = GM*u2;
toc

tic
% Matrix-vector multiplication using convolution and FFT:
sizes = 2*b-2;
disp('Setup Green cube:');
F_GC = Green3DCubePad(b,d);
disp('Matrix*vector using Convolution Theorem:');
F_GC = fftn(F_GC);
F_u1 = fftn(u1,[sizes sizes sizes]);
Sol2 = ifftn(F_GC.*F_u1);
Sol2 = Convert3DTo1D(Sol2,b);
Sol2 = Sol2(1:M);
toc

disp('Solving Ax=b using 2 different methods:');
tic
disp('1. Standard method:');
%S1 = GM\u2;
S1 = gmres(GM,u2,[],1e-3,10);
toc
tic
disp('2. Fast method:');
S2 = gmres(@Afunc,u2,10,1e-3,10);
toc

disp('Difference between 2 methods:');
MatVec = [Sol1 Sol2];
MatVecDiff = norm(Sol1-Sol2)
disp('Solutions of the scalar wave scattering problem solved by the 2 methods:');
SolDiff = norm(S1-S2)

% setpref('Internet','SMTP_Server','mail.math.ksu.edu');
% setpref('Internet','E_mail','nhantran@math.ksu.edu');
% sendmail('nhantran@math.ksu.edu', 'MATLAB scalar wave scattering', 'Solving Wave Scattering Done!');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function Ax = Afunc(x)
        x = Convert1DTo3D(x,size(x,1));
        x = fftn(x,[2*b-2 2*b-2 2*b-2]);
        Ax = ifftn(F_GC.*x);
        Ax = Convert3DTo1D(Ax,b);
        %Ax = Ax(1:M);
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function u0 = InitField(N)
        % Create an inittial field u0 satisfying Helmholtz equation in R^3        
        u0 = zeros(N,N,N,'single');
        
        for m3=1:N
            for m2=1:N
                for m1=1:N
                    %Particle's position in R^3:
                    x = (m1-1)*d;
                    y = (m2-1)*d;
                    z = (m3-1)*d;                    
                    u0(m1,m2,m3) = exp(1i*k*(alpha*[x,y,z]'));
                end
            end
        end         
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function Vec = Convert3DTo1D(Cube,N)        
        Vec = zeros(N^3,1,'single');
                
        for m3=1:N
            for m2=1:N
                for m1=1:N
                    m = Index2Order(m1,m2,m3);
                    Vec(m) = Cube(m1,m2,m3);
                end
            end
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function Cube = Convert1DTo3D(Vec,N) 
        n = round(N^(1/3));
        Cube = zeros(n,n,n,'single');
                
        for m=1:N
            [m1,m2,m3] = Order2Index(m);
            Cube(m1,m2,m3) = Vec(m);
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function m = Index2Order(m1,m2,m3)
        m = (m1-1)+(m2-1)*b+(m3-1)*b^2+1;
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function [m1,m2,m3] = Order2Index(m)
        m3 = floor((m-1)/b^2)+1;
        red1 = mod(m-1,b^2);
        m2 = floor(red1/b)+1;
        m1 = mod(red1,b)+1;
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function [x,y,z] = Particle2Position(m)
        % Set the position for each particle (uniformly distributed)
        % The first particle [x1,y1,z1] is at the left end bottom corner of the
        % cube and is called particle number 1.
        % The cube is in the first octant and the origin is one of its
        % vertex
        [m1,m2,m3] = Order2Index(m);
        x = (m1-1)*d;
        y = (m2-1)*d;
        z = (m3-1)*d;
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function GM = GreenMatrix(N,Step)
        GM = ones(N,N,'single');
        
        for s1=1:N-1
            for s2=s1+1:N
                [m1,m2,m3] = Order2Index(s1);
                [q1,q2,q3] = Order2Index(s2);
                GM(s1,s2) = Green3DF(abs(m1-q1),abs(m2-q2),abs(m3-q3),Step);
                GM(s2,s1) = GM(s1,s2);
            end
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function L = Length(x,y,z)
        % Compute the length of vector (x,y,z)
        L = sqrt(x^2 + y^2 + z^2);
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function G = Green3DF(m1,m2,m3,Step)
        % Create a Green function in 3D
%         if(m1==0 && m2==0 && m3==0)
%             G = 1;
%             return;
%         end
        
        % Distance 
        r = Step*Length(m1,m2,m3);
        G = ha2k*exp(1i*k*r)/(PI4*r);
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function GC = Green3DCube(N,Step)
        % Create a Green cube in 3D
        GC = zeros(N,N,N,'single');
        
        for m1=1:N
            for m2=m1:N
                for m3=1:N
                    GC(m1,m2,m3) = Green3DF(m1-1,m2-1,m3-1,Step);                    
                    GC(m2,m1,m3) = GC(m1,m2,m3);
                end
            end
        end 
        GC(1,1,1) = 1;        
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function GCP = Green3DCubePad(N,Step)
        % Create a Green cube in 3D
        GCP = Green3DCube(N,Step);
        %GCP = ones(sizes,sizes,sizes,'single');
        %GCP(1:N,1:N,1:N) = GC1;
        
        for m3=1:N
            for m2=1:N
                for m1=N+1:sizes
                    GCP(m1,m2,m3) = GCP(2*N-m1,m2,m3);
                end
            end
        end
        
        for m3=1:N
            for m2=N+1:sizes
                for m1=1:sizes
                    GCP(m1,m2,m3) = GCP(m1,2*N-m2,m3);
                end
            end
        end
        
        for m3=N+1:sizes
            for m2=1:sizes
                for m1=1:sizes
                    GCP(m1,m2,m3) = GCP(m1,m2,2*N-m3);
                end
            end
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%     function F = DFT(f)
%         % Discrete Fourier transform of vector f(x,y,z) in 3D. Same
%         % as fftn() of Matlab
%         [s1,s2,s3] = size(f);
%         F = zeros(s1,s2,s3);        
%         coef = -1i*2*pi/s1;
%         
%         for m1=1:s1
%             for m2=1:s2
%                 for m3=1:s3
%                     sum = 0;
%                     for p1=1:s1
%                         for p2=1:s2
%                             for p3=1:s3
%                                 sum = sum + exp(coef*(m1*p1+m2*p2+m3*p3))*f(p1,p2,p3);
%                             end
%                         end
%                     end                    
%                     F(m1,m2,m3) = sum;
%                 end
%             end
%         end
%     end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%     function f = IDFT(F)
%         % Inverse discrete Fourier transform of vector F(x,y,z) in 3D. Same
%         % as ifftn() of Matlab
%         [s1,s2,s3] = size(F);
%         f = zeros(s1,s2,s3);    
%         coef = 1i*2*pi/s1;
%         normal = 1/(s1^3);
%         
%         for m1=1:s1
%             for m2=1:s2
%                 for m3=1:s3
%                     sum = 0;
%                     for p1=1:s1
%                         for p2=1:s2
%                             for p3=1:s3
%                                 sum = sum + exp(coef*(m1*p1+m2*p2+m3*p3))*F(p1,p2,p3);
%                             end
%                         end
%                     end                    
%                     f(m1,m2,m3) = normal*sum;
%                 end
%             end
%         end
%     end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%     function FG = Conv3D(F,G)
%         % Convolution in 3D: Same as convn() of Matlab
%         [s1,s2,s3] = size(F);
%         FG = zeros(2*s1-1,2*s2-1,2*s3-1);    
%         
%         for m1=1:2*s1-1
%             for m2=1:2*s2-1
%                 for m3=1:2*s3-1
%                     sum = 0;
%                     for p1=max(1,m1+1-s1): min(m1,s1)
%                         for p2=max(1,m2+1-s2): min(m2,s2)
%                             for p3=max(1,m3+1-s3): min(m3,s3)
%                                 a1 = (m1-p1+1);
%                                 a2 = (m2-p2+1);
%                                 a3 = (m3-p3+1);
%                                 sum = sum + F(a1,a2,a3)*G(p1,p2,p3);
%                             end
%                         end
%                     end                    
%                     FG(m1,m2,m3) = sum;
%                 end
%             end
%         end
%     end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end