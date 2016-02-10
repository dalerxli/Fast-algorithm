function GF = Green(M)
PI4 = 4*pi;
k = 0;
b = ceil(M^(1/3));
%d = 1/b;
d = 1/(b-1);
GF = zeros(M);
%DistanceMatrix = zeros(M);
u = ones(M,1);

disp('Setup particles:');
tic
[x,y,z] = DistributeParticles();
Position = ParticlePosition();
toc

disp('Setup Green matrix:');
tic
for s=1:M
    for l=1:M
        GF(s,l) = Green3D(s,l);
        %DistanceMatrix(s,l) = Distance3D(s,l);
    end
end
toc
%disp('Distance among particles:');
%DistanceMatrix

disp('Matrix Green*u:');
tic
GFu1 = (GF*u);
GFu1'
toc

disp('Setup Green function:');
tic
f = GreenFunc3D();
toc

disp('Matrix Green*u using convolution:');
tic
BigSize = max(size(f,1),size(u,1));
FGFu2 = fft(f,BigSize).*fft(u,BigSize);
GFu2 = ifft(FGFu2);

sGFu2 = size(GFu2,1);
GFu2 = GFu2(sGFu2-M+1:sGFu2);
GFu2'
toc

disp('Difference between 2 methods:');
norm(GFu1-GFu2)

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

    function pos = ParticlePosition()
        % Return the position in the 3D cube of particle s
        
        pos = zeros(M,3);
        for s1=1:M
            [pos(s1,1),pos(s1,2),pos(s1,3)] = particle2position(s1);
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function [x,y,z] = DistributeParticles()
        % Set the position for each particle (uniformly distributed)
        
        % % Number of particles on a side of a cube of size 1
        % b = ceil(M^(1/3));
        t = -0.5;
        x0 = t-d;
        y0 = x0;
        z0 = x0;
        
        % The first particle [x1,y1,z1] is at the left end bottom corner of the
        % cube and is called particle number 1.
        x = zeros(1,b);
        y = zeros(1,b);
        z = zeros(1,b);
        
        for s1=1:b
            t = d*s1;
            x(s1) = x0 + t;
            y(s1) = y0 + t;
            z(s1) = z0 + t;
        end
        
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function r = Distance3D(s,t)
        % Distance from particle s to particle t in 3D
        r = sqrt((Position(s,1)-Position(t,1))^2 + (Position(s,2)-Position(t,2))^2 + (Position(s,3)-Position(t,3))^2);
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function DF = DistanceFunc3D()
        DF = zeros(1,1);
        sDF = 1;
        
        for j = 1:M
            for t = (j+1):M
                r = Distance3D(j,t);
                if(~IsInArray(DF,r))
                    DF = [DF r];
                    sDF = sDF + 1;
                end
            end
        end
        
        for j = sDF+1:M
            r = Distance3D(1,j);
            DF = [DF r];
        end        
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function In = IsInArray(A,r)
        In = 0;
        for j=1:size(A,2)
            if(A(j)==r)
                In = 1;
                break;
            end
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function G = Green3D(s,t)
        % Create a Green function in 3D
        if(s==t)
            G = 1;
            return;
        end
        
        % Distance from particle s to particle t in 3D
        r = Distance3D(s,t);
        G = exp(1i*k*r)/(PI4*r);
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function GV = GreenFunc3D()
        % Create a vector Green function in 3D
        DistanceFunc = DistanceFunc3D();
        sDF = size(DistanceFunc,2);
        GV1 = zeros(sDF-1,1);        
                
        for j=2:sDF               
            GV1(j-1) = exp(1i*k*DistanceFunc(j))/(PI4*DistanceFunc(j));           
        end        
             
        sGV1 = size(GV1,1);
        GV2 = GV1;
        for j=1:sGV1
            GV2(j) = GV1(sGV1-j+1);
        end                               
        
        GV = [GV2' 1 GV1']';        
    end

end