function GF = Green1D(M)
PI4 = 4*pi;
k = 0;
%GF = zeros(M);
u = ones(M,1);

% disp('Setup Green matrix:');
% tic
% for s=1:M
%     for l=1:M
%         GF(s,l) = Green1D(s,l);        
%     end
% end
% toc
% 
% disp('Matrix Green*u:');
% tic
% GFu1 = GF*u;
% toc

disp('Setup Green function:');
tic
f = GreenFunc1D();
toc

disp('Matrix Green*u using convolution:');
tic
BigSize = max(size(f,1),size(u,1));
FGFu2 = fft(f,BigSize).*fft(u,BigSize);
GFu2 = ifft(FGFu2);

sGFu2 = size(GFu2,1);
GFu2 = GFu2(sGFu2-M+1:sGFu2);
toc

% disp('Difference between 2 methods:');
% norm(GFu1-GFu2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function G = Green1D(s,t)
        % Create a Green function in 3D
        
        if(s==t)
            G = 0;
            return;
        end
        
        % Distance from particle s to particle t in 3D
        r = abs(s-t);
        G = exp(1i*k*r)/(PI4*r);
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function GV = GreenFunc1D()
        % Create a Green function in 3D
        GV1 = zeros(M,1);
        GV2 = zeros(M-1,1);
        
        for s1=2:M
            r = abs(s1-1);
            GV1(s1) = exp(1i*k*r)/(PI4*r);
            GV2(M-s1+1) = GV1(s1);
        end
        
        GV = [GV2' GV1']';
    end

end