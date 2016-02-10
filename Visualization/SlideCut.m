function SlideCut(FilePath,N,SaveFile);

Set.filename = FilePath;
fd = PetscOpenFile(Set.filename);
Set.S = PetscBinaryRead(fd,'complex',true,'indices','int64','precision','float32');
%-vecload_block_size 1
close(fd);

n = round( length(Set.S)^(1/3) );
Set.S = reshape(Set.S, [n,n,n]);

x = zeros(N,1);
step = 1/N;
for j=2:N
    x(j) = x(j-1) + step;
end

S1 = zeros(N,N,N);
step = floor(n/N);
midpoint = round(step/2);
for j=1:N
    for k=1:N
        for l=1:N
            S1(j,k,l) = Set.S(midpoint+(j-1)*step,midpoint+(k-1)*step,midpoint+(l-1)*step);
        end
    end
end
clearvars  Set.S;

fid = fopen(SaveFile,'wt');
step = floor(N/5);
for j=1:step:N
    for k=1:step:N
        for l=1:step:N            
            fprintf(fid,'%f',real(S1(j,k,l)));
            if (imag(S1(j,k,l))>=0)
                fprintf(fid,'+');
            end
            fprintf(fid,'%fi , ',imag(S1(j,k,l)));
        end
        fprintf(fid,'\n');
    end
    fprintf(fid,'\n');
end
fclose(fid);
% dlmwrite(SaveFile,S1);

figure;
S1a = real(S1);
slice(x,x,x,S1a,[ 0.5 ],[],[]);
title ('Slice of wave scattering field: Real part');
axis([0 1 0 1 0 1]);
colorbar('horiz')

figure;
S1a = imag(S1);
slice(x,x,x,S1a,[ 0.5 ],[],[]);
title ('Slice of wave scattering field: Imaginary part'); 
axis([0 1 0 1 0 1]);
colorbar('horiz')

%shading interp;
%colormap(colorcube);
%h = colorbar('vert');
%set(h, 'xlim', [-1 2]);
%index = (-1:0.25:2);
%set(h,'YTickLabel',num2str(index'));



