function PlotField(FilePath,max);

figure;

Set.filename = FilePath;
fd = PetscOpenFile(Set.filename);
Set.S = PetscBinaryRead(fd,'complex',true,'indices','int64','precision','float32');
%-vecload_block_size 1
close(fd);

n = round( length(Set.S)^(1/3) );
Set.S = reshape(Set.S, [n,n,n]);

%max = 5;
step = 1/max;
[x,y,z] = meshgrid(0:step:step*(max-1),0:step:step*(max-1),0:step:step*(max-1));
step = floor(n/max);
u = z;
for j=0:max-1
    for k=0:max-1
        u(j+1,k+1) = Set.S(1,1+j*step,1+k*step); 
    end
end
clearvars  Set.S;
u = abs(u);

hold on;
quiver3(x,y,z,u,u,u);
hold off;
axis([0 1 0 1 0 1]);
title ('Wave scattering field');

