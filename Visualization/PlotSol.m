function PlotSol(FilePath);

figure;

Set.filename = FilePath;
fd = PetscOpenFile(Set.filename);
Set.S = PetscBinaryRead(fd,'complex',true,'indices','int64','precision','float32');
%-vecload_block_size 1
close(fd);

n = round( length(Set.S)^(1/3) );
Set.S = reshape(Set.S, [n,n,n]);
x = zeros(n,1);
step = 1/n;
for j=2:n
    x(j) = x(j-1) + step;
end

Set.S = abs(Set.S);

slice(x,x,x,Set.S,[ 0 ],[],[]);
axis([0 1 0 1 0 1]);
title ('Slide cut of wave scattering field');