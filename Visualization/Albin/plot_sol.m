function plot_sol( dirname, fig );

figure(fig);

S = load( [ dirname, '/sol.dat' ] );

x = S(:,1); y = S(:,2); z = S(:,3);
v = S(:,4) + 1i*S(:,5);
n = round( length(x)^(1/3) );
x = x(1:n);
%[X,Y,Z] = meshgrid(x,x,x);
v = reshape( v, [n,n,n] );

slice( x, x, x, abs(v), [-0.4 0 0.4], [], [] );