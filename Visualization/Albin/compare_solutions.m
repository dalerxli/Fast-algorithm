function compare_solutions_3( dir1, dir2 )

% solution file names
f1 = [ dir1, '/sol.dat' ];
f2 = [ dir2, '/sol.dat' ];

% load the solutions
S1 = load( f1 );
S2 = load( f2 );

% get the data
v1 = S1(:,4) + i*S1(:,5);
v2 = S2(:,4) + i*S2(:,5);

% convert to 3D arrays
n1 = round( size(v1,1)^(1/3) );
n2 = round( size(v2,1)^(1/3) );

v1 = reshape( v1, n1, n1, n1 );
v2 = reshape( v2, n2, n2, n2 );

% compute the divsor
d = n2/n1;

% extend v1
v3 = zeros( size( v2 ) );
for j =1:d
    for k =1:d
        for l = 1:d
            
            v3( d*(0:n1-1) + j, d*(0:n1-1)+k, d*(0:n1-1)+l ) ...
                = v1;
            
        end
    end
end

dV = (1/n2)^3;

% compute norms
inf_norm = max( max( max( abs( v2 ) ) ) );
two_norm = sqrt( sum( sum( sum( abs( v2 ).^2 ) ) ) * dV );

% compute errors
inf_err = max( max( max( abs( v2-v3 ) ) ) ) / inf_norm;
two_err = sqrt( sum( sum( sum( abs( v2-v3 ).^2 ) ) ) * dV ) / two_norm;

fprintf( '\n\n & %1.3e & %1.3e\n\n', inf_err, two_err );
