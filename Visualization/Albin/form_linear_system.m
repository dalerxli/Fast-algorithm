function [A,b] = form_linear_system( Nh, P )

k  = 1;

dx = 1/P;
x = -0.5 + ((1:P)-0.5)*dx;
dV = dx^3;

A = eye( P^3 );
b = zeros( P, 1 );

% loop over targets
for loc_p = 1:P^3
    
    % compute local coordinates
    i_p = [ 1 + mod( loc_p-1, P ), 1 + mod( floor((loc_p-1)/P), P ), ...
        1 + floor((loc_p-1)/P^2) ];
    
    x_p = x(i_p)';
    
    b(loc_p) = exp( -i*k*x_p(3) );
   
    % loop over sources
    for loc_q = 1:P^3
        
        % skip if source is target
        if( loc_p == loc_q )
            continue
        end

        % compute local coordinates
        i_q = [ 1 + mod( loc_q-1, P ), 1 + mod( floor((loc_q-1)/P), P ), ...
            1 + floor((loc_q-1)/P^2) ];
        
         x_q = x(i_q)';
         
         A(loc_p, loc_q) = 4*pi*greens_fn(k,x_p,x_q)*Nh*dV;
            
        
    end
    
end

