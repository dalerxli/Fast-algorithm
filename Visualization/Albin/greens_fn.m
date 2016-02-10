function G = greens_fn( k, x, y )

r = sqrt( (x-y)'*(x-y) );
G = exp( i*k*r )/(4*pi*r);