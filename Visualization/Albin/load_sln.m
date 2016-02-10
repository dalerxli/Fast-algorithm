%
% function P = load_sln( dirname, [suffix] )
%
% Loads the solution given in the directory 'dirname'.
%

function P = load_sln( dirname, varargin )

if( length(varargin) >= 1 )
    suffix = varargin{1};
else 
    suffix = '';
end

% find filenames
fns    = dir( [dirname, '/x', suffix, '_*.dat'] );

% initialize the patch structure
P = [];
c = 1;

% load each patch
for p = 1:length( fns )
    
    match = [ 'x', suffix, '_[0-9]*.dat'];
    if( length(regexp( fns(p).name, match )) > 0 );
        suf = fns(p).name(end-7:end-4);
        P(c).x  = load( [ dirname, '/x', suffix, '_',  suf, '.dat' ] );
        P(c).y  = load( [ dirname, '/y', suffix, '_',  suf, '.dat' ] );
        P(c).z  = load( [ dirname, '/z', suffix, '_',  suf, '.dat' ] );
        P(c).ur = load( [ dirname, '/ur', suffix, '_', suf, '.dat' ] );
        P(c).ui = load( [ dirname, '/ui', suffix, '_', suf, '.dat' ] );
        P(c).ua = abs( P(c).ur + i*P(c).ui );
        c = c + 1;
    end
    
end
