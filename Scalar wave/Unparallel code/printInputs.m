function printInputs(c,f,k,kappa,VolQ,a,d,M,P,alpha,N,n,n0,u0,p,h)
% Display the inputs arguments of wave scattering problem

str = '';
str = strcat(str, sprintf('\nSpeed of light in optics: %e',c));
str = strcat(str, sprintf('\nFrequency in optics: %e',f));
str = strcat(str, sprintf('\nWave number k = 2pi/lambda: %f',k));
str = strcat(str, sprintf('\nPower const with respect to the radius of particles: kappa is in [0,1]: %f',kappa));
str = strcat(str, sprintf('\nVolume of the domain Q that contains all particles: %f',VolQ));
str = strcat(str, sprintf('\nRadius of one particle: %e',a));
str = strcat(str, sprintf('\nDistance between two particles: d = O(a^(1/3)): %e',d));
str = strcat(str, sprintf('\nNumber of particles, M: %d',M));
str = strcat(str, sprintf('\nNumber of small cubes after partitioning the big Q, P: %d',P));
str = strcat(str, sprintf('\nDirection of plane wave, alpha: [%s]',num2str(alpha)));
if(nargin > 10)
    str = strcat(str, sprintf('\nContinuous distribution function of particles: %s',num2str(N)));
    str = strcat(str, sprintf('\nDesired refraction coefficient: %s',num2str(n)));
    str = strcat(str, sprintf('\nOriginal refraction coefficient: %s',num2str(n0)));
    str = strcat(str, sprintf('\nInitial field satisfies Helmholtz equation in R^3: %s',num2str(u0)));
    str = strcat(str, sprintf('\nFunction p(x): %s',num2str(p)));
    str = strcat(str, sprintf('\nFunction h(x): %s',num2str(h)));
end

disp(str);
    
end