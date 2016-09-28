%% Assignment 1 %%

function zero=propagator(E)
hbar2 = 0.076199682;
m=1;
b=0.6;
V0=10;

% defining alpha for allowable region, beta for forbidden region
alpha = sqrt(2.*m*E/hbar2);  
beta  = sqrt(2.*m*(V0-E)/hbar2);

% defining psi for region to the left, x<b psi(x)=C*e^(Bx), choose C=1
psi = 1;
psiPrime = 1*beta;

% elements of the propagator matrix
p11 = @(x,E,b) cos(alpha(E).*(x-b));
p12 = @(x,E,b) (1./alpha(E))*sin(alpha(E).*(x-b));
p21 = @(x,E,b) -alpha(E).*sin(alpha(E).*(x-b));
p22 = @(x,E,b) cos(alpha(E).*(x-b));

% define linspace for the allowable region
xlin = linspace(x,b,1000);








end