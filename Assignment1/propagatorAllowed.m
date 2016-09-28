% Allowed energy propagator for finite square well
% inputs: initial conditions = psi, Energy = E, well width = b,
%           potential height = V0, particle mass =m;
%           switch for initial new==1, allowed ==1/0
% outputs:

function PsiArray=propagator(psi,E,b,V0,allowed,new)
hbar2 = 0.076199682;

% defining alpha for allowable region, beta for forbidden region
alpha = sqrt(2.*m*E/hbar2);
beta = @(E) sqrt(2*m*(V0-E)/hbar2);

if new==1
    psi = [1 ; 1*beta(E)];
end

% allowable region propagator
if allowed==1
    % elements of the propagator matrix
    p11 = @(x) cos(alpha.*(x-b));
    p12 = @(x) (1./alpha)*sin(alpha.*(x-b));
    p21 = @(x) -alpha.*sin(alpha.*(x-b));
    p22 = @(x) cos(alpha.*(x-b));
    
    newPsi = [p11(step) p12(step); p21(step) p22(step)] * psi;
end

% forbidden zone propagator
if allowed==0
    newPsi = [p11(step) p12(step); p21(step) p22(step)] * psi;
end

% check if boundary conditions are satisfied
bcCheck = @(psi) newPsi(1)*beta(E) + newPsi(2);

% return the new Psi as well as the boundary conditions check
PsiArray = [newPsi(1),newPsi(2),bcCheck(newPsi)];
end

