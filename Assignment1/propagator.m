% Propagator Method %% 
% Allowed energy propagator for finite square well
% inputs: initial conditions = psi, Energy = E, well width = w,
%           separation of wells = w, 
%           propArray = [width of well, separation of wells, allowed, new]
%           allowed = 0/1, determines region
%           new = 0/1, determines whether to use psi input
% outputs: newPsi, and a check on the boundary conditions

function PsiArray=propagator(psi,E,propArray)
hbar2 = 0.076199682;
m=1;
V0=10;
propArray = num2cell(propArray);
[w,s,allowedSwitch,newSwitch] = propArray{:};
% defining alpha for allowable region, beta for forbidden region
alpha = sqrt(2.*m*E/hbar2);
beta = @(E) sqrt(2*m*(V0-E)/hbar2);

if newSwitch==1
    psi = [1 ; 1*beta(E)];
end

% allowable region propagator
if allowedSwitch==1
    % elements of the propagator matrix
    pA11 = @(x) cos(alpha.*(w));
    pA12 = @(x) (1./alpha)*sin(alpha.*(w));
    pA21 = @(x) -alpha.*sin(alpha.*(w));
    pA22 = @(x) cos(alpha.*(w));
    
    newPsi = [pA11(w) pA12(w); pA21(w) pA22(w)] * psi;
end

% forbidden zone propagator
if allowedSwitch==0
    pF11 = @(x) cosh(beta(E)*s);
    pF12 = @(x) (1/beta(E))*sinh(beta(E)*s);
    pF21 = @(x) beta(E)*sinh(beta(E)*s);
    pF22 = @(x) cosh(beta(E)*s);
    
    newPsi = [pF11(E) pF12(E); pF21(E) pF22(E)] * psi;
end

% check if boundary conditions are satisfied
bcCheck = @(psi) psi(1)*beta(E) + psi(2);

% return the new Psi as well as the boundary conditions check
PsiArray = [newPsi(1),newPsi(2),bcCheck(newPsi)];
end

