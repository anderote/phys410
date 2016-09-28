
function PsiArray=propAllowed(E,w)
hbar2 = 0.076199682;
m=1;
V0=10;
alpha = sqrt(2.*m*E/hbar2);

% elements of the propagator matrix
pA11 = @(x) cos(alpha.*(w));
pA12 = @(x) (1./alpha)*sin(alpha.*(w));
pA21 = @(x) -alpha.*sin(alpha.*(w));
pA22 = @(x) cos(alpha.*(w));

PsiArray = [pA11(w) pA12(w); pA21(w) pA22(w)];

end