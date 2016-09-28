%% Propagator Method %%
% Forbidden region

function PsiArray=propForbid(E,s)
hbar2 = 0.076199682;
m=1;
V0=10;
beta = @(E) sqrt(2*m*(V0-E)/hbar2);

% forbidden propagator matrix
pF11 = @(x) cosh(beta(E)*s);
pF12 = @(x) (1/beta(E))*sinh(beta(E)*s);
pF21 = @(x) beta(E)*sinh(beta(E)*s);
pF22 = @(x) cosh(beta(E)*s);

PsiArray =  [pF11(E) pF12(E); pF21(E) pF22(E)];
end

