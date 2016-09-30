% beta
function betaResult = beta(e)

hbar2 = 0.076199682;
m=1;
V0=10;
betaResult= sqrt(2*m*(V0-e)/hbar2);
end