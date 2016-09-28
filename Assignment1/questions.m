%% Assignment 1 %%


%% Question 1:
% n_zeros = sum(abs(diff(sign(zeros)))/2);
% define energy linspace to investigate a region of energies
% define propArray to be constants used for propagation method
format long
% define constants and energy space
hbar2 = 0.076199682;
w=0.6;
s=0.6;
m=1;
V0=10;
beta = @(E) sqrt(2*m*(V0-E)/hbar2);
E=linspace(1e-9,10-1e-9,1000);
bcCheck = @(psi,e) psi(1)*beta(e) + psi(2);
% apply the propagator method at each energy, and store the results of a
% boundary condition check in 'results'
index=1;
results=zeros(size(E));

for e=E
    psi = [1 ; beta(e)];
    results(index)=bcCheck(propAllowed(e,w)*psi,e);
    index = index + 1;
end

figure
plot(E,results);
xlabel('Energy');
ylabel('f(E)');
title('Propagator method for root finding');

%% root finding for Q2
n_root = sum(floor(abs(diff(sign(results)))/2),'omitnan');
% find starting guesses by looking for sign changes in f(E), then
% extracting corresponding elements of E
guesses = cell(n_root,1);
CA = sign(results);
guessN = 1;
for k =2:size(CA,2)-1
    if (CA(k) + CA(k+1))==0
        guesses{guessN}(1) = E(k);
        guesses{guessN}(2) = E(k+1);
        guessN = guessN + 1;
    end
end

% perform secant method to find the roots to high tolerance
% x2 = x1 - f(x1)*(x1-x0)/(f(x1)-f(x2))
tol=1e-9;
maxIter=1000;
secant = @(x1,x0,fx1,fx0) (x0*fx1 - x1*fx0)/(fx1-fx0);
roots = [];

for k=1:n_root
    x0 = guesses{k}(1);
    fx0=propAllowed(x1,w)*[1 ; beta(x0)];
    
    x1 = guesses{k}(2);
    fx1=propAllowed(x1,w)*[1 ; beta(x1)];
    
    iter=1;
    change=abs(x1-x0);
    while ((abs(fx1)>tol ) || abs(fx0) > tol ) && (iter<maxIter) &&(change>tol)
        change=abs(x1-x0);
        iter= iter +1;
        c=secant(x1,x0,fx1,fx0);
        
        %update values for next iteration
        x0=x1;
        x1=c;
        fx0=propAllowed(x0,w)*[1 ; beta(x0)];
        fx1=propAllowed(x1,w)*[1 ; beta(x1)];
    end
    roots(k) = x1;
end

% check to see how close the roots are to zeroes;
for k=1:4
    resultTemp = propAllowed(roots(k),w);
    error(k)=resultTemp(3);
end


%% Question 2: 
% double well of width 0.6nm, well separation = 0.2nm
w=0.6;
s=0.2;
m=1;
V0=10;
hbar2 = 0.076199682;

% we now have three boundary conditions, propagate 3 times: 

bcCheck = @(psi) psi(1)*beta(E) + psi(2);

for e=E
    psi = [1;beta(e)];
    psiNew=propAllowed(e,w)*propForbid(e,s)*propAllowed(e,w)*psi;
    resultsQ2=bcCheck(psiNew);
end

n_root = sum(floor(abs(diff(sign(resultsQ2)))/2),'omitnan');
