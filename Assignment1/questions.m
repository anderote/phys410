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
E=linspace(1e-9,10-1e-9,1000);

% apply the propagator method at each energy, and store the results of a
% boundary condition check in 'results'
index=1;
resultsQ1=zeros(size(E));

for e=E
    psi = [1 ; beta(e)];
    resultsQ1(index)=bcCheck(propAllowed(e,w)*psi,e);
    index = index + 1;
end

figure
plot(E,resultsQ1);
xlabel('Energy');
ylabel('f(E)');
title('Propagator method for root finding');

%% root finding for Q2

n_rootsQ1 = rootCounter(resultsQ1);
guessesQ1 = guessFinder(resultsQ1,E);

% perform secant method to find the roots to high tolerance
% x2 = x1 - f(x1)*(x1-x0)/(f(x1)-f(x2))
tol=1e-12;
maxIter=1000;
secant = @(x1,x0,fx1,fx0) (x0*fx1 - x1*fx0)/(fx1-fx0);
rootsQ1 = [];

for k=1:n_rootsQ1
    x0 = guessesQ1{k}(1);
    fx0=bcCheck(propAllowed(x0,w)*[1 ; beta(x0)],x0);
    
    x1 = guessesQ1{k}(2);
    fx1=bcCheck(propAllowed(x1,w)*[1 ; beta(x1)],x1);
    
    iter=1;
    change=abs(x1-x0);
    while ((abs(fx1)>tol ) || abs(fx0) > tol ) && (iter<maxIter) &&(change>tol)
        change=abs(x1-x0);
        iter= iter +1;
        c=secant(x1,x0,fx1,fx0);
        
        %update values for next iteration
        x0=x1;
        x1=c;
        fx0=bcCheck(propAllowed(x0,w)*[1 ; beta(x0)],x0);
        fx1=bcCheck(propAllowed(x1,w)*[1 ; beta(x1)],x1);
    end
    rootsQ1(k) = x1;
end

% check to see how close the roots are to zeroes; %% NOT WORKING SOME
% REASON
for k = 1:n_rootsQ1
    e=rootsQ1(k);
    error(k) = bcCheck(propAllowed(e,w)*[1;beta(e)],e);
end


%% Question 2: 
% double well of width 0.6nm, well separation = 0.2nm
w=0.6;
s=0.2;
m=1;
V0=10;
hbar2 = 0.076199682;
E=linspace(1e-9,10-1e-9,1000);
% we now have three boundary conditions, propagate 3 times: 

index = 1;
for e=E
    psi = [1;beta(e)];
    psiNew=propAllowed(e,w)*propForbid(e,s)*propAllowed(e,w)*psi;
    resultsQ2(index)=bcCheck(psiNew,e);
    index = index +1;
end

n_rootsQ2 = rootCounter(resultsQ2);
guessesQ2 = guessFinder(resultsQ2,E);

rootsQ2 = rootFinderQ2(n_rootsQ2,guessesQ2,w,s);