%% Root finder %%
% Finds roots of the propagator method, corresponding to allowed energy states
% inputs: results = results of propagator on vector E
%         E = linspace of energies
%
% outputs: roots
% count the roots

function roots=root_finder(results,E,type)

beta = @(E) sqrt(2*m*(V0-E)/hbar2);
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
    fx0=matrix*[1 ; beta(x0)];
    
    x1 = guesses{k}(2);
    fx1=matrix*[1 ; beta(x1)];
    
    iter=1;
    change=abs(x1-x0);
    while ((abs(fx1)>tol ) || abs(fx0) > tol ) && (iter<maxIter) &&(change>tol)
        change=abs(x1-x0);
        iter= iter +1;
        c=secant(x1,x0,fx1,fx0);
        
        %update values for next iteration
        x0=x1;
        x1=c;
        fx0=matrix*[1 ; beta(x0)];
        fx1=matrix*[1 ; beta(x1)];
    end
    roots(k) = x1;
end
end
