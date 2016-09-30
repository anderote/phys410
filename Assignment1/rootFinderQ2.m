% rootfinding Q2

function roots = rootFinderQ2(n_roots,guesses,w,s)
tol=1e-12;
maxIter=1000;
secant = @(x1,x0,fx1,fx0) (x0*fx1 - x1*fx0)/(fx1-fx0);

for k=1:n_roots
    x0 = guesses{k}(1);
    fx0=bcCheck(propAllowed(x0,w)*propForbid(x0,s)*propAllowed(x0,w)*[1 ; beta(x0)],x0);
    
    x1 = guesses{k}(2);
    fx1=bcCheck(propAllowed(x1,w)*propForbid(x1,s)*propAllowed(x1,w)*[1 ; beta(x1)],x1);
    
    iter=1;
    change=abs(x1-x0);
    while ((abs(fx1)>tol ) || abs(fx0) > tol ) && (iter<maxIter) &&(change>tol)
        change=abs(x1-x0);
        iter= iter +1;
        c=secant(x1,x0,fx1,fx0);
        
        %update values for next iteration
        x0=x1;
        x1=c;
        fx0=bcCheck(propAllowed(x0,w)*propForbid(x0,s)*propAllowed(x0,w)*[1 ; beta(x0)],x0);
        fx1=bcCheck(propAllowed(x1,w)*propForbid(x1,s)*propAllowed(x1,w)*[1 ; beta(x1)],x1);
    end
    roots(k) = x1;
end

end