% Question 1:

optiCell = zeros(10,3);

for trial=1:20
    format long
    
    % define constants and energy space
    hbar2 = 0.076199682;
    w=0.6;
    s=0.6;
    m=1;
    V0=10;
    Espace = floor(10*log(trial+1));
    E=linspace(1e-9,10-1e-9,Espace);
    
    % apply the propagator method at each energy, and store the results of a
    % boundary condition check in 'results'
    index=1;
    resultsQ1=zeros(size(E));
    
    for e=E
        psi = [1 ; beta(e)];
        resultsQ1(index)=bcCheck(propAllowed(e,w)*psi,e);
        index = index + 1;
    end
    
%     figure
%     plot(E,resultsQ1);
%     xlabel('Energy');
%     ylabel('f(E)');
%     title('Propagator method for root finding');
%     
    % root finding for Q1
    
    n_rootsQ1 = rootCounter(resultsQ1);
    guessesQ1 = guessFinder(resultsQ1,E);
    
    % perform secant method to find the roots to high tolerance
    % x2 = x1 - f(x1)*(x1-x0)/(f(x1)-f(x2))
    tol=1e-12;
    maxIter=1000;
    secant = @(x1,x0,fx1,fx0) (x0*fx1 - x1*fx0)/(fx1-fx0);
    rootsQ1 = zeros(n_rootsQ1,1);
    secantIter = zeros(n_rootsQ1,1);
    
    for k=1:n_rootsQ1
        secantCount = 1;
        
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
            secantCount=secantCount+1;
        end
        rootsQ1(k) = x1;
        secantIterQ1(k)=secantCount;
    end
    
    % check to see how close the roots are to zeroes;
    for k = 1:n_rootsQ1
        e=rootsQ1(k);
        error(k) = bcCheck(propAllowed(e,w)*[1;beta(e)],e);
    end
    error=mean(error);
    
    optiCell(trial,1)=Espace;
    optiCell(trial,2)=mean(secantIterQ1);
    optiCell(trial,3)=error;
end

figure
hold on;
P1=subplot(2,1,1);
plot(optiCell(:,1),optiCell(:,2));
xlabel(P1,'Size of Espace');
ylabel(P1,'Secant iterations');
legend('Secant iterations')
title('Optimizing search space');

P2=subplot(2,1,2);
plot(optiCell(:,1),optiCell(:,3));
xlabel(P2,'Size of Espace');
ylabel(P2,'Error');
legend('Accuracy');
