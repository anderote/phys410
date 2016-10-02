% Question 3 %%

% define constants of the simulation
tic();
run('constants.m');
N=1:30;

% Results are paramaterized by N
Espace = cell(size(N));
tResults = cell(size(N));
resultsQ3 = cell(size(N));
rootsQ3 = cell(size(N));
guessesQ3 = cell(size(N));
n_rootsQ3 = cell(size(N));

%root finding parameters and methods
tol=1e-10
secant = @(x1,x0,fx1,fx0) (x0*fx1 - x1*fx0)/(fx1-fx0);

parfor n=1:30;
    tic()
    Espace{n}=linspace(1e-9,10-1e-9,n*20);
    tResults{n} = zeros(size(Espace{n}));
    k=1;
    
    for e=Espace{n}
        psiNew = nWell(n,e)*[1;beta(e)];
        tResults{n}(k)=bcCheck(psiNew,e);
        k=k+1;
    end
    
    n_rootsQ3{n} = rootCounter(tResults{n});
    guessesQ3{n} = guessFinder(tResults{n},Espace{n});
    rootsQ3{n}=zeros(n_rootsQ3{n},1);
    
    % root finding
    for k=1:n_rootsQ3{n}
        x0 = guessesQ3{n}{k}(1);
        fx0=bcCheck(nWell(n,x0)*[1 ; beta(x0)],x0);
        
        x1 = guessesQ3{n}{k}(2);
        fx1=bcCheck(nWell(n,x1)*[1 ; beta(x1)],x1);
        
        change=abs(x1-x0);
        while ((abs(fx1)>tol ) || abs(fx0) > tol ) &&(change>tol)
            change=abs(x1-x0);
            c=secant(x1,x0,fx1,fx0);
            
            %update values for next iteration
            x0=x1;
            x1=c;
            fx0=bcCheck(nWell(n,x0)*[1 ; beta(x0)],x0);
            fx1=bcCheck(nWell(n,x1)*[1 ; beta(x1)],x1);
        end
        rootsQ3{n}(k) = x1;
    end
    toc()
end

% cleanup the data for plotting in a sensible manner
figure;
hold on
for k=N
    rootData = rootsQ3{k};
    data = zeros(length(rootData),2);
    data(:,1) = k;
    data(:,2) = rootData;
    scatter(data(:,1),data(:,2),'r.');
    title('Lowest energies as function of separation distance');
    xlabel('N wells');
    ylabel('Energy States (eV)');
    
end
toc()
