% Question 3 %%

% define constants of the simulation
tic();
w=0.6;
m=1;
V0=10;
s=0.2;
hbar2 = 0.076199682;
N=1:30;

% Results are paramaterized by N
resultsQ3 = cell(size(N));
rootsQ3 = cell(size(N));
guessesQ3 = cell(size(N));
n_rootsQ3 = cell(size(N));

index=1;
for n=N;
    E=linspace(1e-9,10-1e-9,1000);
    resultsTemp = zeros(size(E));
    k=1;
    
    for e=E
        psi = [1;beta(e)];
        psiNew = nWell(n,e)*psi;
        resultsTemp(k)=bcCheck(psiNew,e);
        k=k+1;
    end
    
    n_rootsTemp = rootCounter(resultsTemp);
    if n_rootsTemp>0
        guessesTemp = guessFinder(resultsTemp,E);
        rootsTemp = rootFinderQ3(n_rootsTemp,guessesTemp,n);
        rootsQ3{index}=rootsTemp;
        n_rootsQ3{index}=n_rootsTemp;
        guessesQ3{index}=guessesTemp;
    end
    
    % assign the results for this trial of S
    resultsQ3{index}=resultsTemp;
    index=index+1;
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
