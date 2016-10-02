% Question 2 %%
% double well of width 0.6nm, separation variable
% INTERESTING: Separations above 1.15 nm yield no energies allowable
tic();
run('constants.m');
sampleSize=300;
E=linspace(1e-9,10-1e-9,sampleSize);
% ensures we capture s=0.2 as part of the sample
S=linspace(1e-9,1,sampleSize);

% Results are paramterized by S
resultsQ2 = cell(size(S));
rootsQ2 = cell(size(S));
guessesQ2 = cell(size(S));
n_rootsQ2 = cell(size(S));
index=1;

for s=S;
    doubleWell = @(e) propAllowed(e,w)*propForbid(e,s)*propAllowed(e,w);
    resultsTemp = zeros(size(E));
    k=1;
    for e=E
        psi = [1;beta(e)];
        psiNew = doubleWell(e)*psi;
        resultsTemp(k)=bcCheck(psiNew,e);
        k=k+1;
    end
    
    n_rootsTemp = rootCounter(resultsTemp);
    if n_rootsTemp > 0
        guessesTemp = guessFinder(resultsTemp,E);
        rootsTemp = rootFinderQ2(n_rootsTemp,guessesTemp,w,s);
    
        % assign the results for this trial of S
        rootsQ2{index}=rootsTemp;
        n_rootsQ2{index}=n_rootsTemp;
        guessesQ2{index}=guessesTemp;
    end

    resultsQ2{index}=resultsTemp;
    index=index+1;
end

% find the lowest energies E1 and E2 from each separation distance
E1=[];
E2=[];

for j=1:size(S,2)
    E1(j)=rootsQ2{j}(1);
    E2(j)=rootsQ2{j}(2);
end

energiesQ2 = rootsQ2{21};

figure;
plot(S,E1,S,E2);
title('Lowest energies as function of separation distance');
xlabel('Separation distance (nm)');
ylabel('Energy (eV)');
legend('E1','E2');
toc()

