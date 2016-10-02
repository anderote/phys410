% guess finder

function G = guessFinder(results,E)

n=rootCounter(results);
CA = sign(results);
guessN = 1;
G=cell(n,1);
for k =2:size(CA,2)-1
    if (CA(k) + CA(k+1))==0
        G{guessN}(1) = E(k);
        G{guessN}(2) = E(k+1);
        guessN = guessN + 1;
    end
end
G=G;
end