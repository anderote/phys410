% guess finder

function guesses = guessFinder(results,E)

CA = sign(results);
guessN = 1;
for k =2:size(CA,2)-1
    if (CA(k) + CA(k+1))==0
        guesses{guessN}(1) = E(k);
        guesses{guessN}(2) = E(k+1);
        guessN = guessN + 1;
    end
end
end