% nWell
% takes in n and e, returns the total propagator matrix for n allowed
% regions

function matrix=nWell(n,e)
w=0.6;
s=0.2;

iter=0;
j=1;
matrix=propAllowed(e,w);
for j=1:n
    if n>1
        if iter==1
            matrix=matrix*propAllowed(e,w);
        end
        
        if iter==0
            matrix=matrix*propForbid(e,s);
        end
        iter=~iter;
        j=j+1;
    end
end
end