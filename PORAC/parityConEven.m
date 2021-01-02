% returns equivalence conditions corresponding to even parity
% oblivious constraints from New Journal of Physics (18) 045003, 2016
function [oblivious,A] = parityConEven(AS,n,d)
%-------------------------------------------------------------------------------
% convert alice's setting to n,dits.
A = zeros(AS,n);
for a = 1:AS
    p=a-1;
    for b = 1:n
        A(a,n-b+1)=rem(p,d);
        p = (p-rem(p,d))/d;
    end
end
%-------------------------------------------------------------------------------
% parity calculations: parity contains n-1 parities for each Alice's input and parityd contains these parity represented in natural numbers.
parity = zeros(AS,n-1);
parityd = zeros(AS,1);
for a = 1:AS
    for p = 1:n-1
        parity(a,p) = mod(A(a,p)+A(a,p+1),d);
    end
    sum = 0;
    for o = 1:n-1
        sum = sum + parity(a,o)*(d^(o-1));
    end
    parityd(a) = sum;
end
%-------------------------------------------------------------------------------
% oblivious has alice's input that belong to same parity.
oblivious = zeros(d-1,d^(n-1));
count=ones(d^(n-1),1);
for a = 1:AS
    oblivious(count(parityd(a)+1),parityd(a)+1) = a;
    count(parityd(a)+1)= count(parityd(a)+1) + 1;
end
%-------------------------------------------------------------------------------