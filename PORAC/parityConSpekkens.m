% returns equivalence conditions corresponding to Spekkens' parity
% oblivious constraints from Phys. Rev. Lett. 102, 010401 (2009)
function [oblivious,A] = parityConSpekkens(AS,n,d)
%-------------------------------------------------------------------------------
% convert alice's setting to n,dits.
A = zeros(AS,n);
for a = 1:AS
   A(a,:)=dec2base(a-1,d,n) -'0';
end
% preparing parity strings
parityX = []; 
for x = 0:2^n -1
    y=dec2base(x,2,n) -'0';
    if sum(y)>=2
        parityX=[parityX;y];
    end
end
% calculating parities for each parity string
[j,k] = size(parityX);
parity = zeros(AS,j);
for a = 1:AS
    for p = 1:j
        parity(a,p) = mod(A(a,:)*parityX(p,:)',d);
    end
end
% collecting indices of states that belong to the same parity value for each parity
% string
oblivious = cell(j,1);
for p = 1:j
    oblivious{p} = zeros(d-1,d);
    count=ones(d,1);
    for a = 1:AS
        oblivious{p}(count(parity(a,p)+1),parity(a,p)+1) = a;
        count(parity(a,p)+1)= count(parity(a,p)+1) + 1;
    end
end