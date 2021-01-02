function S = Success(n,d,A,Rho,M);
AS = d^n;
S = 0;
for a = 1:AS
    for b = 1:(n)
        S = S + trace(Rho{a}*M{b,A(a,b)+1});
    end
end
S = real(S)/(AS*n);

