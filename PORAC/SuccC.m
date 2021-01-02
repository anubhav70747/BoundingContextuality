function S = SuccC(n,d,A,Pe,Pd);
AS = d^n;
S = 0;
for a = 1:AS
    for m = 1:AS
        for b = 1:(n)
            S = S + Pe{a}(m)*Pd{m,b}(A(a,b)+1);
        end
    end
end
S = real(S)/(AS*n);
