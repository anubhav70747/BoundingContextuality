function [S,SB,SE] = SuccessBE(n,d,q,A,Rho,MB,ME,alp);
AS = d^n;
SB = 0;
for a = 1:AS
    for b = 1:(n)
        SB = SB + trace(Rho{a}*kron(MB{b,A(a,b)+1},eye(q)));
end
end
SE=0;
for a = 1:AS
    for e = 1:(n)
        SE = SE + trace(Rho{a}*kron(eye(q),ME{e,A(a,e)+1}));
end
end
SBE = 0;

for a = 1:AS
    for b = 1:n
    for e = 1:(n)
        SBE = SBE + trace(Rho{a}*kron(eye(q),ME{e,A(a,b)+1}));
end
end
end
SBE = real(SBE)/(AS*n^2);
SB=real(SB)/(AS*n);
SE=real(SE)/(AS*n);
S = (alp*SB+(1-alp)*SE);

