% n # no. of Bob's input or # of Alice's dits
% d dimension of Bob's and Eve's output.
% mode = 0 for unconstrained maximization, mode = 1 for maximization given
% Bob's success is fixed to Pb, mode = 2 for security analyiis when success = Pb
% OC = 0 for Spekkens' parity conditions, OC = 1 for even parity oblivious
% conditions
function Pwin = nto1dpoRACBE_Q1BE(n,d,mode,Pb,OC)
AS = d^n; % total # of Alice's settings.
%-------------------------------------------------------------------------------
AS = d^n; % total # of Alice's settings.
%-------------------------------------------------------------------------------
% preparing parity oblivious conditions
[oblivious,A] = parityConSpekkens(AS,n,d);
[np,col]=size(oblivious);
if (OC ==1)
    [oblivious,A] = parityConEven(AS,n,d);
    [np,col]=size(oblivious);
end
%-------------------------------------------------------------------------------
% contraint vector declaration.
F = [];
%-------------------------------------------------------------------------------
% sdpvar declaration.
Gamma = cell(AS,1);
for a = 1:AS
    Gamma{a}=sdpvar(1 + 2*n*(d-1) + n*n*(d-1)*(d-1),1 + 2*n*(d-1) + n*n*(d-1)*(d-1),'hermitian','complex');
    F=[F;Gamma{a}>=0;Gamma{a}(1,1) == 1];
end
%-------------------------------------------------------------------------------
idxB = @(b, B) 1 + (d -1)* b + B + 1;
idxE = @(e, E) 1 + n*(d-1) + (d -1)* e + E + 1;
idxBE = @(b,B,e,E) 1 + 2*n*(d-1) + n*((d-1)^2)*b + ((d-1)^2)*e + (d-1)*B + E + 1;
idxT = @(t, T) 1 + (d-1) * t+ T + 1;
    %-------------------------------------------------------------------------------
    for a = 1:AS
        % constraints for orthogonality of results for the same measurement.
        for tr = 0:(2*n-1)
            for tc = 0:(2*n-1)
                for Tr = 0:(d-2)
                    for Tc = 0:(d-2)
                        r = idxT(tr, Tr);
                        c = idxT(tc, Tc);
                        if r > c
                            if (tr == tc) && (Tr ~= Tc)
                                F = [F; Gamma{a}(r, c) == 0];
                            end
                        end
                    end
                end
                %-------------------------------------------------------------------------------
            end
        end
    %-------------------------------------------------------------------------------
 % constraints for commutation of Bob and Eve's operators.
    for b = 0:(n-1)
        for B = 0:(d-2)
            F = [F; Gamma{a}(idxB(b,B), 1) == Gamma{a}(idxB(b,B), idxB(b,B))];
            F = [F; Gamma{a}(idxE(b,B), 1) == Gamma{a}(idxE(b,B), idxE(b,B))];
            for e = 0:(n-1)
                for E = 0:(d-2)
                    F = [F; Gamma{a}(idxB(b,B),idxE(e,E)) == Gamma{a}(idxE(e,E),idxB(b,B)) ];
                end
            end
        end
    end
%-------------------------------------------------------------------------------    
    for b = 0:(n-1)
        for B = 0:(d-2)
            for e = 0:(n-1)
                for E = 0:(d-2)
                    F=[F; Gamma{a}(idxBE(b,B,e,E), 1) == Gamma{a}(idxB(b,B), idxE(e,E))];
                    F=[F; Gamma{a}(idxBE(b,B,e,E), 1)==Gamma{a}(idxBE(b,B,e,E),idxB(b,B));Gamma{a}(idxBE(b,B,e,E),idxB(b,B))==Gamma{a}(idxB(b,B),idxBE(b,B,e,E))];
                    F=[F; Gamma{a}(idxBE(b,B,e,E), 1)==Gamma{a}(idxBE(b,B,e,E),idxE(e,E));Gamma{a}(idxBE(b,B,e,E),idxE(e,E))==Gamma{a}(idxE(e,E),idxBE(b,B,e,E))];
                    F=[F;Gamma{a}(idxBE(b,B,e,E), idxBE(b,B,e,E))==Gamma{a}(idxB(b,B),idxE(e,E))];
                    for x = 0:(n-1)
                        for X = 0:(d-2)
                            if (x ~= b)
                                F=[F; Gamma{a}(idxBE(b,B,e,E), idxB(x,X)) == Gamma{a}(idxBE(x,X,e,E), idxB(b,B));Gamma{a}( idxB(x,X),idxBE(b,B,e,E)) == Gamma{a}(idxB(b,B),idxBE(x,X,e,E))];
                                F=[F; Gamma{a}(idxBE(b,B,e,E), idxB(x,X)) == Gamma{a}(idxB(b,B),idxBE(x,X,e,E));Gamma{a}(idxB(x,X),idxBE(b,B,e,E)) == Gamma{a}(idxBE(x,X,e,E),idxB(b,B))];
                            end
                            if (x ~= e)
                                F=[F;  Gamma{a}(idxBE(b,B,e,E), idxE(x,X)) == Gamma{a}(idxBE(b,B,x,X), idxE(e,E));Gamma{a}(idxE(x,X),idxBE(b,B,e,E)) == Gamma{a}(idxE(e,E),idxBE(b,B,x,X))];
                                F=[F;  Gamma{a}(idxBE(b,B,e,E), idxE(x,X)) == Gamma{a}(idxE(e,E),idxBE(b,B,x,X));Gamma{a}(idxE(x,X),idxBE(b,B,e,E)) == Gamma{a}(idxBE(b,B,x,X),idxE(e,E))];
                            end
                            for y = 0:(n-1)
                                for Y = 0:(d-2)  
                                    F=[F; Gamma{a}(idxBE(b,B,x,X), idxBE(b,B,y,Y)) == Gamma{a}(idxBE(b,B,x,X),idxE(y,Y)); Gamma{a}(idxBE(x,X,e,E), idxBE(y,Y,e,E)) == Gamma{a}(idxBE(x,X,e,E),idxB(y,Y))];
                                end
                            end
                        end
                    end
                end
            end
        end
    end
%-------------------------------------------------------------------------------    
end
%-------------------------------------------------------------------------------
% Oblivious constraints.
if (OC ==1)
    for p = 1:(d^(n-1)-1)
        sum1 = 0;
        sum2 = 0;
        for o = 1:d
            sum1 = sum1 + Gamma{oblivious(o,p)};
            sum2 = sum2 + Gamma{oblivious(o,p+1)};
        end
        F = [F ; sum1 == sum2];
    end
else
    for p = 1:np
        for vp = 1:d-1
            sum1 = 0;
            sum2 = 0;
            [rowO,colO] = size(oblivious{p});
            for o = 1:rowO
                sum1 = sum1 + Gamma{oblivious{p}(o,vp)};
                sum2 = sum2 + Gamma{oblivious{p}(o,vp+1)};
            end;
            F = [F; sum1 == sum2];
        end
    end
end
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% Making a probability cell.
Pcell = cell(AS,n,n,d-1,d-1);
for a = 1:AS
    for b = 0:(n-1)
        for e = 0:(n-1)
            for B = 0:(d-2)
                for E = 0:(d-2)
                    Pcell{a,b+1,e+1,B+1,E +1 } = Gamma{a}(idxB(b,B), idxE(e,E));
                end
            end
        end
    end
end
%-------------------------------------------------------------------------------
% Making a marginal probability cell for Bob.
PcellB = cell(AS,n,d);
for a = 1:AS
    for b = 0:(n-1)
        sum = 0;
        for B = 0:(d-2)
            PcellB{a,b+1,B+1} = Gamma{a}(idxB(b,B), idxB(b,B));
            sum = sum+ Gamma{a}(idxB(b,B), idxB(b,B));
        end
        PcellB{a,b+1,d} = 1-sum;
    end
end
%-------------------------------------------------------------------------------
% Making a marginal probability cell for Eve.
PcellE = cell(AS,n,d);
for a = 1:AS
    for e = 0:(n-1)
        sum =0;
        for E = 0:(d-2)
            PcellE{a,e+1,E+1 } = Gamma{a}(idxE(e,E), idxE(e,E));
            sum = sum + Gamma{a}(idxE(e,E), idxE(e,E));
        end
        PcellE{a,e+1,d } = 1-sum;
        
    end
end

%-------------------------------------------------------------------------------
% Summing winning probability for Bob
PwinB = 0;
for a = 1:AS
    for b = 0:(n-1)
        PwinB = PwinB + PcellB{a, b+1,A(a,b+1)+1};
    end
end
PwinB=PwinB/(AS*n);
%-------------------------------------------------------------------------------
% Summing winning probability for Eve
PwinE = 0;
for a = 1:AS
    for e = 0:(n-1)
        PwinE = PwinE + PcellE{a,e+1,A(a,e+1)+1};
    end
end
PwinE=PwinE/(AS*n);
%-------------------------------------------------------------------------------
% Summing winning probability for Eve trying to learn Bob's guess upon
% recieving Bob's setting b
PwinEB=0;
for a = 1:AS
    for e = 0:(n-1)
        for b = 0:(n-1)
            PwinEB = PwinEB + PcellE{a,b+1,A(a,b+1)+1};
        end
    end
end
PwinEB=real(PwinEB)/(AS*n^2);
%-------------------------------------------------------------------------------
succ = Pcell{1,1,1,1,1} + Pcell{1,1,2,1,1} +Pcell{1,2,1,1,1} -Pcell{1,2,2,1,1} - PcellB{1,1,1} - PcellE{1,1,1};
%-------------------------------------------------------------------------------
% optimization
Pwin = real(PwinB)
if (mode == 1)
    F = [F; PwinB == Pb];
    Pwin = real(PwinE);
end
if (mode == 2)
    F = [F; PwinB == Pb];
    Pwin = real(PwinEB);
end
diagnostics = optimize([F], -real(Pwin), sdpsettings('solver', 'mosek','verbose',0))
%--------------------------------------------------------------------------------
Pwin = value(Pwin);
