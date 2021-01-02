% n # no. of Bob's input or # of Alice's dits
% d dimension of Bob's and Eve's output.
% mode = 0 for unconstrained maximization, mode = 1 for maximization given
% Bob's success is fixed to Pb, mode = 2 for security analyiis when success = Pb
% OC = 0 for Spekkens' parity conditions, OC = 1 for even parity oblivious
% conditions
function Pwin = nto1dpoRACBE_Q1(n,d,mode,Pb,OC)
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
F0 = []; % contraint for positivity of Gamma.
F1 = []; % constraints for orthogonality of results for the same measurement.
F2 = []; % constraints for summation over results of one measurement.
F3 = []; % constraints for summation over results of all measurements, and normalization.
F4 = []; % oblivious!
F5 = []; % constraints for commutation of Bob and Eve's operators.
F6 = []; % constraints for positivity of probabilities.
%-------------------------------------------------------------------------------
% sdpvar declaration.
Gamma = cell(AS,1);
for a = 1:AS
    Gamma{a}=sdpvar(2*n*d,2*n*d,'hermitian','complex');
    F0=[F0;Gamma{a}>=0];
end
%-------------------------------------------------------------------------------
idxB = @(b, B) d * b + B + 1; % function to return the position of the operator corresponding to Bob's setting (b,B) in a Gamma matrix.
idxE = @(e, E) n * d + d * e + E + 1; % function to return the position of the operator corresponding to Eve's setting (e,E) in a Gamma matrix.
idxBE = @(be, BE) d * be + BE + 1; %  % function to return the position of the operator corresponding to Bob or Eve's setting (be,BE) in a Gamma matrix.
%-------------------------------------------------------------------------------
% General Matrix constraints F1,F2,F3
for a = 1:AS
    for ber = 0:(2*n-1)
        for bec = 0:(2*n-1)
            %-------------------------------------------------------------------------------
            % constraints for orthogonality of results for the same measurement.
            for BEr = 0:(d-1)
                for BEc = 0:(d-1)
                    r = idxBE(ber, BEr);
                    c = idxBE(bec, BEc);
                    if r > c
                        if (ber == bec) && (BEr ~= BEc)
                            F1 = [F1; Gamma{a}(r, c) == 0];
                        end
                    end
                end
            end
            %-------------------------------------------------------------------------------
            % constraints for summation over results of one measurement (ROW).
            for BEr = 0:(d-1)
                PaabbR = 0;
                for BEc = 0:(d-1)
                    r = idxBE(ber, BEr);
                    c = idxBE(bec, BEc);
                    PaabbR = PaabbR + Gamma{a}(r, c);
                end
                F2 = [F2; PaabbR == Gamma{a}(idxBE(ber, BEr), idxBE(ber, BEr))];
            end
            %-------------------------------------------------------------------------------
            % constraints for summation over results of one measurement (COL).
            for BEc = 0:(d-1)
                PaabbC = 0;
                for BEr = 0:(d-1)
                    r = idxBE(ber, BEr);
                    c = idxBE(bec, BEc);
                    PaabbC = PaabbC + Gamma{a}(r, c);
                end
                F2 = [F2; PaabbC == Gamma{a}(idxBE(bec, BEc), idxBE(bec, BEc))];
            end
            %-------------------------------------------------------------------------------
            % constraints for summation over results of all measurements, and normalization.
            PsiPsi = 0;
            for BEc = 0:(d-1)
                for BEr = 0:(d-1)
                    r = idxBE(ber, BEr);
                    c = idxBE(bec, BEc);
                    PsiPsi = PsiPsi + Gamma{a}(r, c);
                end
            end
            F3 = [F3; PsiPsi == 1];
            %-------------------------------------------------------------------------------
        end
    end
    %-------------------------------------------------------------------------------
    % constraints for commutation of Bob and Eve's operators.
    for b = 0:(n-1)
        for e = 0:(n-1)
            for B = 0:(d-1)
                for E = 0:(d-1)
                    F5 = [F5; Gamma{a}(idxB(b,B),idxE(e,E)) == Gamma{a}(idxE(e,E),idxB(b,B)) ];
                end
            end
        end
    end
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
        F4 = [F4 ; sum1 == sum2];
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
            F4 = [F4; sum1 == sum2];
        end
    end
end
%-------------------------------------------------------------------------------
% Making a probability cell.
Pcell = cell(AS,n,n,d,d);
for a = 1:AS
    for b = 0:(n-1)
        for e = 0:(n-1)
            for B = 0:(d-1)
                for E = 0:(d-1)
                    Pcell{a,b+1,e+1,B+1,E+1 } = Gamma{a}(idxB(b,B), idxE(e,E));
                    F6 = [F6 ; Gamma{a}(idxB(b,B), idxE(e,E))>=0];
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
        for B = 0:(d-1)
            PcellB{a,b+1,B+1} = Gamma{a}(idxB(b,B), idxB(b,B));
            F6 = [F6 ; Gamma{a}(idxB(b,B), idxB(b,B))>=0];
        end
    end
end
%-------------------------------------------------------------------------------
% Making a marginal probability cell for Eve.
PcellE = cell(AS,n,d);
for a = 1:AS
    for e = 0:(n-1)
        for E = 0:(d-1)
            PcellE{a,e+1,E+1 } = Gamma{a}(idxE(e,E), idxE(e,E));
            F6 = [F6 ; Gamma{a}(idxE(e,E), idxE(e,E))>=0];
        end
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
% Summing winning probability for Eve trying to learn Bob's guess
PwinEB=0;
for a = 1:AS
    for e = 0:(n-1)
        for b = 0:(n-1)
            for B = 0:(d-1)
                PwinEB = PwinEB + Pcell{a,b+1,e+1,B+1,A(a,b+1)+1};
            end
        end
    end
end
PwinEB=real(PwinEB)/(AS*n^2);
%-------------------------------------------------------------------------------
% test CHSH just for shits and gigles
CHSH = 0;
count =0;
for b = 0:1
    for e = 0:1
        for B = 0:1
            for E = 0:1
                if (b*e) == xor(B,E)
                    count = count + 1;
                    CHSH = CHSH + Gamma{1}(idxB(b,B), idxE(e,E));
                end
            end
        end
    end
end
CHSH = CHSH/4

%-------------------------------------------------------------------------------

% optimization
Pwin = real(PwinE+PwinB)
if (mode == 1)
    F0 = [F0; PwinB == Pb];
    Pwin = real(PwinE);
end
if (mode == 2)
    Pwin = real(PwinB+PwinEB)
    F0 = [F0; PwinB == Pb];
    Pwin = real(PwinEB);
end
diagnostics = optimize([F0;F1;F2;F3;F4;F5;F6], -real(Pwin), sdpsettings('solver', 'sdpt3'))
Pwin = value(Pwin);
%--------------------------------------------------------------------------------