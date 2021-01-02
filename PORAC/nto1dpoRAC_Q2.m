% n # no. of Bob's input or # of Alice's dits
% d dimension of Bob's output.
% mode = 0 for preparaiton equivalence, mode = 1 for measurement equivalence
% OC = 0 for Spekkens' parity conditions, OC = 1 for even parity oblivious
% conditions
function Pwin = nto1dpoRAC_Q2(n,d,mode,OC)
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
F4 = []; % Q2 constraints
F5 = []; % oblivious!
F6 = []; % constraints with identity
%-------------------------------------------------------------------------------
% sdpvar declaration.
Gamma = cell(AS,1);
for a = 1:AS
    Gamma{a}=sdpvar(1 + n*d + n*n*(d-1)*(d-1),1 + n*d + n*n*(d-1)*(d-1),'hermitian','complex');
    F0=[F0;Gamma{a}>=0];
end
%-------------------------------------------------------------------------------
idx = @(b, B) 1 + d * b + B + 1; % function to return the position of the operator corresponding to Bob's setting (b,B) in a Gamma matrix.
idxBB = @(b1, b2, B1, B2) 1 + n * d + b1 * n*(d-1)*(d-1) + b2*(d-1)*(d-1) + B1*(d-1) + B2 + 1; % # function to return the position of the operator corresponding to Bob's setting (b1, b2, B1, B2) in a Q2 Gamma matrix.
%-------------------------------------------------------------------------------
% General Matrix constraints F1,F2,F3
for a = 1:AS
    for br = 0:(n-1)
        for bc = 0:(n-1)
            %-------------------------------------------------------------------------------
            % constraints for orthogonality of results for the same measurement.
            for Br = 0:(d-1)
                for Bc = 0:(d-1)
                    r = idx(br, Br);
                    c = idx(bc, Bc);
                    if r > c
                        if (br == bc) && (Br ~= Bc)
                            F1 = [F1; Gamma{a}(r, c) == 0];
                        end
                    end
                end
            end
            %-------------------------------------------------------------------------------
            % constraints for summation over results of one measurement (ROW).
            for Br = 0:(d-1)
                PaabbR = 0;
                for Bc = 0:(d-1)
                    r = idx(br, Br);
                    c = idx(bc, Bc);
                    PaabbR = PaabbR + Gamma{a}(r, c);
                end
                F2 = [F2; PaabbR == Gamma{a}(idx(br, Br), idx(br, Br))];
            end
            %-------------------------------------------------------------------------------
            % constraints for summation over results of one measurement (COL).
            for Bc = 0:(d-1)
                PaabbC = 0;
                for Br = 0:(d-1)
                    r = idx(br, Br);
                    c = idx(bc, Bc);
                    PaabbC = PaabbC + Gamma{a}(r, c);
                end
                F2 = [F2; PaabbC == Gamma{a}(idx(bc, Bc), idx(bc, Bc))];
            end
            %-------------------------------------------------------------------------------
            % constraints for summation over results of all measurements, and normalization.
            PsiPsi = 0;
            for Bc = 0:(d-1)
                for Br = 0:(d-1)
                    r = idx(br, Br);
                    c = idx(bc, Bc);
                    PsiPsi = PsiPsi + Gamma{a}(r, c);
                end
            end
            F3 = [F3; PsiPsi == 1];
            %-------------------------------------------------------------------------------
            % constraints for double operator BB terms.
            for B1 = 0:(d-2)
                for B2 = 0:(d-2)
                    for b1 = 0:(n-1)
                        for b2 = 0:(n-1)
                            F4 = [F4; Gamma{a}(idx(b1,B1), idx(b2,B2)) == Gamma{a}(1, idxBB(b1, b2, B1, B2))];
                            for B3 = 0:(d-2)
                                for b3 = 0:(n-1)
                                    F4 = [F4; Gamma{a}(idxBB(b1, b2, B1, B2), idxBB(b1, b3, B1, B3)) == Gamma{a}(idxBB(b1, b2, B1, B2), idx(b3, B3))];
                                end
                            end
                        end
                    end
                end
            end
            %-------------------------------------------------------------------------------
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
if (mode == 1)
    for a = 1:AS
        for b = 0:(n-1)
            for B = 0:(d-1)
                for B1 = 0:(d-2)
                    sum1=0; sum2=0; sum3=0;sum4=0;
                    for b1 = 0:(n-1)
                        sum1 = sum1 + Gamma{a}(idx(b, B), idx(b1, B1));
                        sum2 = sum2 + Gamma{a}(idx(b, B), idx(b1, B1+1));
                        sum3 = sum3 + Gamma{a}(idx(b1, B1),idx(b, B));
                        sum4 = sum4 + Gamma{a}(idx(b1, B1+1),idx(b, B));
                    end
                    F5 = [F5; sum1==sum2;sum3==sum4];
                end
            end
        end
    end
    for a = 1:AS
        for B1 = 0:(d-2)
            for B2 = 0:(d-2)
                for b1 = 0:(n-1)
                    for b2 = 0:(n-1)
                        for B3 = 0:(d-2)
                            sum1=0; sum2=0; sum3=0;sum4=0;
                            for b3 = 0:(n-1)
                                sum1 = sum1 + Gamma{a}(idxBB(b1, b2, B1, B2), idx(b3, B3));
                                sum2 = sum2 + Gamma{a}(idxBB(b1, b2, B1, B2), idx(b3, B3+1));
                                sum3 = sum3 + Gamma{a}(idx(b3, B3),idxBB(b1, b2, B1, B2));
                                sum4 = sum4 + Gamma{a}(idx(b3, B3+1),idxBB(b1, b2, B1, B2));
                            end
                            F5 = [F5; sum1==sum2;sum3==sum4];
                        end
                    end
                end
            end
        end
    end
end
%-------------------------------------------------------------------------------
% Making a probability cell.
Pcell = cell(AS,n,d); % a0, a1, b, B
for a = 1:AS
    F6 = [F6; Gamma{a}(1,1) == 1];
    for b = 0:(n-1)
        for B = 0:(d-1)
            Pcell{a,b+1,B+1 } = Gamma{a}(idx(b,B), idx(b,B));
            F6 = [F6; Gamma{a}(idx(b,B), 1) == Gamma{a}(idx(b,B), idx(b,B))];
        end
    end
end
%-------------------------------------------------------------------------------
% Summing winning probability
Pwin = 0;
for a = 1:AS
    for b = 0:(n-1)
        Pwin = Pwin + Pcell{a, b+1, A(a,b+1)+1};
    end
end
Pwin = Pwin/(AS*n);
%-------------------------------------------------------------------------------
% optimization
diagnostics = optimize([F0;F1;F2;F3;F4;F5;F6], -Pwin, sdpsettings('solver', 'mosek'))
%--------------------------------------------------------------------------------
Pwin = value(Pwin);