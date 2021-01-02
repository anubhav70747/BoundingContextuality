% n # no. of Bob's input or # of Alice's dits
% d dimension of Bob's and Eve's output.
% mode = 0 for unconstrained maximization, mode = 1 for maximization given
% Bob's success is fixed to Pb, mode = 2 for security analyiis when success = Pb
% OC = 0 for Spekkens' parity conditions, OC = 1 for even parity oblivious
% conditions
function Pwin = nto1dpoRACBE_OC(n,d,mode,Pb,OC)
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
% function to return the position corresponding to the probability P(E,B|e,b,AS) in the probability vector.
idx = @(AS,b,e,B,E) (n^2 * d^2) * (AS-1) + (n * d^2) * b + (d^2 * e) + (d * B) + E+1
%-------------------------------------------------------------------------------
x=sdpvar( d ^ n * (d ^ 2) * (n ^ 2), 1); % the probability vector as the variable for lin prog.
%-------------------------------------------------------------------------------
con = [0<=x<=1];  %positivity constraint
%-------------------------------------------------------------------------------
for a = 1:AS
    for b = 0:n-1
        for e = 0:n-1
            %-------------------------------------------------------------------------------
            sum = 0;
            for B = 0:d-1
                for E = 0:d-1
                    r=idx(a,b,e,B,E);
                    sum = sum + x(r);
                end
            end
            con = [con;sum==1]; % completeness constraint
            
        end
    end
end
%-------------------------------------------------------------------------------
if (OC == 1)
    for b = 0:n-1
        for e = 0:n-1
            for B = 0:d-1
                for E = 0:d-1
                    for p = 1:(d^(n-1)-1)
                        sum1 = 0;
                        sum2 = 0;
                        for o = 1:d
                            sum1 = sum1 + x(idx(oblivious(o,p),b,e,B,E));
                            sum2 = sum2 + x(idx(oblivious(o,p+1),b,e,B,E));
                        end
                        con = [con ; sum1 == sum2;]; % operational oblivious constraint
                    end
                end
            end
        end
    end
else
    for b = 0:n-1
        for e = 0:n-1
            for B = 0:d-1
                for E = 0:d-1
                    for p = 1:np
                        for vp = 1:d-1
                            sum1 = 0;
                            sum2 = 0;
                            [rowO,colO] = size(oblivious{p});
                            for o = 1:rowO
                                sum1 = sum1 + x(idx(oblivious{p}(o,vp),b,e,B,E));
                                sum2 = sum2 + x(idx(oblivious{p}(o,vp+1),b,e,B,E));
                            end;
                            con = [con; sum1 == sum2];
                        end
                    end
                end
            end
        end
    end
end
%-------------------------------------------------------------------------------
for a = 1:AS
    for b = 0:n-1
        for B = 0:d-1
            for e = 0:n-2
                sum1=0;
                sum2=0;
                for E = 0:d-1
                    sum1 = sum1 +  x(idx(a,b,e,B,E));
                    sum2 = sum2 +  x(idx(a,b,e+1,B,E));
                end
                con =[con;  sum1==sum2]; % no-signaling from Eve to Bob
                sum1=0;
                sum2=0;
                for E = 0:d-1
                    sum1 = sum1 +  x(idx(a,e,b,E,B));
                    sum2 = sum2 +  x(idx(a,e+1,b,E,B));
                end
                con =[con;  sum1==sum2]; % no sigaling from Bob to Eve
            end
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
                    Pcell{a,b+1,e+1,B+1,E+1 } = x(idx(a,b,e,B,E));
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
            sum = 0;
            for e = 0:(n-1)
                for E = 0:(d-1)
                    sum = sum + x(idx(a,b,e,B,E));
                end
            end
            PcellB{a,b+1,B+1} = sum/n;
        end
    end
end
%-------------------------------------------------------------------------------
% Making a marginal probability cell for Eve.
PcellE = cell(AS,n,d); % a0, a1, b, B
for a = 1:AS
    for e = 0:(n-1)
        for E = 0:(d-1)
            sum = 0;
            for b = 0:(n-1)
                for B = 0:(d-1)
                    sum = sum + x(idx(a,b,e,B,E));
                end
            end
            PcellE{a,e+1,E+1 } = sum/n;
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
PwinB=PwinB/(d^n*n);
%-------------------------------------------------------------------------------
% Summing winning probability for Eve
PwinE = 0;
for a = 1:AS
    for e = 0:(n-1)
        PwinE = PwinE + PcellE{a,e+1,A(a,e+1)+1};
    end
end
PwinE=PwinE/(d^n*n);
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
PwinEB=PwinEB/(d^n*n^2);
%-------------------------------------------------------------------------------
% optimization
Pwin = real(PwinEB+PwinB)
if (mode == 1)
    con = [con; PwinB == Pb];
    Pwin = real(PwinE);
end
if (mode == 2)
    con = [con; PwinB == Pb];
    Pwin = real(PwinEB);
end
diagnostics = optimize(con,-Pwin, sdpsettings('solver', 'sdpt3'))
Pwin = value(Pwin);
