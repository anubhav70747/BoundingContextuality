% n # no. of Bob's input or # of Alice's dits
% d dimension of Bob's output.
% mode = 0 for preparaiton equivalence, mode = 1 for measurement equivalence
% OC = 0 for Spekkens' parity conditions, OC = 1 for even parity oblivious
% conditions
function Pwin = nto1dpoRAC_OC(n,d,mode,OC)
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
idx = @(AS,b,B) (n * d) * (AS-1) + (d * b) + B+1
%-------------------------------------------------------------------------------
x=sdpvar( d ^ n * d * n, 1); % the probability vector as the variable for lin prog.
%-------------------------------------------------------------------------------
con = [0<=x<=1];  %positivity constraint
%-------------------------------------------------------------------------------
for a = 1:AS
    for b = 0:n-1
        %-------------------------------------------------------------------------------
        sum = 0;
        for B = 0:d-1
            r=idx(a,b,B);
            sum = sum + x(r);
        end
        con = [con;sum==1]; % completeness constraint
    end
end
%-------------------------------------------------------------------------------
if (OC==1)
    for b = 0:n-1
        for B = 0:d-1
            for p = 1:(d^(n-1)-1)
                sum1 = 0;
                sum2 = 0;
                for o = 1:d
                    sum1 = sum1 + x(idx(oblivious(o,p),b,B));
                    sum2 = sum2 + x(idx(oblivious(o,p+1),b,B));
                end
                con = [con ; sum1 == sum2]; % operational oblivious constraint
            end
        end
    end
else
    for b = 0:n-1
        for B = 0:d-1
            for p = 1:np
                for vp = 1:d-1
                    sum1 = 0;
                    sum2 = 0;
                    [rowO,colO] = size(oblivious{p});
                    for o = 1:rowO
                        sum1 = sum1 + x(idx(oblivious{p}(o,vp),b,B));
                        sum2 = sum2 + x(idx(oblivious{p}(o,vp+1),b,B));
                    end;
                    con = [con; sum1 == sum2];
                end
            end
        end
    end
end
%-------------------------------------------------------------------------------
if (mode==1)
    for a = 1:AS
        for B = 0:d-2
            sum1 = 0; sum2 =0;
            for b = 0:n-1
                sum1 = sum1 + x(idx(a,b,B));
                sum2 = sum2 + x(idx(a,b,B+1));
            end
            con = [con; sum1 == sum2]; % operational measurement equivalence
        end
    end
end
%-------------------------------------------------------------------------------
% Making a probability cell.
Pcell = cell(AS,n,d);
for a = 1:AS
    for b = 0:(n-1)
        for B = 0:(d-1)
            Pcell{a,b+1,B+1} = x(idx(a,b,B));
        end
    end
end
%-------------------------------------------------------------------------------
% Summing winning probability for Bob
PwinB = 0;
for a = 1:AS
    for b = 0:(n-1)
        PwinB = PwinB + Pcell{a, b+1,A(a,b+1)+1};
    end
end
PwinB=PwinB/(d^n*n);

%-------------------------------------------------------------------------------
% optimization
diagnostics = optimize(con,-PwinB, sdpsettings('solver', 'sdpt3'))
Pwin = value(PwinB);