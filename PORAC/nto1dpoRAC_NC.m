% n # no. of Bob's input or # of Alice's dits
% d dimension of Bob's output.
% mode = 0 for preparaiton equivalence, mode = 1 for measurement equivalence
% OC = 0 for Spekkens' parity conditions, OC = 1 for even parity oblivious
% conditions
function [Pwin,Pe] = nto1dpoRAC_NC(n,d,mode,OC)
AS = d^n; % total # of Alice's settings.
iter = 20;
%-------------------------------------------------------------------------------
% preparing parity oblivious conditions
[oblivious,A] = parityConSpekkens(AS,n,d);
[np,col]=size(oblivious);
if (OC ==1)
    [oblivious,A] = parityConEven(AS,n,d);
    [np,col]=size(oblivious);
end
%-------------------------------------------------------------------------------
Fe=[];
Fd=[];
Pe = cell(AS,1);
Pe_LP  = cell(AS,1);
Pd = cell(n,AS);
Pd_LP = cell(AS,n);
for a = 1:AS;
    Pe_LP{a} = sdpvar(AS,1);
    Fe = [Fe;Pe_LP{a}>=0];
    Sum = 0;
    for m = 1:AS
        Sum = Sum + Pe_LP{a}(m);
    end
    Fe = [Fe;Sum==1];
end
for m = 1:AS
    if (OC == 1)
        for p = 1:(d^(n-1)-1)
            sum1 = 0;
            sum2 = 0;
            for o = 1:d
                sum1 = sum1 + Pe_LP{oblivious(o,p)}(m);
                sum2 = sum2 + Pe_LP{oblivious(o,p+1)}(m);
            end
            Fe = [Fe; sum1 == sum2];
        end
    else
        for p = 1:np
            for vp = 1:d-1
                sum1 = 0;
                sum2 = 0;
                [rowO,colO] = size(oblivious{p});
                for o = 1:rowO
                    sum1 = sum1 + Pe_LP{oblivious{p}(o,vp)}(m);
                    sum2 = sum2 + Pe_LP{oblivious{p}(o,vp+1)}(m);
                end
                Fe = [Fe; sum1 == sum2];
            end
        end
    end
end

for m = 1:AS
    for b = 0: n-1
        Pd_LP{m,b+1} = sdpvar(d,1);
        Fd= [Fd;Pd_LP{m,b+1}>=0];
        r=rand(d,1);  r= r/sum(r);
        Pd{m,b+1} = r;
        Sum =0;
        for B = 0:d-1
            Sum = Sum  + Pd_LP{m,b+1}(B+1);
        end
        Fd = [Fd; Sum == 1];
    end
    if (mode ==1)
        for j=1:d-1
            sum1  = 0;
            sum2 = 0;
            for i = 1:n
                sum1 = sum1 + Pd_LP{m,i}(j)
                sum2 = sum2 + Pd_LP{m,i}(j+1)
            end
            Fd = [Fd; sum1 == sum2;];
        end
    end
end
stepM=100;stepS=0;
while (abs(value(stepM-stepS)) >= 0.00001)
    % sdp for states.
    stepS = SuccC(n,d,A,Pe_LP,Pd);
    diagnostics = optimize([Fe;], -stepS, sdpsettings('solver', 'sdpt3','verbose','0'));
    % preparing states for the measurement sdp.
    for i = 1:AS
        Pe{i} = value(Pe_LP{i});
    end
    % sdp for measurements.
    stepM = SuccC(n,d,A,Pe,Pd_LP);
    diagnostics = optimize([Fd;], -stepM, sdpsettings('solver', 'sdpt3','verbose','0'));
    for i = 1:n
        for j = 1:AS
            Pd{j,i} = value(Pd_LP{j,i});
        end
    end
end
Pwin = value(stepM);