% n # no. of Bob's input or # of Alice's dits
% d dimension of Bob's output.
% q dimension of quantum message.
% mode = 0 for preparaiton equivalence, mode = 1 for measurement equivalence
% OC = 0 for Spekkens' parity conditions, OC = 1 for even parity oblivious
% conditions
function [Pwin,Rho,M] = nto1dpoRAC_Q(n,d,q,mode,OC)
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
Fs = []; % constraints for states.
Fm = []; % constraints for measurements.
%-------------------------------------------------------------------------------
% Preparing sdp variables and constraints.
Rho_sdp = cell(AS,1);
Rho = cell(AS,1);
M_sdp = cell(n,d);
M = cell(n,d)
% constraints for positivity and unit trace of the states.
for a = 1:AS
    Rho_sdp{a} = sdpvar(q,q,'hermitian','complex');
    Fs =  [Fs; Rho_sdp{a} >= 0; trace(Rho_sdp{a}) == 1;];
end
% oblivious constraints.
if (OC ==1)
    for p = 1:(d^(n-1)-1)
    sum1 = 0;
    sum2 = 0;
    for o = 1:d
        sum1 = sum1 + Rho_sdp{oblivious(o,p)};
        sum2 = sum2 + Rho_sdp{oblivious(o,p+1)};
    end
    Fs = [Fs; sum1 == sum2];
end
else
    for p = 1:np
        for vp = 1:d-1
            sum1 = 0;
            sum2 = 0;
            [rowO,colO] = size(oblivious{p});
            for o = 1:rowO
                sum1 = sum1 + Rho_sdp{oblivious{p}(o,vp)};
                sum2 = sum2 + Rho_sdp{oblivious{p}(o,vp+1)};
            end;
            Fs = [Fs; sum1 == sum2];
        end
    end
end
% constraints for positivity and summing upto identity on the measurements.
for i = 1:n
    sum = 0;
    for j = 1:d
        M_sdp{i,j} = sdpvar(q,q,'hermitian','complex');
        Fm = [Fm; M_sdp{i,j} >= 0;];
        sum = sum + M_sdp{i,j};
    end
    Fm = [Fm; sum == eye(q);];
end
% constraints for equivalence of measurement effects.
if (mode == 1)
    for j=1:d-1
        sum1  = 0;
        sum2 = 0;
        for i = 1:n
            sum1 = sum1 + M_sdp{i,j}
            sum2 = sum2 + M_sdp{i,j+1}
        end
        Fm = [Fm; sum1 == sum2;];
    end
end

%-------------------------------------------------------------------------------
% Generating random measurements for the state sdp.
for i = 1:n
    R = RandomPOVM(q,d);
    for j = 1:d
        M{i,j} = R{j};
    end
end
vstepM = 100;
vstepS = 0;
%-------------------------------------------------------------------------------
while (abs(vstepM-vstepS) >= 0.000001)
    % sdp for states.
    stepS = Success(n,d,A,Rho_sdp,M);
    diagnostics = optimize([Fs;], -stepS, sdpsettings('solver', 'mosek','verbose','0'));
    % preparing states for the measurement sdp.
    for i = 1:AS
        Rho{i} = value(Rho_sdp{i});
    end
    vstepS = value(stepS)
    % sdp for measurements.
    stepM = Success(n,d,A,Rho,M_sdp);
    diagnostics = optimize([Fm;], -stepM, sdpsettings('solver', 'sdpt3','verbose','0'));
    for i = 1:n
        for j = 1:d
            M{i,j} = value(M_sdp{i,j});
        end
    end
    vstepM = value(stepM);
end
Pwin = vstepM;
%-------------------------------------------------------------------------------




