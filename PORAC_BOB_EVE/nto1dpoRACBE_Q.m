% n # no. of Bob's input or # of Alice's dits,
% d dimension of Bob's and Eve's output,
% OC = 0 for Spekkens' parity conditions, OC = 1 for even parity oblivious
% conditions, alpha is the biasness parameter required for monogamy
% relations and security analysis,
% requires SuccessBE.m for calculation of the objective function
function [Pwin,SB,SE] = nto1dpoRACBE_Q(n,d,q,OC,alpha)
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
FmB = []; % constraints for measurements.
FmE = []; % constraints for measurements.

%-------------------------------------------------------------------------------
% Preparing sdp variables and constraints.
Rho_sdp = cell(AS,1);
Rho = cell(AS,1);
MB_sdp = cell(n,d);
MB = cell(n,d);
ME_sdp = cell(n,d);
ME = cell(n,d);
% constraints for positivity and unit trace of the states.
for a = 1:AS
    Rho_sdp{a} = sdpvar(q^2,q^2,'hermitian','complex');
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
        ME_sdp{i,j} = sdpvar(q,q,'hermitian','complex');
        FmE = [FmE; ME_sdp{i,j} >= 0;];
        sum = sum + ME_sdp{i,j};
    end
    FmE = [FmE; sum == eye(q);];
end
% constraints for equivalence of measurement effects.
for i = 1:n
    sum = 0;
    for j = 1:d
        MB_sdp{i,j} = sdpvar(q,q,'hermitian','complex');
        FmB = [FmB; MB_sdp{i,j} >= 0;];
        sum = sum + MB_sdp{i,j};
    end
    FmB = [FmB; sum == eye(q);];
end

%-------------------------------------------------------------------------------
% Generating random measurements for the state sdp.
for i = 1:n
    R = RandomPOVM(q,d);
    for j = 1:d
        MB{i,j} = R{j};
    end
end
for i = 1:n
    R = RandomPOVM(q,d);
    for j = 1:d
        ME{i,j} = R{j};
    end
end
vstepMB = 100;
vstepME = 100;
vstepS = 0;
%-------------------------------------------------------------------------------
while (abs(vstepME-vstepS) >= 0.000001)
    % sdp for states.
    stepS = SuccessBE(n,d,q,A,Rho_sdp,MB,ME,alpha);
    diagnostics = optimize([Fs;], -stepS, sdpsettings('solver', 'sdpt3','verbose',0));
    % preparing states for the measurement sdp.
    for i = 1:AS
        Rho{i} = value(Rho_sdp{i});
    end
    vstepS = value(stepS);
    
    
    % sdp for measurements.
    stepMB = SuccessBE(n,d,q,A,Rho,MB_sdp,ME,alpha);
    diagnostics = optimize([FmB;], -stepMB, sdpsettings('solver', 'sdpt3','verbose',0));
    for i = 1:n
        for j = 1:d
            MB{i,j} = value(MB_sdp{i,j});
        end
    end
    vstepMB = value(stepMB);
    % sdp for measurements.
    stepME = SuccessBE(n,d,q,A,Rho,MB,ME_sdp,alpha);
    diagnostics = optimize([FmE;], -stepME, sdpsettings('solver', 'sdpt3','verbose',0));
    for i = 1:n
        for j = 1:d
            ME{i,j} = value(ME_sdp{i,j});
        end
    end
    vstepME = value(stepME);
end
Pwin = vstepME;
[S,SB,SE] = SuccessBE(n,d,q,A,Rho,MB,ME,alpha);
%-------------------------------------------------------------------------------




