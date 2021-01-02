function [Pwin,Rho,M] =nonconIneqQ(q,ineq)
AS = 6;
n = 3;
d=2;
iter=20;
Fs = []; % constraints for states.
Fm = []; % constraints for measurements.
% Preparing sdp variables and constraints.
Rho_sdp = cell(AS,1);
Rho = cell(AS,1);
M_sdp = cell(n,d);
M = cell(n,d)
% constraints for positivity and unit trace of the states.
for a = 1:AS+1
    Rho_sdp{a} = sdpvar(q,q,'hermitian','complex');
    Fs =  [Fs; Rho_sdp{a} >= 0; trace(Rho_sdp{a}) == 1;];
end
%Fs = [Fs; Rho_sdp{2}+Rho_sdp{1}== Rho_sdp{3} +Rho_sdp{4}; Rho_sdp{3} +Rho_sdp{4} == Rho_sdp{5} +Rho_sdp{6}];
Fs = [Fs; Rho_sdp{2} == 2*Rho_sdp{7} - Rho_sdp{1}];
Fs = [Fs; Rho_sdp{4} == 2*Rho_sdp{7} - Rho_sdp{3}];
Fs = [Fs; Rho_sdp{6} == 2*Rho_sdp{7} - Rho_sdp{5}];
for i = 1:n
    sum = 0;
    for j = 1:d
        M_sdp{i,j} = sdpvar(q,q,'hermitian','complex');
        Fm = [Fm; M_sdp{i,j} >= 0;];
        sum = sum + M_sdp{i,j};
    end
    Fm = [Fm; sum == eye(q);];
end
Fm = [Fm;M_sdp{1,1}+M_sdp{2,1}+M_sdp{3,1} ==M_sdp{1,2}+M_sdp{2,2}+M_sdp{3,2}   ];
%-------------------------------------------------------------------------------
% Generating random measurements for the state sdp.
for i = 1:n
    R = RandomPOVM(q,d);
    for j = 1:d
        M{i,j} = R{j};
    end
end
stepM=100;stepS=0;
%-------------------------------------------------------------------------------
while abs(value(stepM)-value(stepS))>=0.000001
% sdp for states.    
    stepS = SuccNC(Rho_sdp,M,ineq);
    diagnostics = optimize([Fs;], -stepS, sdpsettings('solver', 'sdpt3','verbose','0'));
% preparing states for the measurement sdp. 
    for i = 1:AS+1
        Rho{i} = value(Rho_sdp{i});
    end
% sdp for measurements.
    stepM = SuccNC(Rho,M_sdp,ineq);
    diagnostics = optimize([Fm;], -stepM, sdpsettings('solver', 'sdpt3','verbose','0'));
    for i = 1:n
        for j = 1:d
            M{i,j} = value(M_sdp{i,j});
        end
    end
end
%-------------------------------------------------------------------------------
Pwin = value(stepM);
