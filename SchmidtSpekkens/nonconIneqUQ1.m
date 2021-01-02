function S = nonconIneqUQ1()
X = 6; % six settings for the preparation device
Y = 3; % three settings for the measurement device
K = 2; % binary outcomes
O = 2*Y*(K-1)+1; % total number of operators in our list
conP = []; % an empty list for constraints on the moment matrix level
conM = []; % an empty list for constraints on the substrate level
Prob = cell(X,Y,K-1); % a cell for probabilities
S = []; % an empty list for upper bounds
G = cell(X,1); % cell for X moment matrices


for x = 0:X-1
    G{x+1}=sdpvar(O,O,'hermitian','complex'); % declaration of our SDP variables
    conP = [conP;G{x+1} >= 0]; % semi-definite constraints
end
conP = [conP; G{1} + G{2} == G{3} + G{4}; G{1} + G{2} == G{5} + G{6}]; % preparation equivalences 


idx = @(y, k, u) 2*(K-1)*y + 2*k + u + 2; % function to return the position of the operators 



for x = 0:X-1
    for y = 0:Y-1
        for k = 0:K-2 
            conM = [conM; G{x+1}(1,idx(y,k,0)) == G{x+1}(idx(y,k,1),1)]; 
            conM = [conM; G{x+1}(idx(y,k,0),1) == G{x+1}(1,idx(y,k,1))];
        end
    end
    for j = 1:O
        conM = [conM; G{x+1}(j,j) == 1]; % unitarity constraints
    end
end

for x = 0:X-1
    for y = 0:Y-1
        for k = 0:K-2
            Prob{x+1,y+1,k+1} = 0.5 + 0.25 * (G{x+1}(1, idx(y,k,0)) + G{x+1}(1, idx(y,k,1)));
        end
    end
end



for x = 0:X-1
    for j = 1:O
        sum1 = 0; sum2 = 0; 
        for y = 0:Y-1
            for k = 0:K-2 
                sum1 = sum1 + G{x+1}(j,idx(y,k,0)) + G{x+1}(j,idx(y,k,1));
                sum2 = sum2 + G{x+1}(idx(y,k,0),j) + G{x+1}(idx(y,k,1),j);
            end
        end
        conM = [conM; sum1 == 0; sum2 == 0]; % measurement equivalences 
    end
end

S1 = real(Prob{1,1,1} + Prob{3,2,1} + Prob{5,3,1});
diagnostics = optimize([conP;conM], -S1, sdpsettings('solver', 'sdpt3'));
S = [S;value(S1)];
%-------------------------------------------------------------------------------  
S2 = real(Prob{1,1,1} + Prob{2,2,1} + Prob{5,3,1});
diagnostics = optimize([conP;conM], -S2, sdpsettings('solver', 'sdpt3'));
S = [S;value(S2)];
%------------------------------------------------------------------------------- 
S3 = real(Prob{1,1,1} - Prob{3,1,1} -2 * Prob{5,1,1} -2 * Prob{2,2,1} + 2 * Prob{3,2,1} + 2 * Prob{5,3,1});
diagnostics = optimize([conP;conM], -S3, sdpsettings('solver', 'sdpt3'));
S = [S;value(S3)];
%------------------------------------------------------------------------------- 
S4 = real(2* Prob{1,1,1} - Prob{2,2,1} +2* Prob{3,2,1}); 
diagnostics = optimize([conP;conM], -S4, sdpsettings('solver', 'sdpt3'));
S = [S;value(S4)];
%------------------------------------------------------------------------------- 
S5 = real(Prob{1,1,1} - Prob{5,1,1} +  Prob{2,2,1} + Prob{3,2,1} + 2 * Prob{5,3,1});
diagnostics = optimize([conP;conM], -S5, sdpsettings('solver', 'sdpt3'));
S = [S;value(S5)];
%------------------------------------------------------------------------------- 
S6 = real(Prob{1,1,1} - Prob{5,1,1} +  2*Prob{2,2,1} + 2 * Prob{5,3,1});
diagnostics = optimize([conP;conM], -S6, sdpsettings('solver', 'sdpt3'));
S = [S;value(S6)];
%------------------------------------------------------------------------------- 
S7 = real(Prob{1,1,1} - Prob{4,1,1} -2 * Prob{5,1,1} -2 * Prob{2,2,1} + 2 * Prob{3,2,1} + 2 * Prob{5,3,1});
diagnostics = optimize([conP;conM], -S7, sdpsettings('solver', 'sdpt3'));
S = [S;value(S7)];
