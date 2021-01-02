function S = SuccNC(Rho,M,ineq);
S1=trace(Rho{1}*M{1,1}+Rho{3}*M{2,1}+Rho{5}*M{3,1});
S2=trace(Rho{1}*M{1,1}+Rho{2}*M{2,1}+Rho{5}*M{3,1});
S3 = trace(Rho{1}*M{1,1}+2*Rho{3}*M{2,1} + 2* Rho{5}*M{3,1} -2*Rho{2}*M{2,1}  -2*Rho{5}*M{1,1}-Rho{3}*M{1,1});
S4=trace(2*Rho{1}*M{1,1}+2*Rho{3}*M{2,1}-Rho{2}*M{2,1});
S5=trace(Rho{1}*M{1,1}+Rho{2}*M{2,1}+Rho{3}*M{2,1}+2*Rho{5}*M{3,1}-Rho{5}*M{1,1}); 
S6=trace(Rho{1}*M{1,1}+2*Rho{2}*M{2,1}+2*Rho{5}*M{3,1}-Rho{5}*M{1,1});
S7= trace(Rho{1}*M{1,1}+2*Rho{3}*M{2,1} + 2* Rho{5}*M{3,1} -2*Rho{2}*M{2,1}  -2*Rho{5}*M{1,1}-Rho{4}*M{1,1});

switch ineq
    case 1
        S = real(S1);
        case 2
        S = real(S2);
        case 3
        S = real(S3);
        case 4
        S = real(S4);
        case 5
        S = real(S5);
        case 6
        S = real(S6);
         case 7
        S = real(S7);
end

