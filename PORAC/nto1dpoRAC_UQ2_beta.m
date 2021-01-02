% n # no. of Bob's input or # of Alice's dits
% d dimension of Bob's output.
% mode = 0 for preparaiton equivalence, mode = 1 for measurement equivalence
% OC = 0 for Spekkens' parity conditions, OC = 1 for even parity oblivious
% conditions
function [Pwin,Pcell] = nto1dpoRAC_UQ2_beta(n,d,mode,OC)

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
F5 = []; % positivity of probabilities.
%-------------------------------------------------------------------------------
for a = 1:AS
  Gamma{a}=sdpvar( 4*n*n*(d-1)*(d-1) + 2*n*(d-1)+1, 4*n*n*(d-1)*(d-1) + 2*n*(d-1)+1,'hermitian','real');
  F0=[F0;Gamma{a}>=0];
end
%-------------------------------------------------------------------------------
idx = @(b, B, U) 2*(d-1) * b + 2*B + U + 2 ; % function to return the position of the operator corresponding to Bob's setting (b,B) in a Gamma matrix.
idx22 = @(b, B, U, b1, B1, U1) 2*n*(d-1)+1 + 2*2*n*(d-1)*(d-1)*b + 2*2*n*(d-1)*B+ 2*n*(d-1)*U + 2*(d-1) * b1 +  2*B1 + U1 + 1 
%-------------------------------------------------------------------------------
% General Matrix constraints F1,F2,F3
for a = 1:AS
    for i = 1:4*n*n*(d-1)*(d-1) + 2*n*(d-1)+1
        F1 = [F1;Gamma{a}(i,i)==1];
  end
end
for a = 1:AS
  for b = 0:(n-1)
    for B = 0:(d-2)
        for U = 0:1
            for b1 = 0:(n-1)
                for B1 = 0:(d-2)
                    for U1 = 0:1
                        
            for br = 0:(n-1)
                for Br = 0:(d-2)
                    for Ur = 0:1
                        if (b == b1) && (B == B1) && (U~=U1)
                            F2=[F2;Gamma{a}(idx(br,Br,Ur),idx22(b,B,U,b1,B1,U1))==Gamma{a}(idx(br,Br,Ur),1)];
                            F2=[F2;Gamma{a}(idx22(b,B,U,b1,B1,U1),idx(br,Br,Ur))==Gamma{a}(1,idx(br,Br,Ur))];
                        elseif (b == br) && (B==Br) && (U==Ur)
                          % F2=[F2;Gamma{a}(idx(br,Br,Ur),idx22(b,B,U,b1,B1,U1))==Gamma{a}(1,idx(b1,B1,U1))];
                          % F2=[F2;Gamma{a}(idx22(b,B,U,b1,B1,U1),idx(br,Br,Ur))==Gamma{a}(idx(b1,B1,U1),1)];
                        elseif (b1 == br) && (B1==Br) && (U1==Ur)
                           F2=[F2;Gamma{a}(idx(br,Br,Ur),idx22(b,B,U,b1,B1,U1))==Gamma{a}(1,idx(b,B,U))];
                           F2=[F2;Gamma{a}(idx22(b,B,U,b1,B1,U1),idx(br,Br,Ur))==Gamma{a}(idx(b,B,U),1)];
                        end
                        for br1 = 0:(n-1)
                for Br1 = 0:(d-2)
                    for Ur1 = 0:1
                        if (br==b) && (Br == B) && (U==Ur)
                            F2=[F2;Gamma{a}(idx22(br,Br,Ur,br1,Br1,Ur1),idx22(b,B,U,b1,B1,U1))==Gamma{a}(idx(br1,Br1,Ur1),idx(b1,B1,U1))];
                            F2=[F2;Gamma{a}(idx22(b,B,U,b1,B1,U1),idx22(br,Br,Ur,br1,Br1,Ur1))==Gamma{a}(idx(b1,B1,U1),idx(br1,Br1,Ur1))];
                        elseif (br1==b1) && (Br1 == B1) && (Ur1==U1)
                            F2=[F2;Gamma{a}(idx22(br,Br,Ur,br1,Br1,Ur1),idx22(b,B,U,b1,B1,U1))==Gamma{a}(idx(br,Br,Ur),idx(b,B,U))];
                            F2=[F2;Gamma{a}(idx22(b,B,U,b1,B1,U1),idx22(br,Br,Ur,br1,Br1,Ur1))==Gamma{a}(idx(b,B,U),idx(br,Br,Ur))];
                        end
                    end
                end
                        end
                end
                end
            end
            
            
                    end
                end
            end
        end
    end
  end
end
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
if (mode == 1)
    for a = 1:AS
         for j = 1:4*n*n*(d-1)*(d-1) + 2*n*(d-1)+1
                   for B1 = 0:(d-3)
                    sum1=0; sum2=0; sum3=0;sum4=0;
                    for b1 = 0:(n-1)
                        sum1 = sum1 +  Gamma{a}(j, idx(b1, B1,0))+Gamma{a}(j, idx(b1, B1,1));
                        sum2 = sum2 + Gamma{a}(j, idx(b1, B1+1,0))+Gamma{a}(j, idx(b1, B1+1,1));
                        sum3 = sum3 + Gamma{a}(idx(b1, B1,0),j)+Gamma{a}(idx(b1, B1,1),j);
                        sum4 = sum4 + Gamma{a}(idx(b1, B1+1,0),j)+Gamma{a}(idx(b1, B1+1,1),j);
                    end
                    F4 = [F4; sum1-sum2==0; d*sum1 == 2*n*(2-d)*Gamma{a}(j,1);d*sum3 == 2*n*(2-d)*Gamma{a}(1,j);sum3==sum4];
                end
                end
    end
end
%-------------------------------------------------------------------------------
% Making a probability cell.
Pcell = cell(AS,n,d); % a0, a1, b, B
for a = 1:AS
    for b = 0:(n-1)
        sum=0;
        for B = 0:(d-2)
            Pcell{a,b+1,B+1 } = 0.5+(Gamma{a}(1, idx(b,B,0))+Gamma{a}(1, idx(b,B,1)))/4;
            F5 = [F5 ; Pcell{a,b+1,B+1 }>=0];
            sum=sum+Pcell{a,b+1,B+1 };
        end
         Pcell{a,b+1,d }=1-sum;
         F5 = [F5 ; Pcell{a,b+1,d }>=0];
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
Pwin = real(Pwin)/(AS*n)
diagnostics = optimize([F0;F1;F2;F3;F4;F5], -Pwin, sdpsettings('solver', 'mosek'))
Pwin = value(Pwin);
%--------------------------------------------------------------------------------
