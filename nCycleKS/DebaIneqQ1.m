function [Pwin,LB] = DebaIneqQ1(n)
AS = 1+3*n;
d=3;
%-------------------------------------------------------------------------------
% contraint vector declaration.
F0 = []; % contraint for positivity of Gamma. 
F1 = []; % constraints for orthogonality of results for the same measurement.
F2 = []; % constraints for summation over results of one measurement.
F3 = []; % constraints for summation over results of all measurements, and normalization.
F4 = []; % oblivious!
F5 = [];
%-------------------------------------------------------------------------------
% sdpvar declaration.
Gamma = cell(AS,1);
for a = 1:AS
  Gamma{a}=sdpvar(n*d,n*d,'hermitian','complex');
  F0=[F0;Gamma{a}>=0];
end
%-------------------------------------------------------------------------------
idx = @(b, B) d * b + B + 1 ; % function to return the position of the operator corresponding to Bob's setting (b,B) in a Gamma matrix.
idxA = @(i,k) 3 * (i-1) + k;
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
    end
  end
end
%-------------------------------------------------------------------------------     
for i = 1:n-1
    sum = 0; sum1=0;
    for k = 1:3
        sum = sum + Gamma{idxA(i,k)}/3;
        sum1 = sum1 + Gamma{idxA(i+1,k)}/3;
    end
    F4 = [F4; sum == sum1];
end

    
%-------------------------------------------------------------------------------     
for a = 1:AS
    for b = 0:(n-1)
        for B = 0:(d-1)
            for b1 = 0:(n-2)
                F5 = [F5;Gamma{a}(idx(b, B), idx(b1, 2))==Gamma{a}(idx(b, B), idx(b1+1, 0))];
                F5 = [F5;Gamma{a}(idx(b1, 2), idx(b, B))==Gamma{a}(idx(b1+1,0), idx(b, B))];
            end
            F5 = [F5;Gamma{a}(idx(b, B), idx(n-1, 2))==Gamma{a}(idx(b, B), idx(0, 0))];
                F5 = [F5;Gamma{a}(idx(n-1, 2), idx(b, B))==Gamma{a}(idx(0,0), idx(b, B))];
         end
    end
end
%-------------------------------------------------------------------------------     
% Making a probability cell.
Pcell = cell(AS,n,d); % a0, a1, b, B
for a = 1:AS
  for b = 0:(n-1)
    for B = 0:(d-1)
      Pcell{a,b+1,B+1 } = Gamma{a}(idx(b,B), idx(b,B));
      F4 = [F4 ; Gamma{a}(idx(b,B), idx(b,B))>=0];
    end
  end
end
%-------------------------------------------------------------------------------  
sum = 0;
for i = 1:n
    for k = 1:3
        sum = sum + Pcell{idxA(i,k),i,k}/n
    end
    sum = sum + Pcell{1+3*n,i,1}
end
Pwin = real(sum)
%-------------------------------------------------------------------------------  
Pwin = real(sum);
diagnostics = optimize([F0;F1;F2;F3;F4;F5], -Pwin, sdpsettings('solver', 'sdpt3','verbose',0));
Pwin = value(Pwin);
if rem(n,2)==0
   LB = 3+n/2;
else
   LB = 3 + (n*cos(pi/n))/(1+cos(pi/n));
end
