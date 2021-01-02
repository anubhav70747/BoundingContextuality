function Pwin = nonconIneqPiQ1()
AS = 6;
n = 3;
d=2;
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
F4 = [F4;Gamma{1}+Gamma{2}==Gamma{3}+Gamma{4};Gamma{1}+Gamma{2}==Gamma{5}+Gamma{6}];
%-------------------------------------------------------------------------------     
for a = 1:AS
    for b = 0:(n-1)
        for B = 0:(d-1)
            F5 = [F5;Gamma{a}(idx(b, B), idx(0, 0))+Gamma{a}(idx(b, B), idx(1, 0))+Gamma{a}(idx(b, B), idx(2,0))==Gamma{a}(idx(b, B), idx(0, 1))+Gamma{a}(idx(b, B), idx(1, 1))+Gamma{a}(idx(b, B), idx(2,1))];
            F5 = [F5;Gamma{a}(idx(0, 0), idx(b, B))+Gamma{a}(idx(1, 0), idx(b, B))+Gamma{a}(idx(2,0), idx(b, B))==Gamma{a}(idx(0, 1), idx(b, B))+Gamma{a}(idx(1, 1), idx(b, B))+Gamma{a}(idx(2,1), idx(b, B))];
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
Pwin = [];
%-------------------------------------------------------------------------------  
Pwin1 = real(Pcell{1,1,1} + Pcell{3,2,1} + Pcell{5,3,1}) ;
diagnostics = optimize([F0;F1;F2;F3;F4;F5;], -Pwin1, sdpsettings('solver', 'sdpt3'))
Pwin = [Pwin;value(Pwin1)];
%-------------------------------------------------------------------------------  
Pwin2 = real(Pcell{1,1,1} + Pcell{2,2,1} + Pcell{5,3,1});
diagnostics = optimize([F0;F1;F2;F3;F4;F5;], -Pwin2, sdpsettings('solver', 'sdpt3'))
Pwin = [Pwin;value(Pwin2)];
%------------------------------------------------------------------------------- 
Pwin3 = real(Pcell{1,1,1} - Pcell{3,1,1} -2 * Pcell{5,1,1} -2 * Pcell{2,2,1} + 2 * Pcell{3,2,1} + 2 * Pcell{5,3,1});
diagnostics = optimize([F0;F1;F2;F3;F4;F5;], -Pwin3, sdpsettings('solver', 'sdpt3'))
Pwin = [Pwin;value(Pwin3)];
%------------------------------------------------------------------------------- 
Pwin4 = real(2* Pcell{1,1,1} - Pcell{2,2,1} +2* Pcell{3,2,1}); 
diagnostics = optimize([F0;F1;F2;F3;F4;F5;], -Pwin4, sdpsettings('solver', 'sdpt3'))
Pwin = [Pwin;value(Pwin4)];
%------------------------------------------------------------------------------- 
Pwin5 = real(Pcell{1,1,1} - Pcell{5,1,1} +  Pcell{2,2,1} + Pcell{3,2,1} + 2 * Pcell{5,3,1});
diagnostics = optimize([F0;F1;F2;F3;F4;F5;], -Pwin5, sdpsettings('solver', 'sdpt3'))
Pwin = [Pwin;value(Pwin5)];
%------------------------------------------------------------------------------- 
Pwin6 = real(Pcell{1,1,1} - Pcell{5,1,1} +  2*Pcell{2,2,1} + 2 * Pcell{5,3,1});
diagnostics = optimize([F0;F1;F2;F3;F4;F5;], -Pwin6, sdpsettings('solver', 'sdpt3'))
Pwin = [Pwin;value(Pwin6)];
%-------------------------------------------------------------------------------  
Pwin7 = real(Pcell{1,1,1} - Pcell{4,1,1} -2 * Pcell{5,1,1} -2 * Pcell{2,2,1} + 2 * Pcell{3,2,1} + 2 * Pcell{5,3,1});
diagnostics = optimize([F0;F1;F2;F3;F4;F5;], -Pwin7, sdpsettings('solver', 'sdpt3'))
Pwin = [Pwin;value(Pwin7)];
