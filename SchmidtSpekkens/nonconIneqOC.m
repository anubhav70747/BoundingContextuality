function [Pwin,Pcellv] = nonconIneqOC()
AS = 6;
n = 3;
d=2;
%-------------------------------------------------------------------------------    
 % function to return the position corresponding to the probability P(E,B|e,b,AS) in the probability vector.
idx = @(a,b,B) (n*d) * (a-1) + d*(b-1) + B 
%-------------------------------------------------------------------------------      
x=sdpvar( AS * n * d , 1); % the probability vector as the variable for lin prog.
%-------------------------------------------------------------------------------  
con = [0<=x<=1];  %positivity constraint
%-------------------------------------------------------------------------------  
for a = 1:AS
    for b = 1:n
        sum = 0;
        for B = 1:d
            sum = sum + x(idx(a,b,B));
        end
        con = [con; sum==1];
    end
end
%-------------------------------------------------------------------------------  
for b = 1:n
    for B = 1:d
        con = [con; x(idx(1,b,B))+x(idx(2,b,B))==x(idx(3,b,B))+x(idx(4,b,B));x(idx(1,b,B))+x(idx(2,b,B))==x(idx(5,b,B))+x(idx(6,b,B)) ];
    end
end
%-------------------------------------------------------------------------------  
for a = 1:AS
      con = [con; x(idx(a,1,1))+x(idx(a,2,1))+x(idx(a,3,1))==x(idx(a,1,2))+x(idx(a,2,2))+x(idx(a,3,2))];
end
%-------------------------------------------------------------------------------  
Pcell = cell(AS,n,d);
for a = 1:AS
    for b = 1:n
        for B = 1:d
            Pcell{a,b,B} = x(idx(a,b,B));
        end
    end
end
%-------------------------------------------------------------------------------  
Pwin = [];
Pwin1 = Pcell{1,1,1} + Pcell{3,2,1} + Pcell{5,3,1} ;
diagnostics = optimize(con,-Pwin1,sdpsettings('solver', 'sdpt3'))
Pwin = [Pwin;value(Pwin1)];
%-------------------------------------------------------------------------------  
Pwin2 = Pcell{1,1,1} + Pcell{2,2,1} + Pcell{5,3,1} ;
diagnostics = optimize(con,-Pwin2,sdpsettings('solver', 'sdpt3'))
Pwin = [Pwin;value(Pwin2)];
%-------------------------------------------------------------------------------  
Pwin3 = Pcell{1,1,1} - Pcell{3,1,1} -2 * Pcell{5,1,1} -2 * Pcell{2,2,1} + 2 * Pcell{3,2,1} + 2 * Pcell{5,3,1}   ;
diagnostics = optimize(con,-Pwin3,sdpsettings('solver', 'sdpt3'))
Pwin = [Pwin;value(Pwin3)];
%-------------------------------------------------------------------------------  
Pwin4 = 2* Pcell{1,1,1} - Pcell{2,2,1} +2* Pcell{3,2,1}; 
diagnostics = optimize(con,-Pwin4,sdpsettings('solver', 'sdpt3'))
Pwin = [Pwin;value(Pwin4)];
%-------------------------------------------------------------------------------  
Pwin5 = Pcell{1,1,1} - Pcell{5,1,1} +  Pcell{2,2,1} + Pcell{3,2,1} + 2 * Pcell{5,3,1};
diagnostics = optimize(con,-Pwin5,sdpsettings('solver', 'sdpt3'))
Pwin = [Pwin;value(Pwin5)];
%-------------------------------------------------------------------------------  
Pwin6 = Pcell{1,1,1} - Pcell{5,1,1} +  2*Pcell{2,2,1} + 2 * Pcell{5,3,1};
diagnostics = optimize(con,-Pwin6,sdpsettings('solver', 'sdpt3'))
Pwin = [Pwin;value(Pwin6)];
%-------------------------------------------------------------------------------
Pwin7 = real(Pcell{1,1,1} - Pcell{4,1,1} -2 * Pcell{5,1,1} -2 * Pcell{2,2,1} + 2 * Pcell{3,2,1} + 2 * Pcell{5,3,1});
diagnostics = optimize([con], -Pwin3, sdpsettings('solver', 'sdpt3'))
Pwin = [Pwin;value(Pwin3)];
%-------------------------------------------------------------------------------  


