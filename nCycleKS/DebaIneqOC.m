function [Pwin,LB] = DebaIneqOC(n)
AS = 1+3*n;
d=3;
%-------------------------------------------------------------------------------    
 % function to return the position corresponding to the probability P(E,B|e,b,AS) in the probability vector.
idx = @(a,b,B) (n*d) * (a-1) + d*(b-1) + B
idxA = @(i,k) 3 * (i-1) + k;
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
        for i = 1:n-1
            sum = 0; sum1=0;
            for k = 1:3
                sum = sum + x(idx(idxA(i,k),b,B))/3;
                sum1 = sum1 + x(idx(idxA(i+1,k),b,B))/3;
            end
            con = [con; sum == sum1];
        end
    end
end
    
%-------------------------------------------------------------------------------     
for a = 1:AS
    for b1 = 1:(n-1)
        con = [con;x(idx(a,b1, 3))==x(idx(a,b1+1, 1))];
    end
        con = [con;x(idx(a,n, 3))==x(idx(a,1, 3))];
end
%-------------------------------------------------------------------------------     
% Making a probability cell.
Pcell = cell(AS,n,d); % a0, a1, b, B
for a = 1:AS
  for b = 0:(n-1)
    for B = 0:(d-1)
      Pcell{a,b+1,B+1 } = x(idx(a,b+1,B+1));
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
diagnostics = optimize([con], -Pwin, sdpsettings('solver', 'mosek','verbose',0));
Pwin = value(Pwin);
if rem(n,2)==0
   LB = 3+n/2;
else
   LB = 3 + (n*cos(pi/n))/(1+cos(pi/n));
end
