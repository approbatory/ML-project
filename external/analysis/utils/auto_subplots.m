function [p,n]=auto_subplots(n)
% function [p,n]=auto_subplots(n)
%
% Calculate how many rows and columns of sub-plots are needed to
% neatly display n subplots. 
%
% Inputs
% n - the desired number of subplots.     
%  
% Outputs
% p - a vector length 2 defining the number of rows and number of
%     columns required to show n plots.     
% [ n - the current number of subplots. This output is used only by
%       this function for a recursive call.]
%
% Rob Campbell - January 2010
% https://www.mathworks.com/matlabcentral/fileexchange/26310-numsubplots-neatly-arrange-subplots

while isprime(n) & n>4, 
    n=n+1;
end

p=factor(n);

if length(p)==1
    p=[1,p];
    return
end


while length(p)>2
    if length(p)>=4
        p(1)=p(1)*p(end-1);
        p(2)=p(2)*p(end);
        p(end-1:end)=[];
    else
        p(1)=p(1)*p(2);
        p(2)=[];
    end    
    p=sort(p);
end


%Reformat if the column/row ratio is too large: we want a roughly
%square design 
while p(2)/p(1)>2.5
    N=n+1;
    [p,n]=numSubplots(N); %Recursive!
end
