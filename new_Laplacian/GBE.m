%posW: is the n*n weighted adjacency matrix of a undirected graph.
%label: is a n*c label indication matrix, with value {-1, 0 1} if label is 
%       a vector,{1,0} if label is a matrix; 0s rows for unlabelled nodes.
% ew: is a [apw,anw] vector for positive and negative added edge weight 
%                         between labeled points; default value is [1,1]. 
function F = GBE(posW,label,ew)
[n,m]=size(label);
if m<=2 % two classes
    if m==1
        Y(:,1)=label;
        Y(:,2)=-label;
        Y(Y<0)=0;
    else
        Y=label;
    end
    F = sign(BCGE(posW,Y,ew));
    
else % more than two classes
    for i=1:m
        Y=[label(:,i),sum(label,2)-label(:,i)];
        V(:,i)=BCGE(posW,Y,ew);
    end
    [C,F] = max(V,[],2);
end