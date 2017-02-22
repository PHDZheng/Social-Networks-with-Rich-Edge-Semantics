%posW: is the n*n weighted adjacency matrix of a undirected graph.
%Y: is a n*2 label indication matrix,  with value { 0 1};
%                                   0s rows for unlabelled nodes.
% ew: is a [apw,anw] vector for positive and negative added edge weight
%                         between labeled points; default value is [1,1].                      ;

function F = BCGE(W,Y,ew)
if nargin < 2
    error('At least two input arguments required.');
end;
if (nargin < 3)
    ew=[1,1];
end;
if W ~= W'
    error('W is not symetric.');
end;
[n,m]=size(Y);
if size(W,2)~=n;
    error('Label is not consistent with the matrix.');
end
if m~=2;
    error('Y must be n*2 matrix.');
end

% build the negatvie matrix with only one undirected edge.
negW=sparse(n+2,n+2);
negW(n+1,n+2)=ew(2);
negW(n+2,n+1)=ew(2);

% build the positve matrix with adjacency matrix and added positive edges.
posW=sparse(n+2,n+2);
posW(1:n,1:n)=W;
posW(1:n,n+1)=ew(1)*Y(:,1);
posW(n+1,1:n)=ew(1)*Y(:,1)';
posW(1:n,n+2)=ew(1)*Y(:,2);
posW(n+2,1:n)=ew(1)*Y(:,2)';

% build the modified total degree matrixes.
Dpos=diag(sum(posW,2));
Dneg=diag(sum(negW,2));
D=diag(sum(W,2));
D(n+1,n+1)=Dpos(n+1,n+1)+Dneg(n+1,n+1);
D(n+2,n+2)=Dpos(n+2,n+2)+Dneg(n+2,n+2);
D2=sparse(n+2,n+2);
for i=1:n+2
    if D(i,i)~=0
        D2(i,i)=1/D(i,i)^.5;
    end
end

% Signed laplacian embedding
L=(Dpos-posW)-(Dneg-negW);
Lsym=max(sum(abs(D2*L*D2),2))*speye(n+2)-D2*L*D2; %- max(sum(abs(D2*L*D2),2))*D3*D3';
Lsym=Lsym+Lsym';
[V,e] = eigs(Lsym,2);

if abs(e(1,1)-2*max(sum(abs(D2*L*D2),2)))<1e-10;
    U=V(:,2);
else
    U=V(:,1);
end
% make the output vector consistent with the labels.
if U(n+1,1)>0
    F=U(1:n,1);
else
    F=-U(1:n,1);
end

