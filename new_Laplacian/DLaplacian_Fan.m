function [varargout] = DLaplacian_Fan(W,epsilon,k)
% Spectral embedding of a directed graph using the method of Fan Chung
%
% W is the weighted adjacency matrix of a directed graph.
% epsilon: is a value between 0 and 1, used to avoid the problem
% of reducibility of directed graph (default is 0)
% k = number of eigenvectors desired
% Vector = the eigenvector matrix for the embedding
% e = the vector of Laplacian eigenvalues
% importance = importance value of the random walk matrix
% varargout = cell array of results
% 1: Vector
% 2: e
% 3: importance
if nargin < 1||nargin >3
    error('At least one input argument required, at most three input arguments required.');
end
if nargin >=2
    if epsilon>=1 || epsilon<0
        error('0 <= epsilon < 1.');
    end
else
    epsilon=0;
end
n=size(W,1);
D=diag(sum(W,2));

    RW = sparse(n,n); %RW is the probability matrix of w
    for i = 1:n
        if D(i,i) ~= 0
            RW(i,:) = W(i,:)./D(i,i);
        else
            RW(i,i) = 1;
        end
    end
    
if epsilon>0
      RW=(1-epsilon)*reshape(RW,n,n)+epsilon/(n-1)*(ones(n,n)-eye(n));
end
    % compute left eigenvector of RW -- the probability
    [pie, eigpie] = eigs(RW.',1);
     importance = abs(pie);

    %build symmetrical directed Laplacian
    imphalf = diag(sqrt(importance));
    impminushalf = sparse(n,n);
    for i = 1:n
        if importance(i) >1e-10;
            impminushalf(i,i) = 1./imphalf(i,i);
        end
    end
    L = speye(n,n) - (imphalf * RW * impminushalf + impminushalf * RW' * imphalf)/2;
    
if nargin<3
    [Vector,E]=eig(L);
    [e,IX] = sort(diag(E));
    Vector=Vector(:,IX);
else
    [Vector,E]=eigs(2*speye(n)-L,k);
    e=2-diag(E);
end
 Vector=impminushalf*Vector;
% for i=1:n
%     Vector(:,i)=Vector(:,i)./norm(Vector(:,i));
% end
varargout{1} = Vector; varargout{2} = e; varargout{3} = importance;


