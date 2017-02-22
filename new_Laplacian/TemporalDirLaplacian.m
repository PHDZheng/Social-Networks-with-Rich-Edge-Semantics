function [varargout] = TemporalLaplacian(W,k,alpha,beta)
% Spectral embedding of a directed graph changing with time
%
% n is the number of nodes, and c is the number of different time
% periods
%
% W: is an n*n*c weighted adjacency matrix of the graph. It can be
% undirected or directed, or a random walk matrix.
% k: the number of eigenvectors corresponding to the k smallest eigenvalues
% desired (default is all);
% alpha: a value between 0 and 1, down-weighting the contribution of
% matrices from previous time periods
% beta: a value between 0 and 1, the total probability of a random
% walk transitioning to one of the other layers in the large constructed
% graph (default is 0.5)
%
% Vout = the eigenvector matrix of out-roles
% Vin = the eigenvector matrix of in-roles
% e = the vector of Laplacian eigenvalue
% bigW = the large constructed adjacency matrix
% for rendering purposes
% varargout = cell array
% 1: Vout
% 2: Vin
% 3: e
% 4: bigW
if nargin < 1
    error('At least one input arguments required.');
end;
[n,m,c]=size(W);
if n~=m || c==1
    error('The first argument have to be a n*n*c matrix. (c>=2)');
end;
if (nargin < 2)
    k=c*n;
end;
if (nargin < 3)
    alpha=0.5;
else
    if alpha>=1 || alpha<0
        error('The range of the third argument (alpha) is [0,1)');
    end
end;

if (nargin < 4)
    beta=0.5;
    else
    if beta>=1 || beta<0
        error('The range of the fourth argument (beta) is [0,1)');
    end
end;

newD=zeros(n,c);

%%%% convert to the aggregate whole graph with alpha forward effection
%A=sparse(n,n,c);
 A(:,:,1)=W(:,:,1); 
 for i=2:c
     A(:,:,i)=alpha*A(:,:,i-1)+W(:,:,i);
 end
 %%%% compute our typed Directed Laplacian
    [Vout,Vin,e,bigW]=TypedDirLaplacian(A,k,beta);
varargout{1} = Vout; varargout{2} = Vin; varargout{3} = e; varargout{4} = bigW;
