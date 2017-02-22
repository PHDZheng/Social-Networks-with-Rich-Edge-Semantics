function [varargout] = TypedDirLaplacian(W,k,LazyRate)
% Spectral embedding of a graph with directed typed edges
%
% n is the number of nodes, and c is the number of different edge
% types.
% W: an n*n*c weighted adjacency matrix of the graph. It can be
% undirected or directed, or a random walk matrix.
% LazyRate: a value between 0 and 1, the probability of moving
% to another layer (default is 0.5)
% k: the number of vectors desired, corresponding to the k smallest
% eigenvalues (default is all)
%
% Vout = the eigenvector matrix of out-roles
% Vin = the eigenvector matrix of in-roles
% e = vector of Laplacian eigenvalues
% bigW = the large constructed matrix
%
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
    LazyRate=0.5;
end;

%%%% step 1: Bind together the versions of each node in the different
%%%% layers to build a cn*cn adjacency matrix
bigW=sparse(c*n,c*n);
for i=1:c
    rowsumd=sum(W(:,:,i),2);
    tempv=sparse(n,1);
    tempv(rowsumd==0)=1;
    W(:,:,i)=W(:,:,i)+diag(tempv);
    rowsumd(rowsumd==0)=1;
    tempD=diag(rowsumd); 
    for j=1:c
        if i==j
            bigW((i-1)*n+1:i*n,(j-1)*n+1:j*n)= (1-LazyRate)*W(:,:,i);
        else
            bigW((i-1)*n+1:i*n,(j-1)*n+1:j*n)= LazyRate/(c-1)*tempD;
        end
    end
end

%%%%step 3: compute our new's Directed Laplacian
[DirOut,DirIn,e] = DirLaplacian(bigW,k);
 Vout=zeros(n,k,c);
  Vin=zeros(n,k,c);
 for i=1:c
     Vout(:,:,i)=DirOut((i-1)*n+1:i*n,:);
     Vin(:,:,i)=DirIn((i-1)*n+1:i*n,:);
 end

varargout{1} = Vout; varargout{2} = Vin; varargout{3} = e; varargout{4} = bigW;
