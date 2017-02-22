function [varargout] = TypedLaplacian(W,k,epsilon,LazyRate)
% Spectral embedding of a typed directed network
%
% n is the number of nodes, and c is the number of different edge
% types
% W: an n*n*c weighted adjacency matrix. It can be
% undirected or directed, or a random walk matrix.
% epsilon: a value between 0 and 1, used to avoid the problem
% of reducibility of directed graphs -- the Google trick (default is 0)
% k: is the number of vectors desired, corresponding to the k smallest
% eigenvalues (default is all)
% LazyRate: a value between 0 and 1, the transition probability to another
% layer (default is 0.5)
%
% Vector = the eigenvector matrix
% e = the vector of Laplacian eigenvalues
% R = the large constructed random walk matrix created
%
% varargout = cell array
% 1: Vector=Eigenvectors
% 2: e= Eigenvalues
% 3: R= random walk matrix of the big connected graph
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
    epsilon=0;
else
    if epsilon>=1 || epsilon<0
        error('The range of the third argument (epsilon) is [0,1)');
    end
end;

if (nargin < 4)
    LazyRate=0.5;
end;

%%%% step 1: convert each layer to a random walk matrix
RW=zeros(n,n,c);
for i=1:c
    for j=1:n
        rs=sum(W(j,:,i));
        if rs==0;
            RW(j,j,i)=1;
        else
            RW(j,:,i)=W(j,:,i)/rs;
        end
    end
end

%%%%%Step 1-2: add epsilon
if epsilon>0
%    RW=(1-epsilon)*RW+epsilon/n;
  for i=1:c    
      RW(:,:,i)=(1-epsilon)*RW(:,:,i)+epsilon/(n-1)*(ones(n,n)-eye(n));
  end
end

%%%%step 2: convert whole graph to a random walk matrix R
R=zeros(c*n,c*n);
for i=1:c
    for j=1:c
        if i==j
            R((i-1)*n+1:i*n,(j-1)*n+1:j*n)= (1-LazyRate)*RW(:,:,i);
        else
            R((i-1)*n+1:i*n,(j-1)*n+1:j*n)= LazyRate/(c-1)*eye(n,n);
        end
    end
end

%%%%step 3: compute Fanchung's Directed Laplacian
    [ Cvector, e] = DLaplacian_Fan(R,0,k);
 Vector=zeros(n,k,c);
 for i=1:c
     Vector(:,:,i)=Cvector((i-1)*n+1:i*n,:);
 end

varargout{1} = Vector; varargout{2} = e; varargout{3} = R;
