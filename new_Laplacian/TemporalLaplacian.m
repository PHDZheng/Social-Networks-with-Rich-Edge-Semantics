function [varargout] = TemporalLaplacian(W,k,alpha,beta,epsilon)
% Spectral embedding of a directed graph with changing weights
% using Fan Chung¡¯s directed embedding technique
%
% n is the number of nodes, and c is the number of different time
% periods
%
% W: an n*n*c weighted adjacency matrix of the graph. It can be
% undirected or directed, or a random walk matrix.
% epsilon: is a value between 0 and 1, used to avoid the problem
% of reducibility of directed graph -- the Google trick (default is 0)
% k: the number of vectors desired, corresponding to the k smallest
% eigenvalues (default is all)
% alpha: a value between 0 and 1, the down-weighting value for
% combining matrices from previous time periods
% beta: a value between 0 and 1, the probability of transitioning out
% of a layer (default is 0.5)
%
% Vector = the matrix of eigenvectors
% e = the vector of Laplacian eigenvalues
% aggregateRW: = the large constructed random walk matrix R
%
% varargout = cell array
% 1: Vector=Eigenvectors
% 2: e= Eigenvalues
% 3: aggregateRW=(cn*cn) matrix
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
if (nargin < 5)
    epsilon=0;
else
    if epsilon>=1 || epsilon<0
        error('The range of the third argument (epsilon) is [0,1)');
    end
end;



aggregateRW=sparse(c*n,c*n);
for i=1:c
    if i==1 %%%% convert to the aggregate whole graph with alpha forward effection
        A=W(:,:,1);
    else
        A=alpha*A+W(:,:,i);
    end
    temprw=sparse(n,n); %% convert aggregate graph to a random walk matrix
    for j=1:n
        rs=sum(A(j,:));
        if rs==0;
            temprw(j,j)=1;
        else
            temprw(j,:)=A(j,:)/rs;
        end
    end
    if epsilon>0
          temprw=(1-epsilon)*temprw+epsilon/(n-1)*(ones(n,n)-eye(n));
    end
    for j=1:c %% minimize the disagreement between layer with beta trade in value
        if i==j
            aggregateRW((i-1)*n+1:i*n,(j-1)*n+1:j*n)= (1-beta)*temprw;
        else
            aggregateRW((i-1)*n+1:i*n,(j-1)*n+1:j*n)= beta/(c-1)*eye(n,n);
        end
    end
    end

%%%% compute Fanchung's Directed Laplacian
    [ Cvector, e] = DLaplacian_Fan(aggregateRW,0,k);
 Vector=zeros(n,k,c);
 for i=1:c
     Vector(:,:,i)=Cvector((i-1)*n+1:i*n,:);
 end

varargout{1} = Vector; varargout{2} = e; varargout{3} = aggregateRW;
