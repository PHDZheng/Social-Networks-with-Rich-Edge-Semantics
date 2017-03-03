function [varargout] = SignedLaplacian(posW,negW,type,k)
% Spectral embedding of a signed adjacency matrix
%
% The matrices posW and negW are the weighted adjacency matrices of the
% positive and negative edges of the network
% type= un or empty: unnormalized graph Laplacian
%               SNS: simple normalized signed graph Laplacian
%               BNS: balanced normalized signed graph Laplacian
%
% k = cluster the graph into k groups (default is disabled)
%
% Vector = eigenvector matrix for embedding
% e = the vector of the Laplacian eigenvalues
%
% varargout = cell array
% 1: Vector
% 2: e
if nargin < 2
    error('At least two input arguments required.');
end
[n,m]=size(posW);
if n~=m 
    error('first argument has to be a n*n matrix.');
elseif n~=size(negW,1) || n~=size(negW,2)
    error('Second argument has to be a n*n matrix as the first.');
end;
if (nargin < 3)
    type='SNS';
end;
if (nargin < 4)
    k=n;
else
    if k>n
       error('Third argument k, the number of smallest eigenvalues and eigenvectors, has to be less than or equal to the dimensions of the adjacency matrix.'); 
    end;
end;

    tempvalue=max(max(posW-posW'));
    if tempvalue>1e-10
        error('The positive adjacency matrix for unnormalized Laplacian decomposition need to be symmetric.');
    end
    if  tempvalue~=0
        posW=(posW+posW')/2;
    end 
    
    tempvalue=max(max(negW-negW'));
    if tempvalue>1e-10
        error('The negative adjacency matrix for unnormalized Laplacian decomposition need to be symmetric.');
    end
    if  tempvalue~=0
        negW=(negW+negW')/2;
    end 
    
Dpos=sparse(1:n,1:n,sum(posW,2));
Dneg=sparse(1:n,1:n,sum(negW,2));
D=Dpos+Dneg;
D2=zeros(n,n);
for i=1:n
    if D(i,i)~=0
        D2(i,i)=1/D(i,i)^.5;
    end
end

%unnormalized graph-Laplacian
if (nargin < 3) || isempty(type)||strcmp(type,'un') || strcmp(type,'Un')||strcmp(type,'UN')
    display('Unnormalized signed Laplacian decomposition');
    L=Dpos-posW-Dneg+negW;
    tempvalue=max(Dpos(:));
    if k==n
        [Vector,eigenvalue]=eig(full(L));
        e=diag(eigenvalue);
        [e,IXY]=sort(e,1,'ascend');
        Vector=Vector(:,IXY);
    else
        [Vector,eigenvalue]=eigs(2*tempvalue*eye(n)-L,k);
        e=2*tempvalue-diag(eigenvalue);
    end
    varargout{1} = Vector; varargout{2} = e;

    
    %simple normalized signed graph Laplacian
elseif strcmp(type,'SNS') || strcmp(type,'Sns')||strcmp(type,'sns')
    display('Simple normalized signed graph Laplacian decomposition');
    L=Dpos-posW-Dneg+negW;
    D=Dpos+Dneg;
D2=zeros(n,n);
for i=1:n
    if D(i,i)~=0
        D2(i,i)=1/D(i,i)^.5;
    end
end
    Lsym=D2*L*D2;  %L=D^(-0.5)*L*D^(-0.5);
    Lsym=(Lsym+Lsym')/2;
    if k==n
        [Vector,eigenvalue]=eig(full(Lsym));
        e=diag(eigenvalue);
        [e,IXY]=sort(e,1,'ascend');
        Vector=Vector(:,IXY);
    else
        [Vector,eigenvalue]=eigs(2*eye(n)-Lsym,k);
        e=2-diag(eigenvalue);
    end
    Vector=D2*Vector;
    varargout{1} = Vector; varargout{2} = e;

    
    % balanced normalized signed graph Laplacian
elseif strcmp(type,'BNS') || strcmp(type,'Bns')||strcmp(type,'bns')
    display('Balanced normalized signed graph Laplacian decomposition');
    L=Dpos-posW+negW;
D=Dpos+Dneg;
D2=zeros(n,n);
for i=1:n
    if D(i,i)~=0
        D2(i,i)=1/D(i,i)^.5;
    end
end
    Lsym=D2*L*D2;  %L=D^(-0.5)*L*D^(-0.5);
    Lsym=(Lsym+Lsym')/2;
    if k==n
        [Vector,eigenvalue]=eig(full(Lsym));
        e=diag(eigenvalue);
        [e,IXY]=sort(e,1,'ascend');
        Vector=Vector(:,IXY);
    else
        [Vector,eigenvalue]=eigs(2*eye(n)-Lsym,k);
        e=2-diag(eigenvalue);
    end
    Vector=D2*Vector;
    varargout{1} = Vector; varargout{2} = e;
    
else
    error('Type can not be indentified.');
end


