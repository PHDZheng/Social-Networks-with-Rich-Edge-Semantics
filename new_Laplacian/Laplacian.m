function [varargout] = Laplacian(W,type,k)
% Spectral embedding of an undirected graph
%
% W is the weighted adjacency matrix of a graph.
% type = un and empty: unnormalized graph Laplacian
%        sym: symmetric normalized graph Laplacian
%        rw: random walk normalized graph Laplacian
% k = number of groups to cluster into (default is disabled)
%
% Vector = eigenvector matrix of the embedding
% e = the vector of Laplacian eigenvalues
%
% varargout = cellarray
% 1: Vector
% 2: e

if nargin < 1
    error('At least one input arguments required.');
end
[n,m]=size(W);
if n~=m 
    error('first argument has to be a n*n matrix.');
else
    tempvalue=max(max(abs(W-W')));
    if tempvalue>1e-5
        error('The matrix must be undirected.');
    end
end;
if (nargin < 2)
    type='RW';
end;
if (nargin < 3)
    k=n;
else
    if k>n
       error('Third argument k, the number of smallest eigenvalues and eigenvectors, has to be less than or equal to the dimensions of the adjacency matrix.'); 
    end;
end;

    tempvalue=max(max(W-W'));
    if tempvalue>1e-10
        error('The adjacency matrix for unnormalized Laplacian decomposition need to be symmetric.');
    end
    if  tempvalue~=0
        W=(W+W')/2;
    end        
    D=sparse(1:n,1:n,sum(W,2));
    L=D-W;
    
%unnormalized graph Laplacian
if (nargin < 2) || isempty(type) || strcmp(type,'un') || strcmp(type,'Un') || strcmp(type,'UN')
    display('Unnormalized Laplacian decomposition');
    tempvalue = max(D(:));
    if k == n
        [Vector,eigenvalue] = eig(full(L));
        e = diag(eigenvalue);
        [e, IXY] = sort(e,1,'ascend');
        Vector = Vector(:,IXY);
    else
        [Vector,eigenvalue] = eigs(2 * tempvalue * eye(n)-L,k);
        e = 2 * tempvalue - diag(eigenvalue);
    end
    varargout{1} = Vector;
    varargout{2} = e;
    % % % symmetric normalized graph-Laplacian
elseif strcmp(type,'sym') || strcmp(type,'Sym') || strcmp(type,'SYM')
    % Normalized spectral clustering according to Ng, Jordan, and Weiss
    % (2002)
    display('Symmetric Laplacian decomposition');
    D2 = sparse(n,n);
    for i=1:n
        if D(i,i) ~= 0
            D2(i,i)=1/D(i,i)^0.5;
        end
    end
    Lsym = D2 * L * D2; % L=D?(-0.5)*L*D?(-0.5);
    Lsym= (Lsym + Lsym')/2;
    if k == n
        [Vector,eigenvalue] = eig(full(Lsym));
        e = diag(eigenvalue);
        [e,IXY] = sort(e,1,'ascend');
        Vector = Vector(:,IXY);
    else
        [Vector,eigenvalue] = eigs(2*eye(n) - Lsym,k);
        e = 2 - diag(eigenvalue);
    end
    varargout{1} = Vector;
    varargout{2} = e;
    % random walk normalized graph-Laplacian
elseif strcmp(type,'rw') || strcmp(type,'Rw') || strcmp(type,'RW')
    display('Random walk Laplacian decomposition');
    D2 = sparse(n,n); % D2= D?(-0.5);
    for i=1:n
        if D(i,i) ~= 0
            D2(i,i) = 1/D(i,i)^.5;
        end
    end
    Lsym = D2 * L * D2; %L=D?(-0.5)*L*D?(-0.5);
    Lsym = (Lsym + Lsym')/2;
    if k == n
        [Vector,eigenvalue] = eig(full(Lsym));
        e = diag(eigenvalue);
    else
        [Vector,eigenvalue] = eigs(2*eye(n)-Lsym,k);
        e=2-diag(eigenvalue);
    end
    Vector = D2*Vector;
    varargout{1} = Vector;
    varargout{2} = e;
else
    error('Type cannot be identified.');
end





