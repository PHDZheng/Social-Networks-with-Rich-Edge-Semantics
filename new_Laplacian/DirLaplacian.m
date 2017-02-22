function [varargout] = DirLaplacian(W,k)
% Creates spectral embedding from a directed adjacency matrix
%
% The matrix W is the weighted adjacency matrix of a directed graph.
% k = number of eigenvectors.
% Vout = the eigenvector of out versions of nodes
% Vin = the eigenvector of in versions of nodes
% e = the vector of Laplacian eigenvalues
% varargout = cell array of results
% 1: Vout
% 2: Vin
% 3: e
if nargin < 1||nargin >3
    error('At least one input argument required, at most three input arguments required.');
end
if nargin < 1
    error('At least one input arguments required.');
end
[n,m]=size(W);
if (nargin < 2)
    k=n;
else
    if k>n
       error('Third argument k, the number of smallest eigenvalues and eigenvectors, has to be less than or equal to the dimensions of the adjacency matrix.'); 
    end;
end;

if n~=m 
    error('first argument has to be a n*n matrix.');
end;
din=sum(W,1)';
dout=sum(W,2);
X=sparse(diag((2*dout+din).^(-0.5)))*(W+diag(din+dout))*sparse(diag((dout+2*din).^(-0.5)));
if k == n
   [ U, E, V] = svd(X);
else
   [ U, E, V] = svds(X,k);
end   
    e=1-diag(E);
Vout=diag((2*dout+din).^(-0.5))*U;
Vin =diag((dout+2*din).^(-0.5))*V;

varargout{1} = Vout; varargout{2} = Vin; varargout{3} = e;


