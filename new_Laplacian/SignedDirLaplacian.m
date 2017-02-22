function [varargout] = SignedDirLaplacian(posW,negW,type,k)
% Spectral embedding of a signed directed network
%
% W an n*n*c weighted adjacency matrix
% with n nodes, and c different edge types
% k is the number of eigenvectors desired, corresponding to the k
% smallest eigenvalues (default is all);
%
% varargout = cell array of output
% 1: PosOut=Eigenvectors of positive Out role
% 2: Negout=Eigenvectors of negative Out role
% 3: PosIn=Eigenvectors of positive in role
% 4: NegIn=Eigenvectors of negative in role
% 5: e=Eigenvalues
% 6: Xpos= the positive adjacency matrix created
% by modelling the signed directed graph
% 6: Xneg= the negative adjacency matrix created
% by modelling the signed directed graph
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

%%%% step 1: Bind together the versions of each node in the 4 different
%%%% copies
Dinpos = sum(posW,1);
Doutpos = sum(posW,2);
Dinneg = sum(negW,1);
Doutneg = sum(negW,2);
crossweights = Dinpos + Doutpos' + Dinneg + Doutneg';
D = diag(crossweights);
Wp = [posW + D, D; D, D];

bigD = sparse(2*n);
bigD(1:n, n+1:2*n)=D;
bigD(n+1:2*n, 1:n)=D;
Xpos = [ bigD, Wp; Wp', bigD]; 
Xneg= sparse(4*n);
Xneg(n+1:2*n, 3*n+1:4*n)=negW;
Xneg(3*n+1:4*n, n+1:2*n)=negW';

%%%%step 2: compute our signed Laplacian
[Vector,e] = SignedLaplacian(Xpos,Xneg,type,k);
PosOut=Vector(1:n,:);
Negout=Vector(n+1:2*n,:);
PosIn=Vector(2*n+1:3*n,:);
NegIn=Vector(3*n+1:4*n,:);
if nargout==1
    varargout{1} = Vector;
elseif nargout==2
    varargout{1} = Vector; varargout{2} = e;
else
    varargout{1} = PosOut; varargout{2} = Negout; varargout{3} = PosIn; varargout{4} = NegIn;  varargout{5} = e; varargout{6} = Xpos; varargout{7} = Xneg;
end



