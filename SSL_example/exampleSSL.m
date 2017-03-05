clc
clear all
close all

path(path, '..\new_Laplacian'); % the path of Library

%%% created a four group adjaceny matrix with the size of each group as
%%% 'subsize' and the probability of connection as Pin and Pbw
subsize=[300; 300; 300; 300]; % sizes of submatrixs [s1; s2; ...];
k=size(subsize,1);
n=sum(subsize);
% probability of connection in submatrix; 0<pin<1;
Pin=[0.1*ones(size(subsize))];
% Pbw=0.02*ones(size(subsize,1)); % probability of connextion between submatrixs; 0<Pbw<Pin;
% only using right up triangle;
Pbw=0.02*  [0 1 1 1 0;
    0 0 1 1 0;
    0 0 0 1 0;
    0 0 0 0 1;
    0 0 0 0 0];

% create an undirected graph adjacency matrix
Y=[];
for i=1 : size(subsize,1)
    Y=[Y;i*ones(subsize(i,1),1)];
end
W=zeros(sum(subsize));
for i=1:sum(subsize)-1
    for j=i+1:sum(subsize)
        if Y(i)==Y(j)
            if rand < Pin(Y(i))
                W(i,j)=1;
                W(j,i)=1;
            end
        else
            if rand < Pbw(Y(i),Y(j))
                W(i,j)=1;
                W(j,i)=1;
            end
        end
    end
end

label=zeros(n,k);
        p = randperm(min(subsize));
        numL=5; %number of labelled data in each class
        label(p(1:numL),1)=1; % random choose numL points in each class as labelled data
        for j=2:k
            label(sum(subsize(1:j-1))+p(1:numL),j)=1; % subnlabel labeled points of class +1
        end
        
 %%% Using our GBE to do semi-supervised learning       
 F=GBE(W,label,[150,5]);
        er=(n-nnz(F==Y))/n; % errror rate

