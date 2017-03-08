%  expects two signed directed network adjacency matrices, both positive
clc
clear all
close all
cpt = 10;%173 nodes;
%cpt = 8;%16 nodes;
path(path, '..\new_Laplacian');
Apostemp = csvread('posdirected_Radicals.csv');
Anegtemp = csvread('negdirected_Radicals.csv');
chosenrows = csvread('chosenrows_Radicals173.csv');
chosennames = csvread('chosennames_Radicals173.csv');
% chosenrows = csvread([int2str(cpt) 'chosenrows_Radicals16.csv']);
% chosennames = csvread([int2str(cpt) 'chosennames_Radicals16.csv']);

posW = Apostemp(chosenrows,chosenrows);
negW = Anegtemp(chosenrows,chosenrows);
posW=posW-diag(diag(posW));
negW=negW-diag(diag(negW));

n=size(negW,2);
posW=posW+posW';
negW=negW+negW';

Dpos=sum(posW,2);
Dneg=sum(negW,2);
D=diag(Dpos+Dneg);
D2=zeros(n,n);
for i=1:n
    if D(i,i)~=0
        D2(i,i)=1/D(i,i)^.5;
    end
end

L=D-posW+negW; % Jerome NCut
Lsym=D2*L*D2;
Lsym=(Lsym+Lsym')/2;
[u1,e1] = eig(Lsym);
u1=D2*u1;
e=diag(e1);
[e,IXY]=sort(e,1,'ascend');
u1=u1(:,IXY);
u2=SignedLaplacian(posW,negW,'sns');
u3=SignedLaplacian(posW,negW,'bns');
%%%%%%%%%% compute AER ANR MER values
 pDinv=zeros(n,1);
    nDinv=zeros(n,1);
    for i=1:n
        if Dpos(i,1)>0;
                pDinv(i,1)=1/Dpos(i,1);
        end
        if Dneg(i,1)>0
                nDinv(i,1)=1/Dneg(i,1);
        end
    end
    posW2=posW;
    posW2(posW2>0)=1;
    negW2=negW;
    negW2(negW2>0)=1;   
    z=zeros(172,9);
for dimensions=1:172
    dis=sparse(n,n);
    for i=1:n-1
        for j=i+1:n
            if posW(i,j)~=0 || negW(i,j)~=0
                dis(i,j)=norm(u1(i,1:dimensions)-u1(j,1:dimensions));
            end
        end
    end    
    dis=dis+dis';
    AER=(sum(sum(posW.*dis)) / sum(Dpos)) / (sum(sum(negW.*dis)) / sum(Dneg));   
    ANR= (sum(sum(posW.*dis,2).*pDinv) / nnz(Dpos)) / (sum(sum(negW.*dis,2).*nDinv) / nnz(Dneg));
    MER= median( nonzeros(dis.*posW2)) / median( nonzeros(dis.*negW2));
    z(dimensions,[1 4 7])=[AER ANR MER];
    
        dis=sparse(n,n);
    for i=1:n-1
        for j=i+1:n
            if posW(i,j)~=0 || negW(i,j)~=0
                dis(i,j)=norm(u2(i,1:dimensions)-u2(j,1:dimensions));
            end
        end
    end    
    dis=dis+dis';
    AER=(sum(sum(posW.*dis)) / sum(Dpos)) / (sum(sum(negW.*dis)) / sum(Dneg));   
    ANR= (sum(sum(posW.*dis,2).*pDinv) / nnz(Dpos)) / (sum(sum(negW.*dis,2).*nDinv) / nnz(Dneg));
    MER= median( nonzeros(dis.*posW2)) / median( nonzeros(dis.*negW2));
    z(dimensions,[2 5 8])=[AER ANR MER];
    
        dis=sparse(n,n);
    for i=1:n-1
        for j=i+1:n
            if posW(i,j)~=0 || negW(i,j)~=0
                dis(i,j)=norm(u3(i,1:dimensions)-u3(j,1:dimensions));
            end
        end
    end    
    dis=dis+dis';
    AER=(sum(sum(posW.*dis)) / sum(Dpos)) / (sum(sum(negW.*dis)) / sum(Dneg));   
    ANR= (sum(sum(posW.*dis,2).*pDinv) / nnz(Dpos)) / (sum(sum(negW.*dis,2).*nDinv) / nnz(Dneg));
    MER= median( nonzeros(dis.*posW2)) / median( nonzeros(dis.*negW2));
    z(dimensions,[3 6 9])=[AER ANR MER];
end
figure;
plot(1:dimensions,z(:,1),'o-',1:dimensions,z(:,2),'*-',1:dimensions,z(:,3),'s-'), legend( '$\overline{L}_{rw}$','$L_{sns}$','$L_{bns}$',2);
h=legend;
set(h, 'interpreter', 'latex');
xlabel('# of used dimentions for embedding')
ylabel('Ratio')
figure;
plot(1:dimensions,z(:,4),'o-',1:dimensions,z(:,5),'*-',1:dimensions,z(:,6),'s-'), legend( '$\overline{L}_{rw}$','$L_{sns}$','$L_{bns}$',2);
h=legend;
set(h, 'interpreter', 'latex');
xlabel('# of used dimentions for embedding')
ylabel('Ratio')
figure;
plot(1:dimensions,z(:,7),'o-',1:dimensions,z(:,8),'*-',1:dimensions,z(:,9),'s-'), legend( '$\overline{L}_{rw}$','$L_{sns}$','$L_{bns}$',2);
h=legend;
set(h, 'interpreter', 'latex');

xlabel('# of used dimentions for embedding')
ylabel('Ratio')