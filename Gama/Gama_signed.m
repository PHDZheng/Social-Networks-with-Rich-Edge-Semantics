clc
clear all
close all
path(path, '..\new_Laplacian');
posW = csvread('gamaPos.csv');
negW = csvread('gamaNeg.csv');
names = textread('name.csv', '%s');

n=size(posW,2);
Dpos=sum(posW,2);
Dneg=sum(negW,2);
D=diag(Dpos+Dneg);
D2=zeros(n,n);
for i=1:n
    if D(i,i)~=0
        D2(i,i)=1/D(i,i)^.5;
    end
end

L=D-posW+negW;
Lsym=D2*L*D2;
[u1,e] = eig(Lsym);
u1=D2*u1;
figure; hold on;
    for i = 1:n
        for j = i+1:n
            if posW(i,j) ~= 0
                plot3(u1(i:j-i:j,1),u1(i:j-i:j,2),u1(i:j-i:j,3),'g-');
            elseif negW(i,j) ~= 0
                plot3(u1(i:j-i:j,1),u1(i:j-i:j,2),u1(i:j-i:j,3),'r--');
            end
        end
    end
    for i = 3:14
        plot3(u1(i,1),u1(i,2),u1(i,3),'k*');
            text(u1(i,1),u1(i,2),u1(i,3),[char(names(i)) ' '],...
                'HorizontalAlignment','right');
    end
        for i = [1 2 15 16]
        plot3(u1(i,1),u1(i,2),u1(i,3),'k*');
            text(u1(i,1),u1(i,2),u1(i,3),[' ' char(names(i))]);
    end
title('Embedding of jerome NCut');
axis off;

u2=SignedLaplacian(posW,negW,'sns',3);
figure; hold on;axis on;
    for i = 1:n
        for j = i+1:n
            if posW(i,j) ~= 0
                plot3(u2(i:j-i:j,1),u2(i:j-i:j,2),u2(i:j-i:j,3),'g-');
            elseif negW(i,j) ~= 0
                plot3(u2(i:j-i:j,1),u2(i:j-i:j,2),u2(i:j-i:j,3),'r--');
            end
        end
    end
    for i = 3:14
        plot3(u2(i,1),u2(i,2),u2(i,3),'k*');
            text(u2(i,1),u2(i,2),u2(i,3),[char(names(i)) ' '],...
                'HorizontalAlignment','right');
    end
        for i = [1 2 15 16]
        plot3(u2(i,1),u2(i,2),u2(i,3),'k*');
            text(u2(i,1),u2(i,2),u2(i,3),[' ' char(names(i))]);
    end
title('Embedding of SNScut');
axis off;

u3=SignedLaplacian(posW,negW,'bns',3);
figure; hold on;axis on;
    for i = 1:n
        for j = i+1:n
            if posW(i,j) ~= 0
                plot3(u3(i:j-i:j,1),u3(i:j-i:j,2),u3(i:j-i:j,3),'g-');
            elseif negW(i,j) ~= 0
                plot3(u3(i:j-i:j,1),u3(i:j-i:j,2),u3(i:j-i:j,3),'r--');
            end
        end
    end
for i = 3:14
        plot3(u3(i,1),u3(i,2),u3(i,3),'k*');
            text(u3(i,1),u3(i,2),u3(i,3),[char(names(i)) ' '],...
                'HorizontalAlignment','right');
    end
        for i = [1 2 15 16]
        plot3(u3(i,1),u3(i,2),u3(i,3),'k*');
            text(u3(i,1),u3(i,2),u3(i,3),[' ' char(names(i))]);
    end
%title('Embedding of our BNScut');
axis off;

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
dimensions=3;
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
    z(1,:)=[AER ANR MER];
    
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
    z(2,:)=[AER ANR MER];
    
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
    z(3,:)=[AER ANR MER];
    z=full(z);
