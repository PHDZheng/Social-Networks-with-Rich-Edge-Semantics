clc
clear all
close all
path(path, '..\new_Laplacian');
pos = csvread('SAMPLK3.csv');
neg = csvread('SAMPDLK.csv');
% pos = csvread('SAMPES.csv');
% neg = csvread('SAMPDES.csv');
names = textread('name.csv', '%s');
classes=['k';'b';'b';'b';'b';'b';'k';'m';'m';'m';'m';'m';'m';'m';'k';'c';'c';'c';];
%%%%% Sampson delineated: m-Young Turks, c-Loyalists, b-Outcasts,
%%%%% k-interstitial group. The social cleavage between the young turks and
%%%%% the Loyalists is evident.
pos=pos(:,1:18);
neg=neg(:,1:18);
%pos(pos>=1)=1;
%neg(neg>=1)=1;
posW=pos+pos';
negW=neg+neg';
n=size(posW,2);


%%%%%%%%% Embedding of jerome NCut
Dpos=sum(posW,2);
Dneg=sum(negW,2);
D=diag(Dpos+Dneg);
D2=zeros(n,n);
for i=1:n
    if D(i,i)~=0
        D2(i,i)=1/D(i,i)^.5;
    end
end
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
    
L=D-posW+negW;
Lsym=D2*L*D2;
[u,e] = eig(Lsym);
u=D2*u;
figure; hold on;axis on;
    for i = 1:n
        for j = i+1:n
            if posW(i,j) ~= 0
                plot3(u(i:j-i:j,1),u(i:j-i:j,2),u(i:j-i:j,3),'g-');
            elseif negW(i,j) ~= 0
                plot3(u(i:j-i:j,1),u(i:j-i:j,2),u(i:j-i:j,3),'r--');
            end
        end
    end
    for i = 1:n
        plot3(u(i,1),u(i,2),u(i,3),'k*');
        text(u(i,1),u(i,2),u(i,3),[' ' char(names(i))],'Color',classes(i,:));
    end
title('Embedding of jerome NCut');
axis off;

dimensions=3;
    dis=sparse(n,n);
    for i=1:n-1
        for j=i+1:n
            if posW(i,j)~=0 || negW(i,j)~=0
                dis(i,j)=norm(u(i,1:dimensions)-u(j,1:dimensions));
            end
        end
    end    
    dis=dis+dis';
    AER=(sum(sum(posW.*dis)) / sum(Dpos)) / (sum(sum(negW.*dis)) / sum(Dneg));   
    ANR= (sum(sum(posW.*dis,2).*pDinv) / nnz(Dpos)) / (sum(sum(negW.*dis,2).*nDinv) / nnz(Dneg));
    MER= median( nonzeros(dis.*posW2)) / median( nonzeros(dis.*negW2));
    z(1,:)=[AER ANR MER];
%%%%%%% The new SNS signed Laplacian embedding
u=SignedLaplacian(posW,negW,'sns');

figure; hold on;axis on;
    for i = 1:n
        for j = i+1:n
            if posW(i,j) ~= 0
                plot3(u(i:j-i:j,1),u(i:j-i:j,2),u(i:j-i:j,3),'g-');
            elseif negW(i,j) ~= 0
                plot3(u(i:j-i:j,1),u(i:j-i:j,2),u(i:j-i:j,3),'r--');
            end
        end
    end
    for i = 1:n
        plot3(u(i,1),u(i,2),u(i,3),'k*');
        text(u(i,1),u(i,2),u(i,3),[' ' char(names(i))],'Color',classes(i,:));
    end
title('Embedding of SNS');
axis off;
    dis=sparse(n,n);
    for i=1:n-1
        for j=i+1:n
            if posW(i,j)~=0 || negW(i,j)~=0
                dis(i,j)=norm(u(i,1:dimensions)-u(j,1:dimensions));
            end
        end
    end    
    dis=dis+dis';
    AER=(sum(sum(posW.*dis)) / sum(Dpos)) / (sum(sum(negW.*dis)) / sum(Dneg));   
    ANR= (sum(sum(posW.*dis,2).*pDinv) / nnz(Dpos)) / (sum(sum(negW.*dis,2).*nDinv) / nnz(Dneg));
    MER= median( nonzeros(dis.*posW2)) / median( nonzeros(dis.*negW2));
    z(2,:)=[AER ANR MER];
%%%%%%% The new BNS signed Laplacian embedding
u=SignedLaplacian(posW,negW,'bns');

figure; hold on;axis on;
    for i = 1:n
        for j = i+1:n
            if posW(i,j) ~= 0
                plot3(u(i:j-i:j,1),u(i:j-i:j,2),u(i:j-i:j,3),'g-');
            elseif negW(i,j) ~= 0
                plot3(u(i:j-i:j,1),u(i:j-i:j,2),u(i:j-i:j,3),'r--');
            end
        end
    end
    for i = 1:n
        plot3(u(i,1),u(i,2),u(i,3),'k*');
        text(u(i,1),u(i,2),u(i,3),[' ' char(names(i))],'Color',classes(i,:));
    end
title('Embedding of BNS');
axis off;
    dis=sparse(n,n);
    for i=1:n-1
        for j=i+1:n
            if posW(i,j)~=0 || negW(i,j)~=0
                dis(i,j)=norm(u(i,1:dimensions)-u(j,1:dimensions));
            end
        end
    end    
    dis=dis+dis';
    AER=(sum(sum(posW.*dis)) / sum(Dpos)) / (sum(sum(negW.*dis)) / sum(Dneg));   
    ANR= (sum(sum(posW.*dis,2).*pDinv) / nnz(Dpos)) / (sum(sum(negW.*dis,2).*nDinv) / nnz(Dneg));
    MER= median( nonzeros(dis.*posW2)) / median( nonzeros(dis.*negW2));
    z(3,:)=[AER ANR MER];
    z=full(z);

%%%%%%%%%%%%%%%%%  the new directed part 4 copies
%%%%%%%%%%%%%%% The two belowing outputs are the same result, where the 
%%%%%%%%%%%%%%% first "newv" is  a 4n*4n matrix, and the second output is
%%%%%%%%%%%%%%% four n*n martix.
[newv,e] = SignedDirLaplacian(pos,neg,'sns',3);
[PosOut, Negout, PosIn, NegIn] = SignedDirLaplacian(pos,neg,'sns',3);

figure;
hold on;
axis off;
plot(newv(1:n,1),newv(1:n,2),'ko');
plot(newv(2*n+1:3*n,1),newv(2*n+1:3*n,2),'k*');
for i = 1:n
    text(newv(i,1),newv(i,2),[' ' char(names(i))],'Color',classes(i,:));
    plot(newv(i:2*n:2*n+i,1),newv(i:2*n:2*n+i,2),'b:','LineWidth',2); 
end
title('Positive part In and Out');
hold off;

figure;
hold on;
axis off;
for i = 1:n
    plot(PosOut(i,1),PosOut(i,2),'o','Color',classes(i,:));
plot(PosIn(i,1),PosIn(i,2),'*','Color',classes(i,:));
    text(PosOut(i,1),PosOut(i,2),[' ' char(names(i))],'Color',classes(i,:));
    plot([PosOut(i,1),PosIn(i,1)],[PosOut(i,2),PosIn(i,2)],'b:','LineWidth',2); 
end
title('Positive part In and Out');
hold off;

figure;
hold on;
axis off;
plot(newv(n+1:2*n,1),newv(n+1:2*n,2),'ko');
plot(newv(3*n+1:4*n,1),newv(3*n+1:4*n,2),'k*');
for i = n+1:2*n
    text(newv(i,1),newv(i,2),[' ' char(names(i-n))],'Color',classes(i-n,:));
    plot(newv(i:2*n:2*n+i,1),newv(i:2*n:2*n+i,2),'b:','LineWidth',2); 
end
title('Negative part In and Out');
hold off;


figure;
hold on;
axis off;
for i = 1:n
    plot(Negout(i,1),Negout(i,2),'o','Color',classes(i,:));
plot(NegIn(i,1),NegIn(i,2),'*','Color',classes(i,:));
    text(Negout(i,1),Negout(i,2),[' ' char(names(i))],'Color',classes(i,:));
    plot([Negout(i,1),NegIn(i,1)],[Negout(i,2),NegIn(i,2)],'b:','LineWidth',2); 
end
title('Negative part In and Out');
hold off;

%%%%%%%% Product of edge length and reciprocal of edge weight; deviations
%%%%%%%% from average indicate nodes with unusual neighborhoods
ds = newv(:,1:2);
for i = 1:n
	% pos in to pos out
	dist(i,1) = norm(ds(2*n+i,:) - ds(i,:),2);
	% neg in to neg out
	dist(i,2) = norm(ds(3*n+i,:) - ds(n+i,:),2);
	% pos in to neg out
	dist(i,3) = norm(ds(2*n+i,:) - ds(n+i,:),2);
	% neg in to pos out
	dist(i,4) = norm(ds(3*n+i,:) - ds(i,:),2);
end
D=diag(sum(posW+negW));
dist2=D*dist; %normalized distance

figure;
imagesc(dist2(:,1:4));
title('Pos-pos; Neg-neg, Pos-neg; Neg-pos');
