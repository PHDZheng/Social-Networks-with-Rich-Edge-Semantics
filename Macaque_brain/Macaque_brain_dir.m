clc
clear all
close all
path(path, 'E:\Queens\research\code\new_Laplacian');
[source,target, weight] =textread('macaque_vistact [Edges].csv', 'n%d n%d %f','headerlines',1,'delimiter', ',');
[node, name, group] =textread('macaque_vistact [Nodes].csv', 'n%d %s %s','headerlines',1,'delimiter', ',');
n=45;
ssub=[30;15];
k=2;
mycolors=[1,0,0;0,0.5,0];
colorID=zeros(n,3);
tempid=0;
for i=1:k
    colorID(tempid+1:tempid+ssub(i),:)=repmat(mycolors(i,:),ssub(i,1),1);
    tempid=tempid+ssub(i);
end
k=4;
W=sparse(source+1,target+1, weight,n,n); 

[Vout,Vin,e] = DirLaplacian(W,k);
firsteig = 2;
 while e(firsteig) < 1e-5
     firsteig = firsteig + 1;

 end

figure; hold on;
for i = 1:n
    plot3([Vout(i,firsteig),Vin(i,firsteig)],[Vout(i,firsteig+1),Vin(i,firsteig+1)],[Vout(i,firsteig+2),Vin(i,firsteig+2)],'b:');
end
for i = 1:n
    plot3(Vout(i,firsteig),Vout(i,firsteig+1),Vout(i,firsteig+2),'ko');
    plot3(Vin(i,firsteig),Vin(i,firsteig+1),Vin(i,firsteig+2),'k*');
     text(Vout(i,firsteig),Vout(i,firsteig+1),Vout(i,firsteig+2),' out','Color',colorID(i,:));
     text(Vin(i,firsteig),Vin(i,firsteig+1),Vin(i,firsteig+2),[' ' name{i} '_{in}'],'Color',colorID(i,:));
end
axis off
title('New directed Embedding');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Fan's DirectedLaplacan
[V2,e2]=DLaplacian_Fan(W,0.2,k);
firsteig = 2;
 while e2(firsteig,1) < 1e-5
     firsteig = firsteig + 1;
 end
 figure; hold on;
xlabel(firsteig);
ylabel(firsteig+1)
zlabel(firsteig+2);
for i = 1:n
    plot3(V2(i,firsteig),V2(i,firsteig+1),V2(i,firsteig+2),'*','Color',colorID(i,:));
     text(V2(i,firsteig),V2(i,firsteig+1),V2(i,firsteig+2),[' ' name{i}],'Color',colorID(i,:));
end
axis off;
title('Embedding of FanChung directed Laplacian');