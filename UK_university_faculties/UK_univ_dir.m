clc
clear all
close all
path(path, '..\new_Laplacian');
[source,target, weight] =textread('univ_dataset_TSPE [Edges].csv', 'n%d n%d %f','headerlines',1,'delimiter', ',');
[node, group] =textread('univ_dataset_TSPE [Nodes].csv', 'n%d %f','headerlines',1,'delimiter', ',');
n=size(node,1);
group=group+1;
k=max(group);
mycolors=hsv(k);
colorID=zeros(n,3);
for i=1:n
    colorID(i,:)=mycolors(group(i),:);
end
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
    plot3(Vout(i,firsteig),Vout(i,firsteig+1),Vout(i,firsteig+2),'o','Color',colorID(i,:));
    plot3(Vin(i,firsteig),Vin(i,firsteig+1),Vin(i,firsteig+2),'*','Color',colorID(i,:));
%      text(Vout(i,firsteig),Vout(i,firsteig+1),Vout(i,firsteig+2),[' out ' num2str(i)]);
%      text(Vin(i,firsteig),Vin(i,firsteig+1),Vin(i,firsteig+2),[' in '  num2str(i)]);
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
%     text(V2(i,firsteig),V2(i,firsteig+1),V2(i,firsteig+2),[num2str(i)]);
end
axis off;
title('Embedding of FanChung directed Laplacian');
