clc
clear all
close all
path(path, 'E:..\new_Laplacian'); % the path of Library

names = textread('Family_names.csv', '%s');
green = csvread('personal.csv');
red = csvread('financial.csv');
red=red - diag(diag(red));
green = green-diag(diag(green));
n = size(red,2);

%%%%%%%%%%%%%%%%%%%%%% Undirected RW embedding
red_undir=red'+red;
green_undir=green'+green;

[vbig,e]=Laplacian(red_undir,'rw');
firsteig=2;
while e(firsteig,1) < 1e-5
    firsteig = firsteig + 1;
end
figure; hold on;
for i=1:n
    for j=i+1:n
        if red_undir(i,j)~=0
            plot3(vbig(i:j-i:j,firsteig),vbig(i:j-i:j,firsteig+1),vbig(i:j-i:j,firsteig+2),'r-');
        end
    end
end
for i=1:n
    plot3(vbig(i,firsteig),vbig(i,firsteig+1),vbig(i,firsteig+2),'k*');
    if i>3 && i<13
        text(vbig(i,firsteig),vbig(i,firsteig+1),vbig(i,firsteig+2),[' ' char(names(i))],'Color', [0.5,0,0.5],'FontSize',12,'FontWeight','bold');
    elseif i==13
        text(vbig(i,firsteig),vbig(i,firsteig+1),vbig(i,firsteig+2),[' ' char(names(i))],'Color',[1,0,1],'FontSize',14,'FontWeight','bold');
    else
        text(vbig(i,firsteig),vbig(i,firsteig+1),vbig(i,firsteig+2),[' ' char(names(i))],'Color','k','FontSize',12);
    end
end
title('The embedding of the Florentine families based on financial relationships and undirected edges');
axis off;

clear e vbig;
[vbig,e]=Laplacian(green_undir,'rw');
firsteig=2;
while e(firsteig,1) < 1e-5
    firsteig = firsteig + 1;
end
figure; hold on;
for i=1:n
    for j=i+1:n
        if green_undir(i,j)~=0
            plot3(vbig(i:j-i:j,firsteig),vbig(i:j-i:j,firsteig+1),vbig(i:j-i:j,firsteig+2),'g-');
        end
    end
end
for i=1:n
    plot3(vbig(i,firsteig),vbig(i,firsteig+1),vbig(i,firsteig+2),'k*');
    if i>3 && i<13
        text(vbig(i,firsteig),vbig(i,firsteig+1),vbig(i,firsteig+2),[' ' char(names(i))],'Color', [0.5,0,0.5],'FontSize',12,'FontWeight','bold');
    elseif i==13
        text(vbig(i,firsteig),vbig(i,firsteig+1),vbig(i,firsteig+2),[' ' char(names(i))],'Color',[1,0,1],'FontSize',14,'FontWeight','bold');
    else
        text(vbig(i,firsteig),vbig(i,firsteig+1),vbig(i,firsteig+2),[' ' char(names(i))],'Color','k','FontSize',12);
    end
end
title('The embedding of the Florentine families based on personal relationships and undirected edges');
axis off;

clear e vbig;
whole=red_undir+green_undir;
[vbig,e]=Laplacian(whole,'rw');
firsteig=2;
while e(firsteig,1) < 1e-5
    firsteig = firsteig + 1;
end
figure; hold on;
for i=1:n
    for j=i+1:n
        if whole(i,j)~=0
            plot3(vbig(i:j-i:j,firsteig),vbig(i:j-i:j,firsteig+1),vbig(i:j-i:j,firsteig+2),'-','color', [0.5 0.5 0.5]);
        end
    end
end
for i=1:n
    plot3(vbig(i,firsteig),vbig(i,firsteig+1),vbig(i,firsteig+2),'k*');
    if i>3 && i<13
        text(vbig(i,firsteig),vbig(i,firsteig+1),vbig(i,firsteig+2),[' ' char(names(i))],'Color', [0.5,0,0.5],'FontSize',12,'FontWeight','bold');
    elseif i==13
        text(vbig(i,firsteig),vbig(i,firsteig+1),vbig(i,firsteig+2),[' ' char(names(i))],'Color',[1,0,1],'FontSize',14,'FontWeight','bold');
    else
        text(vbig(i,firsteig),vbig(i,firsteig+1),vbig(i,firsteig+2),[' ' char(names(i))],'Color','k','FontSize',12);
    end
end
title('The embedding of the Florentine families of undirected edges ignoring edge types');
axis off;

%%%%%%%%%%%%%%%%% our new directed embedding
clear e Vout Vin whole;
whole=red+green;
[Vout,Vin,e] = DirLaplacian(whole);
firsteigeig=2;
while e(firsteigeig,1) < 1e-5
    firsteigeig = firsteigeig + 1;
end
figure; hold on;
xlabel(int2str(firsteigeig));
ylabel(int2str(firsteigeig+1));
zlabel(int2str(firsteigeig+2));
for i=1:n
    for j=1:n
        if whole(i,j)~=0
            plot3([Vout(j,firsteigeig),Vin(i,firsteigeig)],[Vout(j,firsteigeig+1),Vin(i,firsteigeig+1)],[Vout(j,firsteigeig+2),Vin(i,firsteigeig+2)],'-','color', [0.5 0.5 0.5],'LineWidth',1);
        end
    end
end
for i= 1:n
    plot3([Vout(i,firsteigeig),Vin(i,firsteigeig)],[Vout(i,firsteigeig+1),Vin(i,firsteigeig+1)],[Vout(i,firsteigeig+2),Vin(i,firsteigeig+2)],'k:','LineWidth',3);
end
for i=1:n
    plot3(Vin(i,firsteigeig),Vin(i,firsteigeig+1),Vin(i,firsteigeig+2),'k*');
    plot3(Vout(i,firsteigeig),Vout(i,firsteigeig+1),Vout(i,firsteigeig+2),'ko');
    if i>3 && i<13
        text(Vout(i,firsteigeig),Vout(i,firsteigeig+1),Vout(i,firsteigeig+2),[' out'],'Color', [0.5 0 0.5],'FontSize',10,'FontWeight','bold');
        text(Vin(i,firsteigeig),Vin(i,firsteigeig+1),Vin(i,firsteigeig+2),[' ', char(names(i)) '_{in}'],'Color', [0.5 0 0.5],'FontSize',10,'FontWeight','bold');
    elseif i==13
        text(Vout(i,firsteigeig),Vout(i,firsteigeig+1),Vout(i,firsteigeig+2),[' out'],'Color',[1,0,1],'FontSize',12,'FontWeight','bold');
        text(Vin(i,firsteigeig),Vin(i,firsteigeig+1),Vin(i,firsteigeig+2),[' ', char(names(i)) '_{in}'],'Color',[1,0,1],'FontSize',12,'FontWeight','bold');
    else
        text(Vout(i,firsteigeig),Vout(i,firsteigeig+1),Vout(i,firsteigeig+2),[' out'],'Color',[0.1,0.1,0.1],'FontSize',10);
        text(Vin(i,firsteigeig),Vin(i,firsteigeig+1),Vin(i,firsteigeig+2),[' ', char(names(i)) '_{in}'],'Color',[0.1,0.1,0.1],'FontSize',10);
    end
end
title('The embedding of the Florentine families based on our new directed approach');
axis off;

%%%%%%%%%%%%%% / normalized edge lengths of the directed embedding for Florentine families
k=3;
mtp=1000;
Din=sum(whole,1)';
Dout=sum(whole,2);
NormalizedInOut=zeros(n,1);
for i=1:n
    NormalizedInOut(i,1)=norm(Vout(i,firsteigeig:k+firsteigeig-1)-Vin(i,firsteigeig:k+firsteigeig-1))*(Din(i,1)+Dout(i,1)); % normalized edge= distance * edge weight
end

NormalizedDis=zeros(n,n);
for i=1:n
    for j=1:n
        if whole(i,j)~=0
            NormalizedDis(i,j)=norm(Vout(i,firsteigeig:k+firsteigeig-1)-Vin(j,firsteigeig:k+firsteigeig-1))*whole(i,j);
        end
    end
end
AvgNorDis=zeros(n,1);
for i=1:n
    AvgNorDis(i,1)=(sum(NormalizedDis(i,:))+sum(NormalizedDis(:,i))) / (nnz(NormalizedDis(i,:))+nnz(NormalizedDis(:,i)));
end
%     %%%%% write a latex table
%     for i=1:n
%         nodecopies(i,1)={[names{i,1}  ' & ' num2str(int16(NormalizedInOut(i,1)*mtp)) ' & ' num2str(int16(AvgNorDis(i,1)*mtp)) ' & ' num2str(Din(i,1)) ' & ' num2str(Dout(i,1)) ' \\']};
%     end
%     fileID = fopen('table_dir.txt','w');
%     fprintf(fileID,'%s \n','Names  &  in-out  &  average & Din & Dout\\');
%     fprintf(fileID,'%s \n','\hline');
%     for i=1:n
%         fprintf(fileID,'%s \n',nodecopies{i,1:end});
%     end
%     fprintf(fileID,'%s \n',['Means & ' num2str(mean(NormalizedInOut(:,1)*mtp)) ' & ' num2str(mean(AvgNorDis(:,1)*mtp)) ' & ' num2str(mean(Din)) ' & ' num2str(mean(Dout)) ' \\']);
%     fprintf(fileID,'%s \n',['Std & ' num2str(std(NormalizedInOut(:,1)*mtp)) ' & ' num2str(std(AvgNorDis(:,1)*mtp)) ' & ' num2str(std(Din)) ' & ' num2str(std(Dout)) ' \\']);
%     fclose(fileID);

%%%%%%%%%%%%%%%%% FangChung's directed embedding
clear e vbig;
[vbig,e]=DLaplacian_Fan(whole,0.2,4);
firsteig=2;
while e(firsteig,1) < 1e-5
    firsteig = firsteig + 1;
end
figure; hold on;
for i=1:n
    for j=1:n
        if whole(i,j)~=0
            plot3(vbig(i:j-i:j,firsteig),vbig(i:j-i:j,firsteig+1),vbig(i:j-i:j,firsteig+2),'-','color', [0.5 0.5 0.5]);
        end
    end
end
for i=1:n
    plot3(vbig(i,firsteig),vbig(i,firsteig+1),vbig(i,firsteig+2),'k*');
    if i>3 && i<13
        text(vbig(i,firsteig),vbig(i,firsteig+1),vbig(i,firsteig+2),[' ' char(names(i))],'Color', [0.5 0 0.5],'FontSize',10,'FontWeight','bold');
    elseif i==13
        text(vbig(i,firsteig),vbig(i,firsteig+1),vbig(i,firsteig+2),[' ' char(names(i))],'Color',[1,0,1],'FontSize',12,'FontWeight','bold');
    else
        text(vbig(i,firsteig),vbig(i,firsteig+1),vbig(i,firsteig+2),[' ' char(names(i))],'Color',[0.1,0.1,0.1],'FontSize',10);
    end
end
title('The embedding of the Florentine families based on FangChung directed approach');
axis off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% undirected typed embedding
clear e vbig;
big_undir(:,:,1)=red_undir;
big_undir(:,:,2)=green_undir;
[vbig,e] = TypedLaplacian(big_undir,4);
firsteig=2;
figure; hold on;
xlabel(int2str(firsteig));
ylabel(int2str(firsteig+1));
zlabel(int2str(firsteig+2));
for i=1:n
    for j=i+1:n
        if red_undir(i,j)~=0
            plot3(vbig(i:j-i:j,firsteig,1),vbig(i:j-i:j,firsteig+1,1),vbig(i:j-i:j,firsteig+2,1),'r-');
        end
        if green_undir(i,j)~=0
            plot3(vbig(i:j-i:j,firsteig,2),vbig(i:j-i:j,firsteig+1,2),vbig(i:j-i:j,firsteig+2,2),'g-');
        end
    end
end
for i= 1:n
    plot3(reshape(vbig(i,firsteig,:),2,1),reshape(vbig(i,firsteig+1,:),2,1),reshape(vbig(i,firsteig+2,:),2,1),'b:')
end
for i=1:n
    plot3(vbig(i,firsteig,1),vbig(i,firsteig+1,1),vbig(i,firsteig+2,1),'k*');
    plot3(vbig(i,firsteig,2),vbig(i,firsteig+1,2),vbig(i,firsteig+2,2),'ko');
    if i>3 && i<13
        text(vbig(i,firsteig,1),vbig(i,firsteig+1,1),vbig(i,firsteig+2,1),[' ' char(names(i))],'Color', [0.5 0 0.5],'FontSize',10,'FontWeight','bold');
        text(vbig(i,firsteig,2),vbig(i,firsteig+1,2),vbig(i,firsteig+2,2),[' ' lower(char(names(i)))],'Color', [0.5 0 0.5],'FontSize',10,'FontWeight','bold');
    elseif i==13
        text(vbig(i,firsteig,1),vbig(i,firsteig+1,1),vbig(i,firsteig+2,1),[' ' char(names(i))],'Color',[1,0,1],'FontSize',12,'FontWeight','bold');
        text(vbig(i,firsteig,2),vbig(i,firsteig+1,2),vbig(i,firsteig+2,2),[' ' lower(char(names(i)))],'Color',[1,0,1],'FontSize',12,'FontWeight','bold');
    else
        text(vbig(i,firsteig,1),vbig(i,firsteig+1,1),vbig(i,firsteig+2,1),[' ' char(names(i))],'Color',[0.1,0.1,0.1],'FontSize',10);
        text(vbig(i,firsteig,2),vbig(i,firsteig+1,2),vbig(i,firsteig+2,2),[' ' lower(char(names(i)))],'Color',[0.1,0.1,0.1],'FontSize',10);
    end
end
title('The undirected typed embedding of the Florentine families');
axis off;

%%%%%%%%%%%%%%%%%% directed typed embedding with FanChung's directed Laplacian
clear e vbig;
big(:,:,1)=red;
big(:,:,2)=green;
[vbig,e] = TypedLaplacian(big,4,0.2);
firsteig=2;
figure; hold on;
xlabel(int2str(firsteig));
ylabel(int2str(firsteig+1));
zlabel(int2str(firsteig+2));
for i=1:n
    for j=i+1:n
        if red_undir(i,j)~=0
            plot3(vbig(i:j-i:j,firsteig,1),vbig(i:j-i:j,firsteig+1,1),vbig(i:j-i:j,firsteig+2,1),'r-');
        end
        if green_undir(i,j)~=0
            plot3(vbig(i:j-i:j,firsteig,2),vbig(i:j-i:j,firsteig+1,2),vbig(i:j-i:j,firsteig+2,2),'g-');
        end
    end
end
for i= 1:n
    plot3(reshape(vbig(i,firsteig,:),2,1),reshape(vbig(i,firsteig+1,:),2,1),reshape(vbig(i,firsteig+2,:),2,1),'b:')
end
for i=1:n
    plot3(vbig(i,firsteig,1),vbig(i,firsteig+1,1),vbig(i,firsteig+2,1),'k*');
    plot3(vbig(i,firsteig,2),vbig(i,firsteig+1,2),vbig(i,firsteig+2,2),'ko');
    if i>3 && i<13
        text(vbig(i,firsteig,1),vbig(i,firsteig+1,1),vbig(i,firsteig+2,1),[' ' char(names(i))],'Color', [0.5 0 0.5],'FontSize',10,'FontWeight','bold');
        text(vbig(i,firsteig,2),vbig(i,firsteig+1,2),vbig(i,firsteig+2,2),[' ' lower(char(names(i)))],'Color', [0.5 0 0.5],'FontSize',10,'FontWeight','bold');
    elseif i==13
        text(vbig(i,firsteig,1),vbig(i,firsteig+1,1),vbig(i,firsteig+2,1),[' ' char(names(i))],'Color',[1,0,1],'FontSize',12,'FontWeight','bold');
        text(vbig(i,firsteig,2),vbig(i,firsteig+1,2),vbig(i,firsteig+2,2),[' ' lower(char(names(i)))],'Color',[1,0,1],'FontSize',12,'FontWeight','bold');
    else
        text(vbig(i,firsteig,1),vbig(i,firsteig+1,1),vbig(i,firsteig+2,1),[' ' char(names(i))],'Color',[0.1,0.1,0.1],'FontSize',10);
        text(vbig(i,firsteig,2),vbig(i,firsteig+1,2),vbig(i,firsteig+2,2),[' ' lower(char(names(i)))],'Color',[0.1,0.1,0.1],'FontSize',10);
    end
end
title('The directed typed embedding of the Florentine families with FanChung directed Laplacian');
axis off;

%%%%%%%%%%%%%%%%% our new typed directed embedding
clear e Vout Vin;
[Vout,Vin,e, bigW] = TypedDirLaplacian(big,4);
firsteigeig=2;
while e(firsteigeig,1) < 1e-5
    firsteigeig = firsteigeig + 1;
end
figure; hold on;
xlabel(int2str(firsteigeig));
ylabel(int2str(firsteigeig+1));
zlabel(int2str(firsteigeig+2));
for i=1:n
    for j=1:n
        if red(i,j)~=0
            plot3([Vout(j,firsteigeig,1),Vin(i,firsteigeig,1)],[Vout(j,firsteigeig+1,1),Vin(i,firsteigeig+1,1)],[Vout(j,firsteigeig+2,1),Vin(i,firsteigeig+2,1)],'r-');
        end
        if green(i,j)~=0
            plot3([Vout(j,firsteigeig,2),Vin(i,firsteigeig,2)],[Vout(j,firsteigeig+1,2),Vin(i,firsteigeig+1,2)],[Vout(j,firsteigeig+2,2),Vin(i,firsteigeig+2,2)],'g-');
        end
    end
end
for i= 1:n
    plot3([Vout(i,firsteigeig,1),Vin(i,firsteigeig,1)],[Vout(i,firsteigeig+1,1),Vin(i,firsteigeig+1,1)],[Vout(i,firsteig+2,1),Vin(i,firsteig+2,1)],'r--');
    plot3([Vout(i,firsteig,2),Vin(i,firsteig,2)],[Vout(i,firsteig+1,2),Vin(i,firsteig+1,2)],[Vout(i,firsteig+2,2),Vin(i,firsteig+2,2)],'g--');
    plot3([Vout(i,firsteig,1),Vin(i,firsteig,2)],[Vout(i,firsteig+1,1),Vin(i,firsteig+1,2)],[Vout(i,firsteig+2,1),Vin(i,firsteig+2,2)],'b:');
    plot3([Vout(i,firsteig,2),Vin(i,firsteig,1)],[Vout(i,firsteig+1,2),Vin(i,firsteig+1,1)],[Vout(i,firsteig+2,2),Vin(i,firsteig+2,1)],'b:');
end
for i=1:n
    plot3(Vin(i,firsteig,1),Vin(i,firsteig+1,1),Vin(i,firsteig+2,1),'r*');
    plot3(Vout(i,firsteig,1),Vout(i,firsteig+1,1),Vout(i,firsteig+2,1),'ro');
    plot3(Vin(i,firsteig,2),Vin(i,firsteig+1,2),Vin(i,firsteig+2,2),'g*');
    plot3(Vout(i,firsteig,2),Vout(i,firsteig+1,2),Vout(i,firsteig+2,2),'go');
    if i>3 && i<13
        text(Vout(i,firsteig,1),Vout(i,firsteig+1,1),Vout(i,firsteig+2,1),[' OUT'],'Color', [0.5 0 0.5],'FontSize',10,'FontWeight','bold');
        text(Vin(i,firsteig,1),Vin(i,firsteig+1,1),Vin(i,firsteig+2,1),[' ', char(names(i)) '_{IN}'],'Color', [0.5 0 0.5],'FontSize',10,'FontWeight','bold');
        text(Vout(i,firsteig,2),Vout(i,firsteig+1,2),Vout(i,firsteig+2,2),[' out'],'Color', [0.5 0 0.5],'FontSize',10,'FontWeight','bold');
        text(Vin(i,firsteig,2),Vin(i,firsteig+1,2),Vin(i,firsteig+2,2),[' ', lower(char(names(i))) '_{in}'],'Color', [0.5 0 0.5],'FontSize',10,'FontWeight','bold');
    elseif i==13
        text(Vout(i,firsteig,1),Vout(i,firsteig+1,1),Vout(i,firsteig+2,1),[' OUT'],'Color',[1,0,1],'FontSize',12,'FontWeight','bold');
        text(Vin(i,firsteig,1),Vin(i,firsteig+1,1),Vin(i,firsteig+2,1),[' ', char(names(i)) '_{IN}'],'Color',[1,0,1],'FontSize',12,'FontWeight','bold');
        text(Vout(i,firsteig,2),Vout(i,firsteig+1,2),Vout(i,firsteig+2,2),[' out'],'Color',[1,0,1],'FontSize',12,'FontWeight','bold');
        text(Vin(i,firsteig,2),Vin(i,firsteig+1,2),Vin(i,firsteig+2,2),[' ', lower(char(names(i))) '_{in}'],'Color',[1,0,1],'FontSize',12,'FontWeight','bold');
    else
        text(Vout(i,firsteig,1),Vout(i,firsteig+1,1),Vout(i,firsteig+2,1),[' OUT'],'Color',[0.1,0.1,0.1],'FontSize',10);
        text(Vin(i,firsteig,1),Vin(i,firsteig+1,1),Vin(i,firsteig+2,1),[' ', char(names(i)) '_{IN}'],'Color',[0.1,0.1,0.1],'FontSize',10);
        text(Vout(i,firsteig,2),Vout(i,firsteig+1,2),Vout(i,firsteig+2,2),[' out'],'Color',[0.1,0.1,0.1],'FontSize',10);
        text(Vin(i,firsteig,2),Vin(i,firsteig+1,2),Vin(i,firsteig+2,2),[' ', lower(char(names(i))) '_{in}'],'Color',[0.1,0.1,0.1],'FontSize',10);
    end
end
title('The directed typed embedding of the Florentine families with our new directed approach');
axis off;

%%%%%%%%%%%%%% / normalized edge lengths of the typed directed embedding for Florentine families
k=3;
mtp=10000;
%     Din=sum(whole,1)';
%     Dout=sum(whole,2);
bigIn=sum(bigW,1)';
bigOut=sum(bigW,2);
FF=zeros(n,1);
PP=zeros(n,1);
FP=zeros(n,1);
PF=zeros(n,1);
for i=1:n% normalized edge= distance * edge weight
    FF(i,1)=norm(Vout(i,firsteig:k+firsteig-1,1)-Vin(i,firsteig:k+firsteig-1,1))*(bigIn(i,1)+bigOut(i,1)+bigW(i,i));
    PP(i,1)=norm(Vout(i,firsteig:k+firsteig-1,2)-Vin(i,firsteig:k+firsteig-1,2))*(bigIn(i+n,1)+bigOut(i+n,1)+bigW(i,i));
    FP(i,1)=norm(Vout(i,firsteig:k+firsteig-1,1)-Vin(i,firsteig:k+firsteig-1,2))*bigW(i,i+n);
    PF(i,1)=norm(Vout(i,firsteig:k+firsteig-1,2)-Vin(i,firsteig:k+firsteig-1,1))*bigW(i+n,i);
end

NormalizedDis=zeros(2*n,2*n);
for i=1:n
    for j=1:n
        if bigW(i,j)~=0 && i~=j
            NormalizedDis(i,j)=norm(Vout(i,firsteig:k+firsteig-1,1)-Vin(j,firsteig:k+firsteig-1,1))*bigW(i,j);
        end
        if bigW(i+n,j+n)~=0 && i~=j
            NormalizedDis(i+n,j+n)=norm(Vout(i,firsteig:k+firsteig-1,2)-Vin(j,firsteig:k+firsteig-1,2))*bigW(i+n,j+n);
        end
    end
end
AvgNorDis=zeros(n,1);
for i=1:n
    AvgNorDis(i,1)=(sum(NormalizedDis(i,1:n))+sum(NormalizedDis(1:n,i))+sum(NormalizedDis(i+n,n+1:2*n))+sum(NormalizedDis(n+1:2*n,i+n))) / (nnz(NormalizedDis(i,1:n))+nnz(NormalizedDis(1:n,i))+nnz(NormalizedDis(i+n,n+1:2*n))+nnz(NormalizedDis(n+1:2*n,i+n)));
end
    % %%%%% write a latex table
    % for i=1:n
    %     nodecopies(i,1)={[names{i,1}  ' & ' num2str(int16(FF(i,1)*mtp)) ' & ' num2str(int16(PP(i,1)*mtp)) ' & ' num2str(int16(FP(i,1)*mtp)) ' & ' num2str(int16(PF(i,1)*mtp)) ' & ' num2str(int16(AvgNorDis(i,1)*mtp)) ' & ' num2str(Din(i,1)) ' & ' num2str(Dout(i,1)) ' \\']};
    % end
    % fileID = fopen('table_typed_dir.txt','w');
    % fprintf(fileID,'%s \n','Names  &  FF  &  PP  &  FP  &  PF  & Avg & Din & Dout\\');
    % fprintf(fileID,'%s \n','\hline');
    % for i=1:n
    %     fprintf(fileID,'%s \n',nodecopies{i,1:end});
    % end
    % fprintf(fileID,'%s \n',['Means & ' num2str(mean(FF*mtp)) ' & ' num2str(mean(PP*mtp)) ' & ' num2str(mean(FP*mtp)) ' & ' num2str(mean(PF*mtp)) ' & ' num2str(mean(AvgNorDis*mtp)) ' & ' num2str(mean(Din)) ' & ' num2str(mean(Dout)) ' \\']);
    % fprintf(fileID,'%s \n',['Std & ' num2str(std(FF*mtp)) ' & ' num2str(std(PP*mtp)) ' & ' num2str(std(FP*mtp)) ' & ' num2str(std(PF*mtp)) ' & ' num2str(std(AvgNorDis*mtp)) ' & ' num2str(std(Din)) ' & ' num2str(std(Dout)) ' \\']);
    % fclose(fileID);



