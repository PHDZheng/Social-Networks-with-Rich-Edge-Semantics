clc
clear all
close all
posW = csvread('Libya_posdirected.csv');
negW = csvread('Libya_negdirected.csv');
%name=importdata('Libya_names.txt');
name=importdata('Libya_shortnames.txt');

%  posW(posW>0)=1;
%  negW(negW>0)=1;
x=sum(posW+posW',2);
posW=posW(x>0,x>0);
negW=negW(x>0,x>0);
name=name(x>0,1);

tempW=posW+negW+posW'+negW'+eye(size(posW,2));
x=tempW(1,:);
for i=1:10
    x=x*tempW;
end
posW=posW(x>0,x>0);
negW=negW(x>0,x>0);
name=name(x>0,1);
n=size(posW,1);
path(path, '..\new_Laplacian');

[U,e] = SignedLaplacian(posW+posW',negW+negW','sns',3);% undirected
figure; hold on;
    for i = 1:n
        for j = i+1:n
            if posW(i,j) ~= 0
                plot3(U(i:j-i:j,1),U(i:j-i:j,2),U(i:j-i:j,3),'g-');
            elseif negW(i,j) ~= 0
                plot3(U(i:j-i:j,1),U(i:j-i:j,2),U(i:j-i:j,3),'r-');
            end
        end
    end
    for i = 1:n
        plot3(U(i,1),U(i,2),U(i,3),'ko');
            %text(U(i,1),U(i,2),U(i,3),[' ' int2str(i)]);
            text(U(i,1),U(i,2),U(i,3),[' ' name(i)]);
    end
axis off;
title('Undirected signed (SNS)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[newv,e] = SignedDirLaplacian(posW,negW,'sns',3);
figure;
hold on;
axis off;
for i = 1:n
	text(newv(i+n,1),newv(i+n,2),newv(i+n,3),[' ' name(i)]);
end
plot3(newv(1:n,1),newv(1:n,2),newv(1:n,3),'g^');
plot3(newv(2*n+1:3*n,1),newv(2*n+1:3*n,2),newv(2*n+1:3*n,3),'g*');
plot3(newv(n+1:2*n,1),newv(n+1:2*n,2),newv(n+1:2*n,3),'r^'); 
plot3(newv(3*n+1:4*n,1),newv(3*n+1:4*n,2),newv(3*n+1:4*n,3),'r*');
for i = 1:n
    for j = 1:n
        if posW(i,j) ~= 0
            plot3(newv(i:j-i+2*n:j+2*n,1),newv(i:j-i+2*n:j+2*n,2),newv(i:j-i+2*n:j+2*n,3),'g-');
        end

    end
end
for i = 1:n
	plot3(newv(i:2*n:2*n+i,1),newv(i:2*n:2*n+i,2),newv(i:2*n:2*n+i,3),'b--','LineWidth',2); 
end 
for i = 1:n
    for j = 1:n
        if negW(i,j) ~= 0
            plot3(newv(i+n:j-i+2*n:j+3*n,1),newv(i+n:j-i+2*n:j+3*n,2),newv(i+n:j-i+2*n:j+3*n,3),'r-');
        end
    end
end
for i = n+1:2*n
	plot3(newv(i:2*n:2*n+i,1),newv(i:2*n:2*n+i,2),newv(i:2*n:2*n+i,3),'b--','LineWidth',2); 
end
hold off;
title('signed directed');
view(15,90)
[az,el] = view;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
hold on;
axis off;
for i = 1:n
	text(newv(i,1),newv(i,2),newv(i,3),['  ' name{i,1} '_{out}']);
 	text(newv(2*n+i,1),newv(2*n+i,2),newv(2*n+i,3),['  _{in}']);
end
plot3(newv(1:n,1),newv(1:n,2),newv(1:n,3),'g^');
plot3(newv(2*n+1:3*n,1),newv(2*n+1:3*n,2),newv(2*n+1:3*n,3),'g*');
for i = 1:n
    for j = 1:n
        if posW(i,j) ~= 0
            plot3(newv(i:j-i+2*n:j+2*n,1),newv(i:j-i+2*n:j+2*n,2),newv(i:j-i+2*n:j+2*n,3),'g-');
        end
    end
end
for i = 1:n
	plot3(newv(i:2*n:2*n+i,1),newv(i:2*n:2*n+i,2),newv(i:2*n:2*n+i,3),'b--','LineWidth',2); 
end
title('signed directed positive part');
view(az,el)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
hold on;
axis off;
for i = 1:n
	text(newv(n+i,1),newv(n+i,2),newv(n+i,3),['  ' name{i,1} '_{out}']);
	text(newv(3*n+i,1),newv(3*n+i,2),newv(3*n+i,3),['  _{in}']);
end

plot3(newv(n+1:2*n,1),newv(n+1:2*n,2),newv(n+1:2*n,3),'r^'); 
plot3(newv(3*n+1:4*n,1),newv(3*n+1:4*n,2),newv(3*n+1:4*n,3),'r*');
for i = 1:n
    for j = 1:n
        if negW(i,j) ~= 0
            plot3(newv(i+n:j-i+2*n:j+3*n,1),newv(i+n:j-i+2*n:j+3*n,2),newv(i+n:j-i+2*n:j+3*n,3),'r-');
        end
    end
end
for i = n+1:2*n
	plot3(newv(i:2*n:2*n+i,1),newv(i:2*n:2*n+i,2),newv(i:2*n:2*n+i,3),'b--','LineWidth',2); 
end
title('signed directed negative part');
view(az,el)
hold off;
