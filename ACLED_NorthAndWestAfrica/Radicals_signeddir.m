clc
clear all
close all
path(path, '..\new_Laplacian');
 %chosenrows_Radicals 173 has 173 nodes. There are 111 nodes with at least one positive edge.
chosenrows = csvread('chosenrows_Radicals173.csv');
chosennames = csvread('chosennames_Radicals173.csv');
% chosenrows = csvread('chosenrows_Radicals116.csv');
% chosennames = csvread('chosennames_Radicals116.csv');
Apostemp = csvread('posdirected_Radicals.csv');
Anegtemp = csvread('negdirected_Radicals.csv');

Apos = Apostemp(chosenrows,chosenrows);
Aneg = Anegtemp(chosenrows,chosenrows);
Apos=Apos-diag(diag(Apos));
Aneg=Aneg-diag(diag(Aneg));

Dpos=sum(Apos + Apos');

%%%%%%%%%%% only keep nodes with at least one positive edge.
%%%%%%%%%%% This step will only keeps 111 nodes in the analysis for cpt=10;
% Apos=Apos(Dpos>0,Dpos>0);
% Aneg=Aneg(Dpos>0,Dpos>0);
% chosennames=chosennames(Dpos>0);
% chosenrows=chosenrows(Dpos>0);

[n,m] = size(Apos);
[newv,lambda] = SignedDirLaplacian(Apos,Aneg,'sns');

figure;
hold on;
axis off;
for i = 1:n
	text(newv(i,1),newv(i,2),newv(i,3),['  ' int2str(chosennames(i)) '_{out}']);
	text(newv(2*n+i,1),newv(2*n+i,2),newv(2*n+i,3),['  ' int2str(chosennames(i)) '_{in}']);
end
plot3(newv(1:n,1),newv(1:n,2),newv(1:n,3),'ko');
plot3(newv(2*n+1:3*n,1),newv(2*n+1:3*n,2),newv(2*n+1:3*n,3),'k*');

for i = 1:n
    for j = 1:n
        if Apos(i,j) ~= 0
            plot3(newv(i:j-i+2*n:j+2*n,1),newv(i:j-i+2*n:j+2*n,2),newv(i:j-i+2*n:j+2*n,3),'g-');
        end
    end
end
for i = 1:n
	plot3(newv(i:2*n:2*n+i,1),newv(i:2*n:2*n+i,2),newv(i:2*n:2*n+i,3),'k:','LineWidth',2); 
end
%view(-40,-12) % 16 nodes
view(-47,26) % 173 nodes

figure;
hold on;
axis off;
for i = 1:n
	text(newv(n+i,1),newv(n+i,2),newv(n+i,3),['  ' int2str(chosennames(i)) '_{out}']);
	text(newv(3*n+i,1),newv(3*n+i,2),newv(3*n+i,3),['  ' int2str(chosennames(i)) '_{in}']);
end

plot3(newv(n+1:2*n,1),newv(n+1:2*n,2),newv(n+1:2*n,3),'ko');
plot3(newv(3*n+1:4*n,1),newv(3*n+1:4*n,2),newv(3*n+1:4*n,3),'k*');
for i = 1:n
    for j = 1:n
        if Aneg(i,j) ~= 0
            plot3(newv(i+n:j-i+2*n:j+3*n,1),newv(i+n:j-i+2*n:j+3*n,2),newv(i+n:j-i+2*n:j+3*n,3),'r--');
        end
    end
end
for i = n+1:2*n
	plot3(newv(i:2*n:2*n+i,1),newv(i:2*n:2*n+i,2),newv(i:2*n:2*n+i,3),'k:','LineWidth',2); 
end
hold off;
%view(-40,-12) % 16 nodes
view(-47,26) % 173 nodes

ds = newv(:,1:4);
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
D=diag(sum(Apos+Aneg+Apos'+Aneg',2));
dist2=D*dist; %normalized distance
figure;
imagesc(dist2(:,1:4));
title('Pos-pos; Neg-neg, Pos-neg; Neg-pos');
