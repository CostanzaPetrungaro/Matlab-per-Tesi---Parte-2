
clc
clear
close all


%% Grafico m-R1 
mm=linspace(1.1,1.5,30);
RR1=linspace(2010,2035,30);
% RR1=linspace(2010,2040,20);
v=[];
v(6) = 0;
v(2)= 0.616;
v(5) = 0.207;
v(4) = 43.4;

ff=zeros(30,30);
for i=1:30
    v(1)=mm(i);
    for j=1:30
        v(3)=RR1(j);
        f = Funzione_da_minimizzare_Prova35kPa(v)
        ff(i,j)=f;
    end
end

[x,y]= meshgrid(mm,RR1);
[~, minIdx] = min(ff(:)); 
[row,col] = ind2sub(size(ff),minIdx); 
xMin = x(row,col); 
yMin = y(row,col); 

figure (1)
contour(mm, RR1, ff,[0:0.5:5])
colorbar
% hold on
% plot(xMin, yMin, 'rx','MarkerSize', 10)
% title('m - R1')
xlabel('m','FontSize', 14)
ylabel('R1','FontSize', 14)


%% Grafico m-D1

mm=linspace(1.1,1.35,20);
DD1=linspace(38,50,20);

v=[];
v(6) = 0;
v(3)= 2024;
v(2)= 0.616;

v(5) = 0.207;


ff=zeros(20,20);
for i=1:20
    v(1)=mm(i);
    for j=1:20
        v(4)=DD1(j);
        f = Funzione_da_minimizzare_Prova35kPa(v)
        ff(i,j)=f;
    end
end

[x,y]= meshgrid(mm,DD1);
[~, minIdx] = min(ff(:)); 
[row,col] = ind2sub(size(ff),minIdx); 
xMin = x(row,col); 
yMin = y(row,col); 
contour(mm, DD1, ff,[0:0.25:5])
colorbar
% hold on
% plot(xMin, yMin, 'rx', 'MarkerSize', 10)
% title('m - D1')
xlabel('m','FontSize', 14)
ylabel('D1','FontSize', 14) 

%% Grafico m-k

mm=linspace(1.0,1.35,20);
kk=linspace(0.19,0.22,20);

v=[];
v(6) = 0;
v(3)= 2024;
v(2)= 0.616;

v(4) = 43.4;
ff=zeros(20,20);
for i=1:20
    v(1)=mm(i);
    for j=1:20
        v(5)=kk(j);
        f = Funzione35(v)
        ff(i,j)=f;
    end
end

[x,y]= meshgrid(mm,kk);
[~, minIdx] = min(ff(:)); 
[row,col] = ind2sub(size(ff),minIdx); 
xMin = x(row,col); 
yMin = y(row,col); 
contour(mm, kk, ff,[0:0.2:5])
colorbar
% hold on
% plot(xMin, yMin, 'rx', 'MarkerSize', 10)
% title('m - k')
xlabel('m','FontSize', 14)
ylabel('k','FontSize', 14) 


%% Grafico m-Tau
% 
mm=linspace(1.0,1.5,20);
CCriticalStress=linspace(-2,3,20);

v=[];
v(3)= 2024;
v(2)= 0.616;
v(5) = 0.207;
v(4) = 43.4;

ff=zeros(20,20);
for i=1:20
    v(1)=mm(i);
    for j=1:20
        v(6)=CCriticalStress(j);
        f = Funzione35(v)
        ff(i,j)=f;
    end
end

[x,y]= meshgrid(mm,CCriticalStress);
[~, minIdx] = min(ff(:)); 
[row,col] = ind2sub(size(ff),minIdx); 
xMin = x(row,col); 
yMin = y(row,col); 
contour(mm, CCriticalStress, ff,[0:0.25:5])
colorbar
% hold on
% plot(xMin, yMin, 'rx')
% title('m - CriticalStress')
xlabel('m','FontSize', 14)
ylabel('Critical Stress','FontSize', 14) 


%% Grafico n-R1

nn=linspace(0.4,0.8,20);
RR1=linspace(2010,2030,20);

v=[];
v(6) = 0;
v(1) = 1.1852;
v(5) = 0.207;
v(4) = 43.4;

ff=zeros(20,20);
for i=1:20
    v(2)=nn(i);
    for j=1:20
        v(3)=RR1(j);
        f = Funzione35(v)
        ff(i,j)=f;
    end
end

[x,y]= meshgrid(nn,RR1);
% [~, minIdx] = min(ff(:)); 
% [row,col] = ind2sub(size(ff),minIdx); 
% xMin = x(row,col); 
% yMin = y(row,col); 
contour(nn, RR1, ff,[0:0.2:10])
% contour(nn, RR1, ff)
colorbar
% hold on
% plot(xMin, yMin, 'rx')
% title('n - R1')
xlabel('n','FontSize', 14)
ylabel('R1','FontSize', 14) 


%% Grafico n-D1
