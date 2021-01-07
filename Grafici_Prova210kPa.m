clc
clear
close all

%% Grafico n-m

mm=linspace(1.2,1.37,20);
nn=linspace(0.5,0.9,20);

v=[];
v(3) = 3130;
v(4) = 331.2;
v(5) = 0.124;
v(6) = 1.250;

ff=zeros(20,20);
for i=1:20
    v(1)=mm(i);
    for j=1:20
        v(2)=nn(j);
        f =Funzione_da_minimizzare_Prova210kPa(v)
        ff(i,j)=f;
    end
end

[x,y]= meshgrid(mm,nn);
[~, minIdx] = min(ff(:)); 
[row,col] = ind2sub(size(ff),minIdx); 
xMin = x(row,col); 
yMin = y(row,col); 
contour(mm, nn, ff,[0:0.2:10])
colorbar
% hold on
% plot(xMin, yMin, 'rx','MarkerSize', 5)
% title('m - n')
xlabel('m','FontSize', 14)
ylabel('n','FontSize', 14)

figure (1)
contour(mm, nn, ff,[50:5:150])
colorbar
hold on
plot(xMin, yMin, 'rx')


%% Grafico m-R1 - SI
mm=linspace(1.05,1.6,20);
RR1=linspace(2800,3600,20);

v=[];
v(2) = 0.77;
v(4) = 331.2;
v(5) = 0.1251;
v(6) = 1.256;

ff=zeros(20,20);
for i=1:20
    v(1)=mm(i);
    for j=1:20
        v(3)=RR1(j);
        f = Funzione_da_minimizzare_Prova210kPa(v)
        ff(i,j)=f;
    end
end

[x,y]= meshgrid(mm,RR1);
[~, minIdx] = min(ff(:)); 
[row,col] = ind2sub(size(ff),minIdx); 
xMin = x(row,col); 
yMin = y(row,col); 
% 

figure (2)
contour(mm, RR1, ff,[0:0.1:5])
colorbar
% hold on
% plot(xMin, yMin, 'rx','MarkerSize', 10)
% title('m - R1')
xlabel('m','FontSize', 14)
ylabel('R1','FontSize', 14)


%% Grafico m-D1

mm=linspace(1.05,1.6,20);
DD1=linspace(250,400,20);

v=[];
v(3)= 3131;
v(2) = 0.77;
v(5) = 0.1251;
v(6) = 1.256;

ff=zeros(20,20);
for i=1:20
    v(1)=mm(i);
    for j=1:20
        v(4)=DD1(j);
        f = Funzione_da_minimizzare_Prova210kPa(v)
        ff(i,j)=f;
    end
end

[x,y]= meshgrid(mm,DD1);
[~, minIdx] = min(ff(:)); 
[row,col] = ind2sub(size(ff),minIdx); 
xMin = x(row,col); 
yMin = y(row,col); 
contourf(mm, DD1, ff,[0:2:30])
colorbar
% hold on
% plot(xMin, yMin, 'rx', 'MarkerSize', 10)
% title('m - D1')
xlabel('m','FontSize', 14)
ylabel('D1','FontSize', 14) 

%% Grafico m-k

mm=linspace(1.1,1.6,20);
kk=linspace(0.1,0.16,20);

v=[];
v(3)= 3130;
v(2) = 0.77;
v(4) = 331.2;
v(6) = 1.256;

ff=zeros(20,20);
for i=1:20
    v(1)=mm(i);
    for j=1:20
        v(5)=kk(j);
        f = Funzione_da_minimizzare_Prova210kPa(v)
        ff(i,j)=f;
    end
end

[x,y]= meshgrid(mm,kk);
[~, minIdx] = min(ff(:)); 
[row,col] = ind2sub(size(ff),minIdx); 
xMin = x(row,col); 
yMin = y(row,col); 
contourf(mm, kk, ff,[0:2:20])
colorbar
% hold on
% plot(xMin, yMin, 'rx', 'MarkerSize', 10)
% title('m - k')
xlabel('m','FontSize', 14)
ylabel('k','FontSize', 14) 


%% Grafico m-Tau

mm=linspace(1.15,1.6,20);
CCriticalStress=linspace(-15,30,20);

v=[];
v(3)= 3130;
v(2) = 0.77;
v(4) = 331.2;
v(5) = 0.1251;

ff=zeros(20,20);
for i=1:20
    v(1)=mm(i);
    for j=1:20
        v(6)=CCriticalStress(j);
        f = Funzione_da_minimizzare_Prova210kPa(v)
        ff(i,j)=f;
    end
end

[x,y]= meshgrid(mm,CCriticalStress);
[~, minIdx] = min(ff(:)); 
[row,col] = ind2sub(size(ff),minIdx); 
xMin = x(row,col); 
yMin = y(row,col); 
contourf(mm, CCriticalStress, ff,[0:0.3:20])
colorbar
% hold on
% plot(xMin, yMin, 'rx')
% title('m - CriticalStress')
xlabel('m','FontSize', 14)
ylabel('Critical Stress','FontSize', 14) 


%% Grafico n-R1

nn=linspace(0.5,0.9,20);
RR1=linspace(3100,3150,20);

v=[];
v(4)= 331.2;
v(1) = 1.313;
v(5) = 0.125;
v(6) = 1.25;

ff=zeros(20,20);
for i=1:20
    v(2)=nn(i);
    for j=1:20
        v(3)=RR1(j);
        f = Funzione_da_minimizzare_Prova210kPa(v)
        ff(i,j)=f;
    end
end

[x,y]= meshgrid(nn,RR1);
% [~, minIdx] = min(ff(:)); 
% [row,col] = ind2sub(size(ff),minIdx); 
% xMin = x(row,col); 
% yMin = y(row,col); 
contour(nn, RR1, ff,[0:0.22:20])
colorbar
% hold on
% plot(xMin, yMin, 'rx')
% title('n - R1')
xlabel('n','FontSize', 14)
ylabel('R1','FontSize', 14) 


