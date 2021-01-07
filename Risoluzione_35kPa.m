clc 
close all 
clear

options = optimset('PlotFcns','optimplotfval','MaxFunEvals',10000000, 'Maxiter',20, 'TolFun',0.0001);

%Minimizzo la funzione e trovo i 6 parametri 
% [v,fval] = fminsearch(@Funzione_da_minimizzare_Prova35kPa,[1.188099001912;0.610804930181348;2028.00880840107;43.142740488595;0.208315639175806;0.000496686277694060],options)

%Dal vettore che ottengo alla fine della minimizzazione (che chiamo
%v0),Calcolo la matrice hessiana per differenze finite e i rispettivi
%autovalori ed autovettori

v0 =[1.1802;0.617;2024.013;43.4069;0.207;0.0084]
options = optimset('PlotFcns','optimplotfval','MaxFunEvals',100000000, 'Maxiter',60, 'TolFun',0.1);
[v,fval] = fminsearch(@Funzione_da_minimizzare_Prova35kPa,[1.0852;0.6159;2023.9876;43.417;0.20755;0.0005],options)
options = optimoptions('fmincon','Algorithm','interior-point','SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true,'ConstraintTolerance',0.00000000001,'FiniteDifferenceStepSize',sqrt(eps),'CheckGradients', true,'HessianApproximation','finite-difference','FiniteDifferenceType','central','FunctionTolerance', 0.0000000001,'MaxFunctionEvaluations', 100000, 'MaxIterations', 10000, 'OptimalityTolerance', 0.000000000001, 'FiniteDifferenceStepSize', 0.000000000001);
[x,fval,exitflag,output,lambda,grad,hessian] = fmincon(@Funzione_da_minimizzare_Prova35kPa,v0)
[Autovettori35, Autovalori35] = eig(hessian)
ConditionNumber35 = cond(hessian)

%DA QUI GRAFICI

%% f-m
[1.18526590880511;0.615928796610345;2023.98764579424;43.4174955102268;0.207558444406141;0.000572569332164720] 
mm=linspace(1.1,1.7,100);
v=[];
v(3)= 2024;
v(2) = 0.616;
v(4) = 43.4;
v(5) = 0.207;
v(6) = 0;
for i=1:100
    v(1)=mm(i);
    f = Funzione_da_minimizzare_Prova35kPa(v)
    ff(i) =f;
end
figure (1)
plot(mm, ff)
% title('m - f')
x = mm;
y = ff;
[yMin, indexAtMin] = min(y);
xMin = x(indexAtMin)
hold on;
plot(xMin, yMin, 'rv', 'LineWidth', 2, 'MarkerSize', 5);
xlabel('m','FontSize', 14)
ylabel('function','FontSize', 12)

%% f-n
nn=linspace(0.5,0.9,200);
v=[];
v(1) = 1.1852;
v(3)= 2024;
v(4) = 43.4;
v(5) = 0.207;
v(6) = 0;

for i=1:200
    v(2)=nn(i);
    f = Funzione_da_minimizzare_Prova35kPa(v)
    ff(i) =f;
end
x = nn;
y = ff;
[yMin, indexAtMin] = min(y);
xMin = x(indexAtMin)
figure (2)
plot(nn, ff)
% title('n - f')
hold on;
plot(xMin, yMin, 'rv', 'LineWidth', 2, 'MarkerSize', 5);
xlabel('n','FontSize', 14)
ylabel('function','FontSize', 12)

%% f-R1

RR1=linspace(1500,2200,200);
v=[];
v(2)= 0.616;
v(1) = 1.1852;
v(4) = 43.4;
v(5) = 0.207;
v(6) = 0;
for i=1:200
    v(3)=RR1(i);
    f = Funzione_da_minimizzare_Prova35kPa(v)
    ff(i) =f;
end
x = RR1;
y = ff;
[yMin, indexAtMin] = min(y);
xMin = x(indexAtMin)
figure (3)
plot(RR1, ff)
% title('R1 - f')
% Plot a red marker there
hold on;
plot(xMin, yMin, 'rv', 'LineWidth', 2, 'MarkerSize', 5);
xlabel('R1','FontSize', 14)
ylabel('function','FontSize', 12)

%% f-D1

DD1=linspace(10,100,200);
v=[];
v(3)= 2024;
v(2)= 0.616;
v(1) = 1.1852;
v(5) = 0.207;
v(6) = 0;
for i=1:200
    v(4)=DD1(i);
    f = Funzione_da_minimizzare_Prova35kPa(v)
    ff(i) =f;
end
figure (4)
plot(DD1, ff)
% title('D1 - f')
x = DD1;
y = ff;
[yMin, indexAtMin] = min(y);
xMin = x(indexAtMin)
hold on;
plot(xMin, yMin, 'rv', 'LineWidth', 2, 'MarkerSize', 5);
xlabel('D1','FontSize', 14)
ylabel('function','FontSize', 12)

%% f-k

kk=linspace(0.15,0.25,200);
v=[];
v(3)= 2024;
v(2)= 0.616;
v(1) = 1.1852;

v(6) = 0;
v(4) = 43.4;
for i=1:200
    v(5)=kk(i);
    f = Funzione_da_minimizzare_Prova35kPa(v)
    ff(i) =f;
end
figure (5)
plot(kk, ff)
% title('k - f')
x = kk;
y = ff;
[yMin, indexAtMin] = min(y);
xMin = x(indexAtMin)
hold on;
plot(xMin, yMin, 'rv', 'LineWidth', 2, 'MarkerSize', 5);
xlabel('k','FontSize', 14)
ylabel('function','FontSize', 12)


%% f-cc

cc=linspace(-10,20,200);
v=[];
v(3)= 2024;
v(2)= 0.616;
v(1) = 1.1852;
v(5) = 0.207;
v(4) = 43.4;
for i=1:200
    v(6)=cc(i);
    f = Funzione_da_minimizzare_Prova35kPa(v)
    ff(i) =f;
end
figure (6)
plot(cc, ff)
% title('Critical Stress - f')
x = cc;
y = ff;
[yMin, indexAtMin] = min(y);
xMin = x(indexAtMin)
hold on;
plot(xMin, yMin, 'rv', 'LineWidth', 2, 'MarkerSize', 5);
xlabel('Critical Stress','FontSize', 14)
ylabel('function','FontSize', 12)




