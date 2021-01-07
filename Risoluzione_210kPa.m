clc 
close all 
clear

% Minimizzazione
options = optimset('PlotFcns','optimplotfval','MaxFunEvals',100000000, 'Maxiter',60, 'TolFun',0.1);
[v,fval] = fminsearch(@Funzione_da_minimizzare_Prova210kPa,[1.31305776499711;0.763877996412080;3129.42987241491;328.150022815106;0.121921485236332;1.28258934854589],options)

%Risultato ottenuto
v0 =[1.3129;0.7721;3131.078;331.200;0.1247;1.2510]

%Calcolo della matrice Hessiana, autovalori, autovettori e numero di
%condizionamento
options=optimset('LargeScale','off','GradObj','off','HessUpdate','steepdesc');
options = optimoptions('fmincon','Algorithm','interior-point','SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true,'ConstraintTolerance',0.00000000001,'FiniteDifferenceStepSize',sqrt(eps),'CheckGradients', true,'HessianApproximation','finite-difference','FiniteDifferenceType','central','FunctionTolerance', 0.0000000001,'MaxFunctionEvaluations', 100000, 'MaxIterations', 10000, 'OptimalityTolerance', 0.000000000001, 'FiniteDifferenceStepSize', 0.000000000001);
[x,fval,exitflag,output,lambda,grad,hessian] = fmincon(@Funzione_da_minimizzare_Prova210kPa,v0)
[Autovettori210, Autovalori210] = eig(hessian)
ConditionNumber210 = cond(hessian)

%Da qui grafici

%% f-m
mm=linspace(1.0,1.8,400);
v=[];
v(3)= 3131;
v(2) = 0.77;
v(4) = 331.2;
v(5) = 0.1251;
v(6) = 1.256;
for i=1:400
    v(1)=mm(i);
    f = Funzione_da_minimizzare_Prova210kPa(v)
    ff(i) =f;
end
figure(1)
plot(mm, ff)
% title('m - f')
x = mm;
y = ff;
[yMin, indexAtMin] = min(y);
xMin = x(indexAtMin)
hold on;
plot(xMin, yMin, 'rv', 'LineWidth', 2, 'MarkerSize', 5);
xlabel('m')
ylabel('function')

%% f-n
nn=linspace(0.5,0.9,400);
v=[];
v(3)= 3131;
v(1) = 1.313;
v(4) = 331.2;
v(5) = 0.1251;
v(6) = 1.256;

for i=1:400
    v(2)=nn(i);
    f = Funzione_da_minimizzare_Prova210kPa(v)
    ff(i) =f;
end
x = nn;
y = ff;
[yMin, indexAtMin] = min(y);
xMin = x(indexAtMin)
figure(2)
plot(nn, ff)
% title('n - f')
hold on;
plot(xMin, yMin, 'rv', 'LineWidth', 2, 'MarkerSize', 5);
xlabel('n')
ylabel('function')

%% f-R1

RR1=linspace(2000,4500,400);
v=[];
v(2)= 0.77;
v(1) = 1.313;
v(4) = 331.2;
v(5) = 0.1251;
v(6) = 1.256;
for i=1:400
    v(3)=RR1(i);
    f = Funzione_da_minimizzare_Prova210kPa(v)
    ff(i) =f;
end
x = RR1;
y = ff;
[yMin, indexAtMin] = min(y);
xMin = x(indexAtMin)
figure(3)
plot(RR1, ff)
% title('R1 - f')
% Plot a red marker there
hold on;
plot(xMin, yMin, 'rv', 'LineWidth', 2, 'MarkerSize', 5);
xlabel('R1')
ylabel('function')

%% f-D1

DD1=linspace(100,500,20);
v=[];
v(3)= 3130;
v(1) = 1.313;
v(2) = 0.77;
v(5) = 0.1251;
v(6) = 1.256;
for i=1:20
    v(4)=DD1(i);
    f = Funzione_da_minimizzare_Prova210kPa(v)
    ff(i) =f;
end
figure(4)
plot(DD1, ff)
% title('D1 - f')
x = DD1;
y = ff;
[yMin, indexAtMin] = min(y);
xMin = x(indexAtMin)
hold on;
plot(xMin, yMin, 'rv', 'LineWidth', 2, 'MarkerSize', 5);
xlabel('D1')
ylabel('function')

%% f-k

kk=linspace(0,0.3,200);
v=[];
v(3)= 3130;
v(1) = 1.313;
v(2) = 0.77;
v(4) = 331.2;
v(6) = 0;
for i=1:200
    v(5)=kk(i);
    f = Funzione_da_minimizzare_Prova210kPa(v)
    ff(i) =f;
end
figure(5)
plot(kk, ff)
% title('k - f')
x = kk;
y = ff;
[yMin, indexAtMin] = min(y);
xMin = x(indexAtMin)
hold on;
plot(xMin, yMin, 'rv', 'LineWidth', 2, 'MarkerSize', 5);
xlabel('k')
ylabel('function')


%% f-cc

cc=linspace(-50,100,200);
v=[];
v(3)= 3130;
v(1) = 1.313;
v(2) = 0.77;
v(4) = 331.2;
v(5) = 0.125;
for i=1:200
    v(6)=cc(i);
    f = Funzione_da_minimizzare_Prova210kPa(v)
    ff(i) =f;
end
figure(6)
plot(cc, ff)
% title('Critical Stress - f')
x = cc;
y = ff;
[yMin, indexAtMin] = min(y);
xMin = x(indexAtMin)
hold on;
plot(xMin, yMin, 'rv', 'LineWidth', 2, 'MarkerSize', 5);
xlabel('Critical Stress')
ylabel('function')
