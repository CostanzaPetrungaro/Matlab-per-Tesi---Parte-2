clc
clear
close all



%% VALORI SPERIMENTALI
SigmaZSperimentale = [-210.00
-403.99
-587.47
-755.22
-849.56
-980.58
-1085.40
-1200.70
-1315.98
-1415.54
-1488.89
-1583.20
-1703.66
-1777.01
-1855.46
-1892.06
-1928.69
-1933.84
-1965.13
-1959.74
-1964.86
-1954.19
-1948.79
-1953.88
-1937.88
-1927.13
-1895.43
-1874.19
-1837.17
-1805.45
-1763.23
-1736.75
]';

EpsilonZSperimentale_Sigma = [0.0000
-0.0013
-0.0029
-0.0045
-0.0058
-0.0076
-0.0090
-0.0107
-0.0127
-0.0143
-0.0159
-0.0177
-0.0208
-0.0224
-0.0266
-0.0288
-0.0304
-0.0322
-0.0354
-0.0383
-0.0407
-0.0441
-0.0470
-0.0501
-0.0552
-0.0600
-0.0646
-0.0694
-0.0752
-0.0801
-0.0852
-0.0901
]';

%Valori Sperimentali Curva SigmaZ - EpsilonZ
EpsilonVSperimentale = [-0.00000001
-0.000947439
-0.0013849
-0.0018119
-0.0022385
-0.0025371
-0.0028357
-0.003006
-0.0030907
-0.0029191
-0.0027475
-0.0025334
-0.0023191
-0.0020621
-0.0017198
-0.0013776
-0.0010354
-0.0004792
-0.0000944
0.0001197
0.0004619
0.0009322
0.0011889
0.0014455
0.0017448
0.0020871
0.0023864
0.002643
0.0030278
0.0033274
0.0036694
0.0040114
 ]';
% 

EpsilonZSperimentale_V = [0.0000
-0.0022
-0.0032
-0.0043
-0.0062
-0.0079
-0.0092
-0.0111
-0.0133
-0.0152
-0.0169
-0.0181
-0.0200
-0.0215
-0.0227
-0.0236
-0.0246
-0.0265
-0.0270
-0.0280
-0.0291
-0.0297
-0.0304
-0.0309
-0.0315
-0.0325
-0.0330
-0.0337
-0.0342
-0.0352
-0.0356
-0.0361
]';


%% DATI 

PoissonRatio = 0.287; %Ricavato
E0 = 130477;            %Ricavato
p0 = -210;              %Noto

m = 1.29;
n = 0.77;
R1 =3131;
D1 = 331.2;
k = 0.1251;
CriticalStress =1.25;

K0 = 1/3*E0/(1-2*PoissonRatio);
alpha = 0:0.01:0.999999;
R = zeros(size(alpha));
S = zeros(size(alpha));

%% CALCOLI - FUNZIONI 

for i = 1:length(alpha)
    R(i) = R1 *((1-alpha(i))^ m)/(alpha(i)^n); 
    RFirstDerivative(i) = - R1*((m - n)*alpha(i) + n)*((1 - alpha(i))^(m - 1))/((alpha(i))^(n + 1));
    RSecondDerivative(i)  = R1*((m-n)*(m-n-1)*alpha(i)^2 + 2*n*(m - n - 1)*alpha(i) + n*(n+1))*((1-alpha(i))^(m-2))/(alpha(i)^(n+2));
    S(i)= 1/R1* (alpha(i) ^ n)/((1 - alpha(i))^m);
    SFirstDerivative(i) = 1/R1 * ((m - n)* alpha(i) + n)/(alpha(i)^(1-n)*(1 - alpha(i))^(m + 1));
    SSecondDerivative(i) = 1/R1*((m-n)*(m-n+1)*alpha(i)^2 + 2*n*(m - n + 1)*alpha(i) - n*(1-n))/(alpha(i)^(2-n)*(1-alpha(i))^(m+2));
    Fi = SSecondDerivative./RSecondDerivative.*((-RFirstDerivative).^(3/2))./((SFirstDerivative).^(3/2));    
    Trp(i)= sqrt((6*k^2*D1)/(-RFirstDerivative(i)));    
    SigmaZ(i) = -3*CriticalStress/(1-k) - p0*(1+2*k)/(1-k) - sqrt(6*D1/(SFirstDerivative(i)))/(1-k);    
    EpsilonZ (i) = SigmaZ(i)/E0 - (1/k - 1)* Trp(i)/3 + 2*PoissonRatio*p0/E0;    
    EpsilonV (i) = (SigmaZ(i) - 2*p0)/(3*K0) + Trp (i);   
   end

%% INTERPOLAZIONE Alpha = f(epsilonZ)

p = polyfit(EpsilonZ,alpha,20);
alpha_P = polyval(p,EpsilonZ);

%Grafico
figure(1)
plot(EpsilonZ,alpha,'b', EpsilonZ, alpha_P, 'b')
hold on 
plot(EpsilonZ, alpha_P, '--r','LineWidth',2)
grid on
grid minor
title('\epsilon_z = f(\alpha)', 'FontSize', 14)
ylabel('\alpha', 'FontSize', 12)
xlabel('\epsilon_z',  'FontSize', 12)
figure(1)
plot(EpsilonZ,alpha,'r')
grid on
grid minor
title('\alpha = f(\epsilon_z)', 'FontSize', 14)
ylabel('\alpha', 'FontSize', 14)
xlabel('\epsilon_z',  'FontSize', 14)
xlim([-0.4  0])

%% INTERPOLAZIONE TrP = f(Alpha)= f(f(epsilonZ))

r=polyfit(alpha, Trp,20);
Trp_P =polyval(r,alpha);

figure(2)
plot(Trp_P,alpha,'--r','LineWidth',2)
hold on
plot(Trp, alpha, 'b')
grid on
grid minor
title('Trp = f(\epsilon_z)', 'FontSize', 14)
ylabel('Tr(p)', 'FontSize', 12)
xlabel('\epsilon_z',  'FontSize', 12)

%% Funzione SigmaZ = f(EpsilonZ)
SigmaZ_P = (EpsilonZ*E0 +(1/k - 1)*E0* (polyval(r,polyval(p,EpsilonZ)))/3 -2*PoissonRatio*p0);

%Grafico
figure(3)
plot(EpsilonZ,SigmaZ_P,'--r','LineWidth',2)
hold on
plot(EpsilonZ, SigmaZ,'b','LineWidth',0.5)
% plot(EpsilonZ,SigmaZ_P,'LineWidth',1,'r', EpsilonZ,SigmaZ,'LineWidth',1,'b')
grid on
grid minor
title('\sigma_z = f (\epsilon_z)','FontSize', 14)
ylabel('\sigma_z', 'FontSize', 14)
xlabel('\epsilon_z',  'FontSize', 14)
xlim([-0.25  0])
ylim([-2000  0])
legend({'Funzione Approssimata', 'Funzione Reale'});


%% Funzione EpsilonV = f(EpsilonZ) e Funzione EpsilonV* = f(EpsilonZ*)
EpsilonV_P = (((EpsilonZ*E0 +(1/k - 1)*E0* (polyval(r,polyval(p,EpsilonZ)))/3 -2*PoissonRatio*p0)) - 2*p0)/(3*K0) + polyval(r,polyval(p,EpsilonZ));

%Grafico
figure(4)
plot(EpsilonZ,EpsilonV_P,'--r','LineWidth',2)
hold on
plot(EpsilonZ,EpsilonV,'b','LineWidth',0.5)
grid on
grid minor
title('\epsilon_v = f (\epsilon_z)','FontSize', 14)
ylabel('\epsilon_v', 'FontSize', 12)
xlabel('\epsilon_z',  'FontSize', 12)
xlim([-0.07  0])
legend({'Funzione Approssimata', 'Funzione Reale'});
% Black line x = 0
plot([0 0], ylim, 'k', 'LineWidth', 0.5)
% Black line y = 0
plot(xlim, [0 0], 'k', 'LineWidth', 0.5)

%% Esprimere SigmaZ e EpsilonV in funzione solo dei parametri
%p1 = polyfit(EpsilonZ,alpha,25);
%p1 = polyfit((SigmaZ/E0 - (1/k - 1)* Trp/3 + 2*PoissonRatio*p0/E0),alpha,25);
%p1 = polyfit(((-3*CriticalStress/(1-k) - p0*(1+2*k)/(1-k) - sqrt(6*D1./SFirstDerivative)/(1-k))/E0 - (1/k - 1)*(sqrt((6*k^2*D1)./(-RFirstDerivative)))/3 + 2*PoissonRatio*p0/E0),alpha,25);
%p1 = polyfit(((-3*CriticalStress/(1-k) - p0*(1+2*k)/(1-k) - sqrt(6*D1./(1/R1.* ((m - n).* alpha + n)./(alpha.^(1-n).*(1 - alpha).^(m + 1))))/(1-k))/E0 - (1/k - 1)*(sqrt((6*k^2*D1)./(-RFirstDerivative)))/3 + 2*PoissonRatio*p0/E0),alpha,25);
 p1 = polyfit(((-3*CriticalStress/(1-k) - p0*(1+2*k)/(1-k) - sqrt(6*D1./(1/R1.* ((m - n).* alpha + n)./(alpha.^(1-n).*(1 - alpha).^(m + 1))))/(1-k))/E0 - (1/k - 1)*(sqrt((6*k^2*D1)./(-(-R1.*((m - n).*alpha + n).*((1 - alpha).^(m - 1))./((alpha).^(n + 1))))))/3 + 2*PoissonRatio*p0/E0),alpha,10);

%r1=polyfit(alpha, Trp,25);
%r1=polyfit(alpha, (sqrt((6*k^2*D1)./(-RFirstDerivative))),25);
 r1=polyfit(alpha, (sqrt((6*k^2*D1)./(-(- R1.*((m - n).*alpha+ n).*((1 - alpha).^(m - 1))./((alpha).^(n + 1)))))),10);

%Calcolo il valore medio 
MediaSigmaZ = abs((sum(SigmaZSperimentale))/length(SigmaZSperimentale));
MediaEpsilonV = abs(sum(EpsilonVSperimentale)/length(EpsilonVSperimentale));

v(1)=m;
v(2)=n;
v(3)=R1;
v(4)=D1;
v(5)=k;
v(6)=CriticalStress;

SIGMAZSPE = (EpsilonZSperimentale_Sigma*E0 +(1/v(5) - 1)*E0* (polyval((polyfit(alpha, (sqrt((6*v(5)^2*v(4))./(-(- v(3).*((v(1) - v(2)).*alpha+ v(2)).*((1 - alpha).^(v(1) - 1))./((alpha).^(v(2) + 1)))))),20)),polyval((polyfit(((-3*v(6)/(1-v(5)) - p0*(1+2*v(5))/(1-v(5)) - sqrt(6*v(4)./(1/v(3).* ((v(1) - v(2)).* alpha + v(2))./(alpha.^(1-v(2)).*(1 - alpha).^(v(1) + 1))))/(1-v(5)))/E0 - (1/v(5) - 1)*(sqrt((6*v(5)^2*v(4))./(-(-v(3).*((v(1) - v(2)).*alpha + v(2)).*((1 - alpha).^(v(1) - 1))./((alpha).^(v(2) + 1))))))/3 + 2*PoissonRatio*p0/E0),alpha,20)),EpsilonZSperimentale_Sigma)))/3 +2*PoissonRatio*p0);

%GRAFICO I RISULTATI di SigmaZ
figure(5)
plot(-EpsilonZSperimentale_Sigma, -SigmaZSperimentale,'ro', -EpsilonZSperimentale_Sigma, -SIGMAZSPE, 'b')
grid on
grid minor
xlim([0   0.1])
ylim([0   2250])
title('\sigma_z - \epsilon_Z','FontSize', 15)
xlabel('\epsilon_z', 'FontSize', 14)
ylabel('\sigma_z',  'FontSize', 14) 
hold on
% Black line x = 0
plot([0 0], ylim, 'k', 'LineWidth', 0.5)
% Black line y = 0
plot(xlim, [0 0], 'k', 'LineWidth', 0.5)
legend({'Punti sperimentali', 'Valori stimati'});
 
EPSILONVSPE =(((((EpsilonZSperimentale_V*E0 +(1/v(5) - 1)*E0* (polyval((polyfit(alpha, (sqrt((6*v(5)^2*v(4))./(-(- v(3).*((v(1) - v(2)).*alpha+ v(2)).*((1 - alpha).^(v(1) - 1))./((alpha).^(v(2) + 1)))))),20)),polyval((polyfit(((-3*v(6)/(1-v(5)) - p0*(1+2*v(5))/(1-v(5)) - sqrt(6*v(4)./(1/v(3).* ((v(1) - v(2)).* alpha + v(2))./(alpha.^(1-v(2)).*(1 - alpha).^(v(1) + 1))))/(1-v(5)))/E0 - (1/v(5) - 1)*(sqrt((6*v(5)^2*v(4))./(-(-v(3).*((v(1) - v(2)).*alpha + v(2)).*((1 - alpha).^(v(1) - 1))./((alpha).^(v(2) + 1))))))/3 + 2*PoissonRatio*p0/E0),alpha,20)),EpsilonZSperimentale_V)))/3 +2*PoissonRatio*p0)) - 2*p0)/(3*K0) + polyval((polyfit(alpha, (sqrt((6*v(5)^2*v(4))./(-(- v(3).*((v(1) - v(2)).*alpha+ v(2)).*((1 - alpha).^(v(1) - 1))./((alpha).^(v(2) + 1)))))),20)),polyval((polyfit(((-3*v(6)/(1-v(5)) - p0*(1+2*v(5))/(1-v(5)) - sqrt(6*v(4)./(1/v(3).* ((v(1) - v(2)).* alpha + v(2))./(alpha.^(1-v(2)).*(1 - alpha).^(v(1) + 1))))/(1-v(5)))/E0 - (1/v(5) - 1)*(sqrt((6*v(5)^2*v(4))./(-(-v(3).*((v(1) - v(2)).*alpha + v(2)).*((1 - alpha).^(v(1) - 1))./((alpha).^(v(2) + 1))))))/3 + 2*PoissonRatio*p0/E0),alpha,20)),EpsilonZSperimentale_V))));

%GRAFICO I RISULTATI di EpsilonV 
figure(6)
plot(EpsilonZSperimentale_V, EpsilonVSperimentale,'ro', EpsilonZSperimentale_V, EPSILONVSPE, 'b')
grid on
grid minor
xlim([-0.038  0])
ylim([-0.004  0.004])
title('\epsilon_v - \epsilon_Z','FontSize', 15)
xlabel('\epsilon_z', 'FontSize', 14)
ylabel('\epsilon_v',  'FontSize', 14) 
hold on
% Black line x = 0
plot([0 0], ylim, 'k', 'LineWidth', 0.5)
% Black line y = 0
plot(xlim, [0 0], 'k', 'LineWidth', 0.5)
legend({'Punti sperimentali', 'Valori stimati'});


%% controllo curva epsilonZ-EpsilonV punto iniziale

EPSILONVSPE00 = 0;
EPSILONZSPE00 = -(2*(1+PoissonRatio)*(-2*PoissonRatio+1)*p0)/E0/(-2*PoissonRatio-1);

figure(7)
plot(-EpsilonZSperimentale_V, EpsilonVSperimentale,'ro', -EpsilonZSperimentale_V, EPSILONVSPE, 'b')
grid on
grid minor
xlim([EPSILONZSPE00  0.038])
ylim([-0.004  0.004])
title('\epsilon_v - \epsilon_Z','FontSize', 15)
xlabel('\epsilon_z', 'FontSize', 14)
ylabel('\epsilon_v',  'FontSize', 14) 
hold on
% Black line x = 0
plot([0 0], ylim, 'k', 'LineWidth', 0.5)
% Black line y = 0
plot(xlim, [0 0], 'k', 'LineWidth', 0.5)
legend({'Punti sperimentali', 'Valori stimati'});

