function [f] = Funzione35(v)

SigmaZSperimentale = [-35.00
-96.85
-131.25
-165.62
-213.74
-268.71
-296.23
-337.48
-364.96
-396.39
-443.50
-480.15
-516.81
-548.22
-569.10
-595.30
-610.97
-637.13
-647.56
-668.45
-663.13
-668.30
-657.74
-647.16
-641.83
-631.21
-625.80
-625.70
-615.06
-599.06
-572.59
-546.09
]';

EpsilonZSperimentale_Sigma = [0.0000
-0.0016
-0.0025
-0.0034
-0.0047
-0.0061
-0.0069
-0.0080
-0.0087
-0.0096
-0.0112
-0.0125
-0.0136
-0.0146
-0.0161
-0.0170
-0.0181
-0.0193
-0.0206
-0.0222
-0.0237
-0.0251
-0.0264
-0.0282
-0.0298
-0.0322
-0.0354
-0.0373
-0.0401
-0.0452
-0.0499
-0.0552
]';

%Valori Sperimentali Curva SigmaZ - EpsilonZ
EpsilonVSperimentale = [0
-0.000189
-0.0003809
-0.0005726
-0.0007219
-0.0008284
-0.000828
-0.0008277
-0.0006564
-0.0005278
-0.0003565
-0.0001853
-0.0000142
0.0002212
0.0004137
0.000649
0.0007987
0.0010554
0.0013334
0.0016328
0.0018894
0.0021246
0.0024026
0.0026165
0.0029373
0.0032366
0.0034718
0.0037497
0.0040064
0.004263
0.0044768
0.0046907
 ]';
% 

EpsilonZSperimentale_V = [0.0000
-0.0014
-0.0026
-0.0040
-0.0050
-0.0060
-0.0070
-0.0079
-0.0089
-0.0098
-0.0109
-0.0117
-0.0122
-0.0130
-0.0138
-0.0144
-0.0147
-0.0155
-0.0162
-0.0168
-0.0174
-0.0179
-0.0185
-0.0191
-0.0198
-0.0204
-0.0208
-0.0212
-0.0219
-0.0224
-0.0228
-0.0234
]';
 
%% DATI 
%Dati
PoissonRatio = 0.426;  %Ricavato
E0 = 38014;            %Ricavato
p0 = -35;              %Noto

m = v(1);
n = v(2);
R1 = v(3);
D1 = v(4);
k = v(5);
CriticalStress = v(6);

K0 = 1/3*E0/(1-2*PoissonRatio);
alpha = 0:0.01:1;

% CALCOLI - FUNZIONI 
for i = 1:length(alpha)
    R(i) = R1 *((1-alpha(i))^ m)/(alpha(i)^n); 
    RFirstDerivative(i) = - R1*((m - n)*alpha(i) + n)*((1 - alpha(i))^(m - 1))/((alpha(i))^(n + 1));
    RSecondDerivative(i)  = R1*((m-n)*(m-n-1)*alpha(i)^2 + 2*n*(m - n - 1)*alpha(i) + n*(n+1))*((1-alpha(i))^(m-2))/(alpha(i)^(n+2));
    S(i)= 1/R1* (alpha(i) ^ n)/((1 - alpha(i))^m);
    SFirstDerivative(i) = 1/R1 * ((m - n)* alpha(i) + n)/(alpha(i)^(1-n)*(1 - alpha(i))^(m + 1));
    SSecondDerivative(i) = 1/R1*((m-n)*(m-n+1)*alpha(i)^2 + 2*n*(m - n + 1)*alpha(i) - n*(1-n))/(alpha(i)^(2-n)*(1-alpha(i))^(m+2));
    
    Fi = SSecondDerivative./RSecondDerivative.*((-RFirstDerivative).^(3/2))./((SFirstDerivative).^(3/2));
    FiApprossimata(i) = (m+1)/(m-1)*(1-alpha(i))^m;
    Fi0(i) = -(1+n)/(1-n)*alpha(i)^(-n);
    Trp(i)= sqrt((6*k^2*D1)/(-RFirstDerivative(i)));    
    SigmaZ(i) = -3*CriticalStress/(1-k) - p0*(1+2*k)/(1-k) - sqrt(6*D1/(SFirstDerivative(i)))/(1-k);    
    EpsilonZ (i) = SigmaZ(i)/E0 - (1/k - 1)* Trp(i)/3 + 2*PoissonRatio*p0/E0;
    EpsilonV (i) = (SigmaZ(i) - 2*p0)/(3*K0) + Trp (i); 
   
end

EpsilonV(101)=[];
EpsilonZ(101)=[];
alpha(101) =[];
SigmaZ(101)  =[];
Trp(101)  =[];
SFirstDerivative(101) = [];
RFirstDerivative(101) = [];
SSecondDerivative(101) = [];
RSecondDerivative(101) = [];
S(101) = [];
R(101) = [];

MediaSigmaZ = abs((sum(SigmaZSperimentale))/length(SigmaZSperimentale));
MediaEpsilonV = abs(sum(EpsilonVSperimentale)/length(EpsilonVSperimentale));

f = sum((((EpsilonVSperimentale -(((((EpsilonZSperimentale_V*E0 +(1/v(5) - 1)*E0* (polyval((polyfit(alpha, (sqrt((6*v(5)^2*v(4))./(-(- v(3).*((v(1) - v(2)).*alpha+ v(2)).*((1 - alpha).^(v(1) - 1))./((alpha).^(v(2) + 1)))))),20)),polyval((polyfit(((-3*v(6)/(1-v(5)) - p0*(1+2*v(5))/(1-v(5)) - sqrt(6*v(4)./(1/v(3).* ((v(1) - v(2)).* alpha + v(2))./(alpha.^(1-v(2)).*(1 - alpha).^(v(1) + 1))))/(1-v(5)))/E0 - (1/v(5) - 1)*(sqrt((6*v(5)^2*v(4))./(-(-v(3).*((v(1) - v(2)).*alpha + v(2)).*((1 - alpha).^(v(1) - 1))./((alpha).^(v(2) + 1))))))/3 + 2*PoissonRatio*p0/E0),alpha,20)),EpsilonZSperimentale_V)))/3 +2*PoissonRatio*p0)) - 2*p0)/(3*K0) + polyval((polyfit(alpha, (sqrt((6*v(5)^2*v(4))./(-(- v(3).*((v(1) - v(2)).*alpha+ v(2)).*((1 - alpha).^(v(1) - 1))./((alpha).^(v(2) + 1)))))),20)),polyval((polyfit(((-3*v(6)/(1-v(5)) - p0*(1+2*v(5))/(1-v(5)) - sqrt(6*v(4)./(1/v(3).* ((v(1) - v(2)).* alpha + v(2))./(alpha.^(1-v(2)).*(1 - alpha).^(v(1) + 1))))/(1-v(5)))/E0 - (1/v(5) - 1)*(sqrt((6*v(5)^2*v(4))./(-(-v(3).*((v(1) - v(2)).*alpha + v(2)).*((1 - alpha).^(v(1) - 1))./((alpha).^(v(2) + 1))))))/3 + 2*PoissonRatio*p0/E0),alpha,20)),EpsilonZSperimentale_V)))))/MediaEpsilonV).^2 +((SigmaZSperimentale -(EpsilonZSperimentale_Sigma*E0 +(1/v(5) - 1)*E0* (polyval((polyfit(alpha, (sqrt((6*v(5)^2*v(4))./(-(- v(3).*((v(1) - v(2)).*alpha+ v(2)).*((1 - alpha).^(v(1) - 1))./((alpha).^(v(2) + 1)))))),20)),polyval((polyfit(((-3*v(6)/(1-v(5)) - p0*(1+2*v(5))/(1-v(5)) - sqrt(6*v(4)./(1/v(3).* ((v(1) - v(2)).* alpha + v(2))./(alpha.^(1-v(2)).*(1 - alpha).^(v(1) + 1))))/(1-v(5)))/E0 - (1/v(5) - 1)*(sqrt((6*v(5)^2*v(4))./(-(-v(3).*((v(1) - v(2)).*alpha + v(2)).*((1 - alpha).^(v(1) - 1))./((alpha).^(v(2) + 1))))))/3 + 2*PoissonRatio*p0/E0),alpha,20)),EpsilonZSperimentale_Sigma)))/3 +2*PoissonRatio*p0))./SigmaZSperimentale).^2));

end 

