function [f ] = Funzione210(v)

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
EpsilonVSperimentale = [-0.00000
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

PoissonRatio = 0.287;   %Ricavato
E0 = 130477;            %Ricavato
p0 = -210;              %Noto

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
MediaEpsilonV = abs(sum(abs(EpsilonVSperimentale))/length(EpsilonVSperimentale));

% f = sum((((EpsilonVSperimentale -(((((EpsilonZSperimentale_V*E0 +(1/v(5) - 1)*E0* (polyval((polyfit(alpha, (sqrt((6*v(5)^2*v(4))./(-(- v(3).*((v(1) - v(2)).*alpha+ v(2)).*((1 - alpha).^(v(1) - 1))./((alpha).^(v(2) + 1)))))),20)),polyval((polyfit(((-3*v(6)/(1-v(5)) - p0*(1+2*v(5))/(1-v(5)) - sqrt(6*v(4)./(1/v(3).* ((v(1) - v(2)).* alpha + v(2))./(alpha.^(1-v(2)).*(1 - alpha).^(v(1) + 1))))/(1-v(5)))/E0 - (1/v(5) - 1)*(sqrt((6*v(5)^2*v(4))./(-(-v(3).*((v(1) - v(2)).*alpha + v(2)).*((1 - alpha).^(v(1) - 1))./((alpha).^(v(2) + 1))))))/3 + 2*PoissonRatio*p0/E0),alpha,20)),EpsilonZSperimentale_V)))/3 +2*PoissonRatio*p0)) - 2*p0)/(3*K0) + polyval((polyfit(alpha, (sqrt((6*v(5)^2*v(4))./(-(- v(3).*((v(1) - v(2)).*alpha+ v(2)).*((1 - alpha).^(v(1) - 1))./((alpha).^(v(2) + 1)))))),20)),polyval((polyfit(((-3*v(6)/(1-v(5)) - p0*(1+2*v(5))/(1-v(5)) - sqrt(6*v(4)./(1/v(3).* ((v(1) - v(2)).* alpha + v(2))./(alpha.^(1-v(2)).*(1 - alpha).^(v(1) + 1))))/(1-v(5)))/E0 - (1/v(5) - 1)*(sqrt((6*v(5)^2*v(4))./(-(-v(3).*((v(1) - v(2)).*alpha + v(2)).*((1 - alpha).^(v(1) - 1))./((alpha).^(v(2) + 1))))))/3 + 2*PoissonRatio*p0/E0),alpha,20)),EpsilonZSperimentale_V)))))/MediaEpsilonV).^2*10 +((SigmaZSperimentale -(EpsilonZSperimentale_Sigma*E0 +(1/v(5) - 1)*E0* (polyval((polyfit(alpha, (sqrt((6*v(5)^2*v(4))./(-(- v(3).*((v(1) - v(2)).*alpha+ v(2)).*((1 - alpha).^(v(1) - 1))./((alpha).^(v(2) + 1)))))),20)),polyval((polyfit(((-3*v(6)/(1-v(5)) - p0*(1+2*v(5))/(1-v(5)) - sqrt(6*v(4)./(1/v(3).* ((v(1) - v(2)).* alpha + v(2))./(alpha.^(1-v(2)).*(1 - alpha).^(v(1) + 1))))/(1-v(5)))/E0 - (1/v(5) - 1)*(sqrt((6*v(5)^2*v(4))./(-(-v(3).*((v(1) - v(2)).*alpha + v(2)).*((1 - alpha).^(v(1) - 1))./((alpha).^(v(2) + 1))))))/3 + 2*PoissonRatio*p0/E0),alpha,20)),EpsilonZSperimentale_Sigma)))/3 +2*PoissonRatio*p0))./SigmaZSperimentale).^2));
f = sum((((EpsilonVSperimentale -(((((EpsilonZSperimentale_V*E0 +(1/v(5) - 1)*E0* (polyval((polyfit(alpha, (sqrt((6*v(5)^2*v(4))./(-(- v(3).*((v(1) - v(2)).*alpha+ v(2)).*((1 - alpha).^(v(1) - 1))./((alpha).^(v(2) + 1)))))),20)),polyval((polyfit(((-3*v(6)/(1-v(5)) - p0*(1+2*v(5))/(1-v(5)) - sqrt(6*v(4)./(1/v(3).* ((v(1) - v(2)).* alpha + v(2))./(alpha.^(1-v(2)).*(1 - alpha).^(v(1) + 1))))/(1-v(5)))/E0 - (1/v(5) - 1)*(sqrt((6*v(5)^2*v(4))./(-(-v(3).*((v(1) - v(2)).*alpha + v(2)).*((1 - alpha).^(v(1) - 1))./((alpha).^(v(2) + 1))))))/3 + 2*PoissonRatio*p0/E0),alpha,20)),EpsilonZSperimentale_V)))/3 +2*PoissonRatio*p0)) - 2*p0)/(3*K0) + polyval((polyfit(alpha, (sqrt((6*v(5)^2*v(4))./(-(- v(3).*((v(1) - v(2)).*alpha+ v(2)).*((1 - alpha).^(v(1) - 1))./((alpha).^(v(2) + 1)))))),20)),polyval((polyfit(((-3*v(6)/(1-v(5)) - p0*(1+2*v(5))/(1-v(5)) - sqrt(6*v(4)./(1/v(3).* ((v(1) - v(2)).* alpha + v(2))./(alpha.^(1-v(2)).*(1 - alpha).^(v(1) + 1))))/(1-v(5)))/E0 - (1/v(5) - 1)*(sqrt((6*v(5)^2*v(4))./(-(-v(3).*((v(1) - v(2)).*alpha + v(2)).*((1 - alpha).^(v(1) - 1))./((alpha).^(v(2) + 1))))))/3 + 2*PoissonRatio*p0/E0),alpha,20)),EpsilonZSperimentale_V)))))/MediaEpsilonV).^2 +((SigmaZSperimentale -(EpsilonZSperimentale_Sigma*E0 +(1/v(5) - 1)*E0* (polyval((polyfit(alpha, (sqrt((6*v(5)^2*v(4))./(-(- v(3).*((v(1) - v(2)).*alpha+ v(2)).*((1 - alpha).^(v(1) - 1))./((alpha).^(v(2) + 1)))))),20)),polyval((polyfit(((-3*v(6)/(1-v(5)) - p0*(1+2*v(5))/(1-v(5)) - sqrt(6*v(4)./(1/v(3).* ((v(1) - v(2)).* alpha + v(2))./(alpha.^(1-v(2)).*(1 - alpha).^(v(1) + 1))))/(1-v(5)))/E0 - (1/v(5) - 1)*(sqrt((6*v(5)^2*v(4))./(-(-v(3).*((v(1) - v(2)).*alpha + v(2)).*((1 - alpha).^(v(1) - 1))./((alpha).^(v(2) + 1))))))/3 + 2*PoissonRatio*p0/E0),alpha,20)),EpsilonZSperimentale_Sigma)))/3 +2*PoissonRatio*p0))./SigmaZSperimentale).^2));

end 

