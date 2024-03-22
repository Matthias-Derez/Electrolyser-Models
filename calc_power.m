function [Power] = calc_power(jrange,Trange)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here


% Trange = Tmin:1:Tmax; %K
% jrange = jmin:(jmax-jmin)/(Tmax-Tmin):jmax; %A/m²

R = 8.314; %J/(mol K)
F = 96485.3; %sA/mol
P = 1; %bar
pH2 = 0.5*P; %bar
pH2O = 0.5*P; %bar
A = 0.0016; %m² (22)
pO2 = 1*P; %Simplification, not taken into account: air as sweep gas (79% N2) and O2 production
gamma_c = 1.344*10^10; %A/m² (84)
gamma_a = 2.051*10^9; %A/m²(84)
E_actc = 10^5; %J/mol(84)
E_acta = 1.2*10^5; %J/mol(84)
d_c = 50*10^(-5); %m waarde van (84)
d_a = 50*10^(-6); %m waarde van (84)
d_e = 50*10^(-6); %m waarde van (84)
sigma_a = 8.4*10^3; %1/(Ohm*m) (64)
sigma_c = 80*10^3; %1/(Ohm*m) (64)
eps_H2k = 59.7; %K (84)
eps_H2Ok = 809.1; %K (84)
sigma_H2O = 2.641; %Angstrom (84)
sigma_H2 = 2.827; %Angstrom (84)
sigma_H2OH2 = (sigma_H2+sigma_H2O)/2; %Angstrom (84)
M_H2O = 18.01528; %kg/kmol
M_H2 = 2.016; %kg/kmol
r_p = 5*10^(-7); %m (84)
chi = 5; %- (84)
eps = 0.4; %- (84) (0.3 in 22 en 64)    
U_total = zeros(length(jrange), length(Trange));
U_rev = zeros(length(jrange), length(Trange));
U_act = zeros(length(jrange), length(Trange));
U_ohm = zeros(length(jrange), length(Trange));
U_conc = zeros(length(jrange), length(Trange));
Power = zeros(length(jrange), length(Trange));
index_current = 0;
for j = jrange
    index_current = index_current + 1;
    index_temp = 0;
    for T = Trange
        index_temp = index_temp + 1;
        j_0a = gamma_a*exp(-E_acta/(R*T)); %A/m² (84)
        j_0c = gamma_c*exp(-E_actc/(R*T)); %A/m² (84)        
        sigma_e = 33.4*10^3*exp(-10.3*10^3/T); %1/(Ohm*m) (64)
        %U_ohmtest = 2.99*10^(-5)*j*d_e*exp(10300/T);        
        Tstar = T/sqrt(eps_H2k*eps_H2Ok); %- (84)
        OmegaD = 1.06036/Tstar^0.1561 + 0.193/exp(0.47635*Tstar) + 1.03587/exp(1.52996*Tstar) + 1.76474/exp(3.89411*Tstar); %dimensionless (84)
        DH2OH2 = 0.00188*sqrt(1/M_H2 + 1/M_H2O)*(T^(3/2))./(P*10^(-5)*sigma_H2OH2^2*OmegaD); %cm²/s (84)
        DH2Ok = 9700*r_p*sqrt(T/M_H2O)*10^2; %cm²/s (84)
        D_effH2O = (chi/eps*(1/DH2OH2+1/DH2Ok))^(-1); %cm²/s (84) 
        tau = T/1000; %(84)
        mu = -1.6918 + 889.75*tau - 892.79*tau^2 + 905.98*tau^3 - 598.36*tau^4 + 221.64*tau^5 -34.754*tau^6; %10^(-7) Pa s (84)
        B_g = eps^3/(72*chi*(1-eps)^2)*(2*r_p)^2; %m² (84)
        U_rev(index_current,index_temp) = 1.253-2.4516*10^(-4)*T+R*T/(2*F)*log(pH2*(pO2)^(1/2)/pH2O);
        U_actc = R*T/F*(asinh(j/(2*j_0c)));
        U_acta = R*T/F*(asinh(j/(2*j_0a)));
        U_act(index_current,index_temp) = U_acta + U_actc;
        U_ohm(index_current,index_temp) = j*(d_c/sigma_c+d_e/sigma_e+d_a/sigma_a);
        U_concc = R*T/(2*F)*log((1+(j*R*T*d_c)/(2*F*(pH2*10^5)*(D_effH2O*10^(-4))))/(1-(j*R*T*d_c)/(2*F*(pH2O*10^5)*(D_effH2O*10^(-4)))));
        U_conca = R*T/(4*F)*log(sqrt((pO2*10^5)^2+(j*R*T*mu*d_a/(2*F*B_g)))/(pO2*10^5));
        U_conc(index_current,index_temp) = U_conca + U_concc;
        U_total(index_current,index_temp) = U_rev(index_current,index_temp) + U_act(index_current,index_temp) + U_ohm(index_current,index_temp) + U_conc(index_current,index_temp);
        Power(index_current,index_temp) = j*A*U_total(index_current,index_temp);
        
    end
end

end