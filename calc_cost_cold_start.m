
%% SOEC
iMinSOEC = 2000 ; %Minimal current density A/m²
iMaxSOEC = 10000; %Maximal current density A/m²
U_tnSOEC = 1.2955; %Thermonteutral voltage V (paper 73)
T_sSOEC = 1273; %Normal working temperature K
C_HSSOEC = 0; %Estimated hot start cost €
C_CSSOEC = 98; %Estimated cold start cost €
T_0SOEC = 1273; %Initial temperature of electrolyser K
C_hSOEC = 141210000; % J/K Electrolyser heat capacity 
R_tSOEC = 0.00325; %K/W Thermal resistance
n_cSOEC = 4707; %Number of cells
TminSOEC = 1073; %K Minimal operating temperature

pi_H = 4; %€/kg
pi_e_average = 4.324e-5; %€/Wh
pi_e_min = 3.702e-5; %€/Wh
m_H2 = 0.0020166669*n_cSOEC*iMaxSOEC*0.21/2/96485.3;
T = 20:1:799;
duration = (800-T)/120; %h duration to heat up electrolyser to min operating temp
Q_CS_e = zeros(length(T),1);
Cost1a_e = zeros(length(T),1);

% Without nuclear heat
for j = 1 :length(T)
    Q_CS_e(j) = C_hSOEC*(800-T(j))/(duration(j)*3600) + ((800+T(j))/2-20)/R_tSOEC;
    Cost1a_e(j) = Q_CS_e(j)*duration(j)/0.95*pi_e_average;
    % dit is de enige cost voor SOEC, geen extra kosten in rekening brengen
    % voor waterstof niet geproduceerd ofzo want 
end
figure
plot(T,Cost1a_e)
title("Cost to heat up SOEC to operating temperature using electricity")
xlabel('Starting temperature [°C]')
ylabel('Cost [€]')

Q_CS_h = zeros(length(T),1);
Cost1a_h = zeros(length(T),1);

% With nuclear heat
for j = 1:length(T)
    Q_CS_h(j) = C_hSOEC*(800-T(j))/(duration(j)*3600) + ((800+T(j))/2-20)/R_tSOEC;
    Cost1a_h(j) = Q_CS_h(j)*duration(j)*0.45*pi_e_average;
    % dit is de enige cost voor SOEC, geen extra kosten in rekening brengen
    % voor waterstof niet geproduceerd ofzo want 
end
figure
plot(T,Cost1a_h)
title("Cost to heat up SOEC to operating temperature with nuclear heat")
xlabel('Starting temperature [°C]')
ylabel('Cost [€]')


% With heat factor 2.33 goedkoper
figure
plot(T,Cost1a_h./Cost1a_e)
title("Cost using heat / cost using electricity")
xlabel('Starting temperature [°C]')
ylabel('Cost [€]')

figure
plot(T,Cost1a_e - Cost1a_h)
title("Cost using electricitiy - cost using heat")
xlabel('Starting temperature [°C]')
ylabel('Cost [€]')

%% AWE without nuclear heat 
iMinAWE = 1000; %Minimal current density A/m²
iMaxAWE = 4000; %Maximal current density A/m²
U_tnAWE = 1.4813; %Thermonteutral voltage V (paper 95)
T_sAWE = 363; %Normal working temperature K
C_HSAWE = 0; %Estimated hot start cost €
C_CSAWE = 98; %Estimated cold start cost €
T_0AWE = 363; %Initial temperature of electrolyser K
C_hAWE = 200070000; % J/K Electrolyser heat capacity 
R_tAWE = 0.0001133; %K/W Thermal resistance
n_cAWE = 6669; %Number of cells
TminAWE = 313; %K Minimal operating temperature

% Heel belangrijk hier: assumpties gemaatk dat het in 20min gebeurt en altijd
% opwarmen van 20°C naar 40°C
% Er wordt echter tijdens startup geen winst gemaakt maar verlies. Hierdoor wordt voor 
% iets hogere electriciteitsprijzen de start up cost negatief
% De startup cost wordt negatiever voor hoge electriciteitsprijzen, lage waterstofprijzen
% en wordt negatiever bij lager vermogen direct bij het opstarten wat eigenlijk een beslissingsvariabele is van het model (zie ander matlab bestand)

% Hier oppervlakkig toegepast voor half en full power, in visual studio code wel als beslissingsvariabele!!


pi_e = (0:1:150)*10^-6;
pMax20deg = 2368.32;
pMax40deg = 2206.68;
eta_fMax20deg = 0.966522;
eta_fMax40deg = 0.96617;

Q_CS = C_hAWE*(40-20)/1200 + (30-20)/R_tAWE; % in watt!
m_H2_cost = eta_fMax20deg*0.0020166669*n_cAWE*iMaxAWE*0.21/2/96485.3;
m_H2_revenue = eta_fMax40deg*0.0020166669*n_cAWE*iMaxAWE*0.21/2/96485.3;
m_H2O = m_H2_revenue*0.01801528/0.0020166669;


for pi_H = 2:7 %€/kg
    Cost1_e(pi_H-1,:) = Q_CS*1200/3600/0.95*pi_e;
    Cost2(pi_H-1,:) = m_H2_cost*3600*pi_H - pMax20deg*n_cAWE*pi_e;
    Revenue(pi_H-1,:) = m_H2_revenue*(3600-1200)*pi_H - pMax40deg*2/3*n_cAWE*pi_e - m_H2O*4184*(40-20)/0.95*2/3*pi_e;
end
Cost_CS = Cost1_e + Cost2 - Revenue;
figure
plot(pi_e*10^6, Cost_CS)
xlabel("Electricity price [€/MWh]")
ylabel("Cost cold start up")
legend("Hydrogen price = 2 [€/kg]","Hydrogen price = 3 [€/kg]","Hydrogen price = 4 [€/kg]","Hydrogen price = 5 [€/kg]","Hydrogen price = 6 [€/kg]","Hydrogen price = 7 [€/kg]")
title("Cost cold start up AWE, immediately at maximal power")

pi_e = (0:1:150)*10^-6;
iHalf = 1985.71; %A/m²
pHalf20deg = 1054.18;
pHalf40deg = 981.747;
eta_fHalf20deg = 0.914043;
eta_fHalf40deg = 0.913372;

Q_CS = C_hAWE*(40-20)/1200 + (30-20)/R_tAWE; % in watt!
m_H2_cost = eta_fHalf20deg*0.0020166669*n_cAWE*iHalf*0.21/2/96485.3;
m_H2_revenue = eta_fHalf40deg*0.0020166669*n_cAWE*iHalf*0.21/2/96485.3;
m_H2O = m_H2_revenue*0.01801528/0.0020166669;


for pi_H = 2:7 %€/kg
    Cost1_e(pi_H-1,:) = Q_CS*1200/3600/0.95*pi_e;
    Cost2(pi_H-1,:) = m_H2_cost*3600*pi_H - pHalf20deg*n_cAWE*pi_e;
    Revenue(pi_H-1,:) = m_H2_revenue*(3600-1200)*pi_H - pHalf40deg*2/3*n_cAWE*pi_e - m_H2O*4184*(40-20)/0.95*2/3*pi_e;
end
Cost_CS = Cost1_e + Cost2 - Revenue;
figure
plot(pi_e*10^6, Cost_CS)
xlabel("Electricity price [€/MWh]")
ylabel("Cost cold start up")
legend("Hydrogen price = 2 [€/kg]","Hydrogen price = 3 [€/kg]","Hydrogen price = 4 [€/kg]","Hydrogen price = 5 [€/kg]","Hydrogen price = 6 [€/kg]","Hydrogen price = 7 [€/kg]")
title("Cost cold start up AWE, immediately at half power")


%% AWE with nuclear heat 
pi_H = 4; %€/kg
pi_e_average = 4.324e-5; %€/Wh
pi_e_min = 3.702e-5; %€/Wh


Q_CS = C_hAWE*(40-20)/600 + (30-20)/R_tAWE;
m_H2_cost = 0.9665*0.0020166669*n_cAWE*iMaxAWE*0.21/2/96485.3;
m_H2_revenue = 0.9662*0.0020166669*n_cAWE*iMaxAWE*0.21/2/96485.3;
m_H2O = m_H2_revenue*0.01801528/0.0020166669;

Cost1a_e = Q_CS*600/3600*0.45*pi_e_average;
Cost2a = m_H2_cost*3600*pi_H - 2412.72*n_cAWE*pi_e_average;
Revenue_a = m_H2_revenue*3000*pi_H - 2249.27*5/6*n_cAWE*pi_e_average - m_H2O*4184*(40-20)*0.45*5/6*pi_e_average;

Cost1b = Q_CS*600/3600*0.45*pi_e_min;
Cost2b = m_H2_cost*3600*pi_H - 2412.72*n_cAWE*pi_e_min;
Revenue_b = m_H2_revenue*3000*pi_H - 2249.27*5/6*n_cAWE*pi_e_min - m_H2O*4184*(40-20)*0.45*5/6*pi_e_min;

Cost_CSa = Cost1a_e + Cost2a - Revenue_a;
Cost_CSb = Cost1b + Cost2b - Revenue_b;

%% PEM
% cold startup cost = 0?!

iMinPEM = 600; %Minimal current density A/m²
iMaxPEM = 20000; %Maximal current density A/m²
U_tnPEM = 1.4813; %Thermonteutral voltage V (paper 66)
T_sPEM = 373; %Normal working temperature K
C_HSPEM = 0; %Estimated hot start cost €
C_CSPEM = 98; %Estimated cold start cost €
T_0PEM = 373; %Initial temperature of electrolyser K
C_hPEM = 35280000; % J/K Electrolyser heat capacity 
R_tPEM = 0.00010666; %K/W Thermal resistance
n_cPEM = 1176; %Number of cells
TminPEM = 293; %K Minimal operating temperature


T = 100 + 3600/C_hPEM*(n_cPEM*1355.0517022112194 - n_cPEM*934.2627453262709*U_tnPEM  - (100-20)/R_tPEM)



