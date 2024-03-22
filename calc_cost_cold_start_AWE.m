% Belangrijke opmerking: hier berekend zoals in visual studio code met waarden voor j = 2000A/m²
% en verdere assumpties van T 20-> 40°C

% eta_f_product_I (@ 20°C, 1958.71A/m²) = 387.118
a_f = -0.352484;
%  power (@20°C, 1958.71A/m²) = 1087.39W
a = -3.945498;

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
M_H2 = 0.0020166669;
F = 96485.3;


pi_H = 7; %€/kg
pi_e = (0:1:150)*10^-6;


% eta_f_product_I (@ 20°C, 1958.71A/m²) = 387.118
a_f = -0.352484;

m_H2_20 = 387.118*M_H2*n_cAWE/(2*F);
m_H2_40 = (387.118+a_f*20)*M_H2*n_cAWE/(2*F);


%  power (@20°C, 1958.71A/m²) = 1087.39W
a = -3.945498;

p_20 = 1087.39;
p_40 = p_20+a*20;



Q_CS = C_hAWE*(40-20)/600 + (30-20)/R_tAWE; % in watt!

Cost1_e = Q_CS*600/3600/0.95*pi_e;
Cost1_h = Q_CS*600/3600*0.45*pi_e;
Cost2 = m_H2_20*3600*pi_H - p_20*n_cAWE*pi_e;
Revenue = m_H2_40*3000*pi_H - p_40*3000/3600*n_cAWE*pi_e; 

Cost_CS_e = Cost1_e + Cost2 - Revenue;
Cost_CS_h = Cost1_h + Cost2 - Revenue;
figure, hold on
plot(pi_e*10^6, Cost_CS_e)
plot(pi_e*10^6, Cost_CS_h)
xlabel("Electricity price [€/MWh]")
ylabel("Cost cold start up [€]")
legend("Heating up with electricity", "Heating up with nuclear heat")



pi_H = 7; %€/kg
pi_e = (0:1:150)*10^-6;


% eta_f_product_I (@ 20°C, 1958.71A/m²) = 387.118
a_f = -0.352484;

m_H2_20 = 387.118*M_H2*n_cAWE/(2*F);
m_H2_40 = (387.118+a_f*20)*M_H2*n_cAWE/(2*F);


%  power (@20°C, 1958.71A/m²) = 1087.39W
a = -3.945498;

p_20 = 1087.39;
p_40 = p_20+a*20;



Q_CS = C_hAWE*(40-20)/1800 + (30-20)/R_tAWE; % in watt!

Cost1_e = Q_CS*1800/3600/0.95*pi_e;
Cost1_h = Q_CS*1800/3600*0.45*pi_e;
Cost2 = m_H2_20*3600*pi_H - p_20*n_cAWE*pi_e;
Revenue = m_H2_40*3000*pi_H - p_40*3000/3600*n_cAWE*pi_e; 

Cost_CS_e = Cost1_e + Cost2 - Revenue;
Cost_CS_h = Cost1_h + Cost2 - Revenue;
figure, hold on
plot(pi_e*10^6, Cost_CS_e)
plot(pi_e*10^6, Cost_CS_h)
xlabel("Electricity price [€/MWh]")
ylabel("Cost cold start up [€]")
legend("Heating up with electricity", "Heating up with nuclear heat")

%% Hot start duurt 20 min (122)


pi_H = 7; %€/kg
pi_e = (0:1:150)*10^-6;


% eta_f_product_I (@ 20°C, 1958.71A/m²) = 387.118
a_f = -0.352484;

m_H2_20 = 387.118*M_H2*n_cAWE/(2*F);
m_H2_40 = (387.118+a_f*20)*M_H2*n_cAWE/(2*F);


%  power (@20°C, 1958.71A/m²) = 1087.39W
a = -3.945498;

p_20 = 1087.39;
p_40 = p_20+a*20;


Cost2 = m_H2_20*3600*pi_H - p_20*n_cAWE*pi_e;
Revenue = m_H2_40*3000*pi_H - p_40*3000/3600*n_cAWE*pi_e; 

Cost_CS_e = Cost1_e + Cost2 - Revenue;
Cost_CS_h = Cost1_h + Cost2 - Revenue;
figure, hold on
plot(pi_e*10^6, Cost_CS_e)
plot(pi_e*10^6, Cost_CS_h)
xlabel("Electricity price [€/MWh]")
ylabel("Cost cold start up [€]")
legend("Heating up with electricity", "Heating up with nuclear heat")



