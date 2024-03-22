%% POGING OPSTARTCOST EXOGEEN AAN MODEL

% @40°C
p_1985 = 981.747*6669;
p_4000 = 2206.68*6669;

m_H2_1985 = 3.97908e-6*6669;
m_H2_4000 = 8.47875e-6*6669;

C_HS_1985 = (m_H2_1985*3600*pi_H-p_1985*pi_e)/6;
C_HS_4000 = (m_H2_4000*3600*pi_H-p_4000*pi_e)/6;
figure, hold on
plot(pi_e*10^6, C_HS_1985)
plot(pi_e*10^6, C_HS_4000)
xlabel("Electricity price [€/MWh]")
ylabel("Cost hot start up [€]")
legend("Half power", "full power")

% @65°C
p_1985 = 904.903*6669;
p_4000 = 2045.67*6669;

m_H2_1985 = 3.89582e-6*6669;
m_H2_4000 = 8.39021e-6*6669;

pi_H = 2.6; %€/kg
pi_e = (0:1:150)*10^-6;

C_HS_1985 = (m_H2_1985*3600*pi_H-p_1985*pi_e)/6;
C_HS_4000 = (m_H2_4000*3600*pi_H-p_4000*pi_e)/6;
figure, hold on
plot(pi_e*10^6, C_HS_1985)
plot(pi_e*10^6, C_HS_4000)
xlabel("Electricity price [€/MWh]")
ylabel("Cost hot start up [€]")
legend("Half power", "full power")

% @90°C
p_1985 = 838.01*6669;
p_4000 = 1911.39*6669;

m_H2_1985 = 3.72951e-6*6669;
m_H2_4000 = 8.39021e-6*6669;

pi_H = 2.6; %€/kg
pi_e = (0:1:150)*10^-6;

C_HS_1985 = (m_H2_1985*3600*pi_H-p_1985*pi_e)/6;
C_HS_4000 = (m_H2_4000*3600*pi_H-p_4000*pi_e)/6;
figure, hold on
plot(pi_e*10^6, C_HS_1985)
plot(pi_e*10^6, C_HS_4000)
xlabel("Electricity price [€/MWh]")
ylabel("Cost hot start up [€]")
legend("Half power", "full power")
