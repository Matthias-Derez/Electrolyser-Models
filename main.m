%% General variables
close all
clearvars

Area = 0.21; %m²

% n_c = 3800; 
%% SOEC: Calculating overpotentials
jminSOEC = 2000; %A/m²
jmaxSOEC = 10000; %A/m²
TminSOEC = 800+273; %K
TmaxSOEC = 1000+273; %K
[U_totalSOEC,U_revSOEC, U_actSOEC, U_ohmSOEC, U_concSOEC,PowerSOEC, j_0aSOEC, j_0cSOEC] = calc_overpotentials_SOEC(jminSOEC,jmaxSOEC,TminSOEC,TmaxSOEC,Area,0.5,0);
%% SOEC: Calculating overpotentials using ASR empirical model
[PowerASR] = calc_overpotentials_ASR(jminSOEC,jmaxSOEC,TminSOEC,TmaxSOEC,PowerSOEC,Area, 0.5,0);
%% SOEC: Calculating mass rate of produced hydrogen per cell
[mH2_SOEC, eff_total_SOEC] = calc_mH2_SOEC(jminSOEC, jmaxSOEC, TminSOEC, TmaxSOEC, PowerSOEC, Area,0.5,1);

%% SOEC: Linearizing power curve
errorvectorSOEC = zeros(80,1);
rel_errorvectorSOEC = zeros(80,1);
for i = 1:80
    [coeffsSOEC, error, relative_error] = linearize_power(jminSOEC,jmaxSOEC,TminSOEC,TmaxSOEC,PowerSOEC,i, 'SOEC',Area, 0.5,0);
    errorvectorSOEC(i) = error; 
    rel_errorvectorSOEC(i) = relative_error;
end
figure(101)
plot(errorvectorSOEC(1:end))
xlabel('Number of segments')
ylabel('Mean error')
title('Mean linearization error SOEC')
figure(102)
plot(rel_errorvectorSOEC(1:end))
xlabel('Number of segments')
ylabel('Mean relative error')
title('Mean relative linearization error SOEC')
figure(99)
plot(log10(rel_errorvectorSOEC(1:end)))
xlabel('Number of segments')
ylabel('log(Mean relative error)')
title('Mean relative linearization error SOEC')
%% SOEC: Linearizing power curve 2
errorvectorSOEC = zeros(10,1);
rel_errorvectorSOEC = zeros(10,1);
for i = 1:10
    [coeffsSOEC, error, relative_error] = linearize_power2(jminSOEC,jmaxSOEC,TminSOEC,TmaxSOEC,PowerSOEC,3, 'SOE',Area, 0.5,1);
    errorvectorSOEC(i) = error; 
    rel_errorvectorSOEC(i) = relative_error;
end
num_segments = 1:10;
figure(115)
plot(num_segments.^2, errorvectorSOEC)
xlabel('Number of segments')
ylabel('Mean error')
title('Mean linearization error SOEC')
figure(116)
plot(num_segments.^2, rel_errorvectorSOEC)
xlabel('Number of segments')
ylabel('Mean relative error [-]')
title('Mean relative linearization error for the power curve of a SOEC')
figure(117)
semilogy(num_segments.^2, rel_errorvectorSOEC)
xlabel('Number of segments', FontSize=12)
ylabel('Mean relative error [-]', FontSize=12)
title(sprintf('Mean relative linearization error\n for the power curve of an SOE'), FontSize=14)
%% AWE: Calculating overpotentials
jminAWE = 1000; %A/m²
jmaxAWE = 4000; %A/m²
TminAWE = 40+273; %K
TmaxAWE = 90+273; %K
[U_totalAWE,U_revAWE, U_actAWE, U_ohmAWE, U_concAWE,PowerAWE] = calc_overpotentials_AWE(jminAWE,jmaxAWE,TminAWE,TmaxAWE,Area,0.5,0);

%% AWE: Calculating overpotentials using Ulleberg
[U_totalUlle, U_revUlle, Power_densityUlle, PowerUlle] = calc_overpotentials_Ulleberg(jminAWE,jmaxAWE,TminAWE,TmaxAWE, PowerAWE, Area,0.5,1);

%% AWE: Faraday efficiency: 5 parameter model
eff_Farad = calc_eff_Farad(jminAWE,jmaxAWE,TminAWE,TmaxAWE,0.5,1);

%% AWE: Calculating mass rate of produced hydrogen per cell
[mH2_AWE, eff_total_AWE] = calc_mH2_AWE(jminAWE,jmaxAWE,TminAWE,TmaxAWE,eff_Farad, PowerAWE,Area,0.5,1);

%% AWE: Linearizing power cur
errorvectorAWE = zeros(30,1);
rel_errorvectorAWE = zeros(30,1);
for i = 1:30
    [coeffsAWE, error, relative_error] = linearize_power(jminAWE,jmaxAWE,TminAWE,TmaxAWE,PowerAWE, i, 'AWE',Area, 0.5,0);
    errorvectorAWE(i) = error;
    rel_errorvectorAWE(i) = relative_error;
end

figure(103)
plot(errorvectorAWE(1:end))
xlabel('Number of segments')
ylabel('log(Mean error)')
title('Mean linearization error AWE')
figure(104)
plot(rel_errorvectorAWE(1:end))
xlabel('Number of segments')
ylabel('Mean relative error')
title('Mean relative linearization error AWE')
figure(105)
semilogy(rel_errorvectorAWE)
xlabel('Number of segments')
ylabel('Mean relative error [-]')
title('Mean relative linearization error for the power curve of a alkaline electrolyser')
%% AWE: Linearizing power curve 2
errorvectorAWE = zeros(10,1);
rel_errorvectorAWE = zeros(10,1);
for i = 1:10
    [coeffsAWE, error, relative_error] = linearize_power2(jminAWE,jmaxAWE,TminAWE,TmaxAWE,PowerAWE, i, 'AWE',Area, 0.5,0);
    errorvectorAWE(i) = error;
    rel_errorvectorAWE(i) = relative_error;
end
num_segments = 1:10;
figure(118)
plot(errorvectorAWE(1:end))
xlabel('Number of segments')
ylabel('log(Mean error)')
title('Mean linearization error AWE')
figure(119)
plot(rel_errorvectorAWE(1:end))
xlabel('Number of segments')
ylabel('Mean relative error')
title('Mean relative linearization error AWE')
figure(120)
semilogy(num_segments.^2, rel_errorvectorAWE)
xlabel('Number of segments', FontSize=12)
ylabel('Mean relative error [-]', FontSize=12)
title(sprintf('Mean relative linearization error\n for the power curve of an AWE'), FontSize=14)
%% AWE: Linearizing Faraday efficiency
errorvectorFARAD = zeros(10,1);
rel_errorvectorFARAD = zeros(10,1);
for i = 1:10
    [coeffs_eff_Farad, error, relative_error] = linearize_eff_Farad(jminAWE,jmaxAWE,TminAWE,TmaxAWE, eff_Farad, 1,0.5,1);
    errorvectorFARAD(i) = error;
    rel_errorvectorFARAD(i) = relative_error;
end
num_segments = 1:10;

figure(106)
plot(errorvectorFARAD(1:end))
xlabel('Number of segments')
ylabel('Mean error')
title('Mean linearization error Faraday efficiency')
figure(107)
plot(rel_errorvectorFARAD(1:end))
xlabel('Number of segments')
ylabel('Mean relative error')
title('Mean relative linearization error Faraday efficiency')
figure(108)
semilogy(num_segments.^2, rel_errorvectorFARAD)
xlabel('Number of segments')
ylabel('Mean relative error [-]')
title({'Mean relative linearization error'},{'for \eta_F*I of a alkaline electrolyser'})

%% AWE: Linearizing product (Faraday efficiency)*Current 2
errorvectorFARAD = zeros(10,1);
rel_errorvectorFARAD = zeros(10,1);
for i = 1:10
    [coeffs_eff_Farad, error, relative_error] = linearize_eff_Farad_current2(jminAWE,jmaxAWE,TminAWE,TmaxAWE, eff_Farad, 1,Area, 0.5,1);
    errorvectorFARAD(i) = error;
    rel_errorvectorFARAD(i) = relative_error;
end
num_segments = 1:10;

figure(109)
plot(errorvectorFARAD(1:end))
xlabel('Number of segments')
ylabel('Mean error')
title('Mean linearization error Faraday efficiency')
figure(110)
plot(rel_errorvectorFARAD(1:end))
xlabel('Number of segments')
ylabel('Mean relative error')
title('Mean relative linearization error Faraday efficiency')
figure(111)
semilogy(num_segments.^2, rel_errorvectorFARAD)
xlabel('Number of segments')
ylabel('Mean relative error [-]')
title({'Mean relative linearization error'},{'for \eta_F*I of a alkaline electrolyser'})
%% PEM: Calculating overpotentials
jminPEM = 1500; %A/m²
jmaxPEM = 20000; %A/m²
TminPEM = 20+273; %K
TmaxPEM = 100+273; %K
[U_totalPEM,U_revPEM, U_actPEM, U_ohmPEM, U_concPEM, PowerPEM] = calc_overpotentials_PEM(jminPEM,jmaxPEM,TminPEM,TmaxPEM,Area,0.5,1);

%% PEM: Calculating overpotentials using Ulleberg
[U_totalUlle, U_revUlle, Power_densityUlle, PowerUlle] = calc_overpotentials_Ulleberg(jminPEM,jmaxPEM,TminPEM,TmaxPEM, PowerPEM, Area,0.5,1);

%% PEM: Calculating mass rate of produced hydrogen per cell
[mH2_PEM, eff_total_PEM] = calc_mH2_PEM(jminPEM, jmaxPEM, TminPEM, TmaxPEM, PowerPEM, Area,0.5,1);

%% PEM: Linearizing power curve
errorvectorPEM = zeros(30,1);
rel_errorvectorPEM = zeros(30,1);
for i = 1:30
    [coeffsPEM, error, relative_error] = linearize_power(jminPEM,jmaxPEM,TminPEM,TmaxPEM,PowerPEM,i, 'PEM',Area, 0.5,0);
    errorvectorPEM(i) = error;
    rel_errorvectorPEM(i) = relative_error;
end
figure(112)
plot(errorvectorPEM(1:end))
xlabel('Number of segments')
ylabel('Mean error')
title('Linearization error PEM')
figure(113)
plot(rel_errorvectorPEM(1:end))
xlabel('Number of segments')
ylabel('Mean relative error')
title('Mean relative linearization error PEM')
figure(114)
semilogy(rel_errorvectorPEM(1:end))
xlabel('Number of segments')
ylabel('log(Mean relative error)')
title('Mean relative linearization error PEM')
%% PEM: Linearizing power curve 2
errorvectorPEM = zeros(10,1);
rel_errorvectorPEM = zeros(10,1);
for i = 1:10
    [coeffsPEM, error, relative_error] = linearize_power2(jminPEM,jmaxPEM,TminPEM,TmaxPEM,PowerPEM,i, 'PEM',Area, 0.5,0);
    errorvectorPEM(i) = error;
    rel_errorvectorPEM(i) = relative_error;
end
num_segments = 1:10;

figure(121)
plot(errorvectorPEM(1:end))
xlabel('Number of segments')
ylabel('Mean error')
title('Linearization error PEM')
figure(122)
plot(rel_errorvectorPEM(1:end))
xlabel('Number of segments')
ylabel('Mean relative error')
title('Mean relative linearization error PEM')
figure(123)
semilogy(num_segments.^2, rel_errorvectorPEM)
xlabel('Number of segments', FontSize=12)
ylabel('Mean relative error [-]', FontSize=12)
title(sprintf('Mean relative linearization error\n for the power curve of a PEM electrolyser'), FontSize=14)