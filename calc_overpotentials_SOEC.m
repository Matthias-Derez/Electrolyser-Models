function [U_total,U_rev, U_act, U_ohm, U_conc, Power, j_0a, j_0c] = calc_overpotentials_SOEC(jmin,jmax,Tmin,Tmax, Area, Tstep, figures)
%Hier voorlopig eenheden van 84 gevolgd, opletten, Urev in bar, voor Uconc
%moet je met Pa werken dus maal 10^5
% in verslag overal in Pa en by Urev*10^-5, bij Uconc geen *10^5

Trange = Tmin:Tstep:Tmax; %K
jrange = jmin:(jmax-jmin)/((Tmax-Tmin)/Tstep):jmax; %A/m²
R = 8.314; %J/(mol K)
F = 96485.3; %sA/mol
P = 1; %bar
pH2 = 0.5*P; %bar
pH2O = 0.5*P; %bar
A = Area; %m² (22)
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
    index_temp = 0;
    index_current = index_current + 1;
    for T = Trange
        index_temp = index_temp + 1;
        j_0a = gamma_a*exp(-E_acta/(R*T)); %A/m² (84) WAARDES VOOR EXCHANGE CURRENT DENSITIES IN PAPER (39)
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

if figures 
    Tmin_celcius = Tmin -273;
    Tmax_celcius = Tmax - 273;
    
    Trange_celcius = Tmin_celcius:Tstep:Tmax_celcius; %°C
    jrange = jmin:(jmax-jmin)/((Tmax-Tmin)/Tstep):jmax; %A/m²
    % Power = calc_power(jrange,Trange);
    figure(1)
    h = surf(Trange_celcius, jrange, Power); %Beetje vreemd surf functie plot jrange op y as en Trange op x as terwijl Power(index_current, index_temp)
    xlabel("Temperature [°C]", FontSize=12)
    ylabel("Current density [A/m²]",FontSize=12)
    zlabel("Power [W]",FontSize=10)
    set(h,'LineStyle','none')
    title("Power consumption per cell of an SOE",FontSize=14)
    view(30,40)
    grid on
    yticklabels(num2str(get(gca, 'YTick').'));
    xh = get(gca,'XLabel'); % Handle of the x label
    set(xh, 'Units', 'Normalized')
    pos = get(xh, 'Position');
    set(xh, 'Position',pos.*[1,-0.0,1],'Rotation',-10)
    yh = get(gca,'YLabel'); % Handle of the y label
    set(yh, 'Units', 'Normalized')
    pos = get(yh, 'Position');
    set(yh, 'Position',pos.*[1.12,-0.7,1],'Rotation',40)
    print -depsc Power_SOEC.eps
    
    counter_current = 1:100:length(Trange_celcius);
    counter_temp = 1:length(Trange_celcius);
    
    figure(2),hold on
    plot(Trange_celcius, U_total(counter_current,counter_temp), 'LineWidth', 1)
    axis tight
    xlabel('Temperature  [°C]')
    ylabel('U_{total} [Volt]')
    legend('Current density=' + string(jrange(counter_current(1))) + ' A/m²','Current density=' + string(jrange(counter_current(2))) + ' A/m²','Current density=' + string(jrange(counter_current(3))) + ' A/m²','Current density=' + string(jrange(counter_current(4))) + ' A/m²','Current density=' + string(jrange(counter_current(5))) + ' A/m²');
    hold off
    print -depsc U_tot.eps
    
    figure(3),hold on
    plot(Trange_celcius, U_rev(counter_current,counter_temp), 'LineWidth', 1)
    axis tight
    xlabel('Temperature  [°C]')
    ylabel('U_{rev} [Volt]')
    legend('Current density=' + string(jrange(counter_current(1))) + ' A/m²','Current density=' + string(jrange(counter_current(2))) + ' A/m²','Current density=' + string(jrange(counter_current(3))) + ' A/m²','Current density=' + string(jrange(counter_current(4))) + ' A/m²','Current density=' + string(jrange(counter_current(5))) + ' A/m²');
    hold off
    print -depsc U_rev.eps
    
    figure(4),hold on
    plot(Trange_celcius, U_act(counter_current,counter_temp), 'LineWidth', 1)
    axis tight
    xlabel('Temperature  [°C]')
    ylabel('U_{act} [Volt]')
    legend('Current density=' + string(jrange(counter_current(1))) + ' A/m²','Current density=' + string(jrange(counter_current(2))) + ' A/m²','Current density=' + string(jrange(counter_current(3))) + ' A/m²','Current density=' + string(jrange(counter_current(4))) + ' A/m²','Current density=' + string(jrange(counter_current(5))) + ' A/m²');
    hold off
    print -depsc U_act.eps
    
    figure(5),hold on
    plot(Trange_celcius, U_ohm(counter_current,counter_temp), 'LineWidth', 1)
    axis tight
    xlabel('Temperature  [°C]')
    ylabel('U_{ohm} [Volt]')
    legend('Current density=' + string(jrange(counter_current(1))) + ' A/m²','Current density=' + string(jrange(counter_current(2))) + ' A/m²','Current density=' + string(jrange(counter_current(3))) + ' A/m²','Current density=' + string(jrange(counter_current(4))) + ' A/m²','Current density=' + string(jrange(counter_current(5))) + ' A/m²');
    hold off
    print -depsc U_ohm.eps
    
    
    figure(6), hold on
    plot(Trange_celcius, U_conc(counter_current,counter_temp), 'LineWidth', 1)
    axis tight
    xlabel('Temperature  [°C]')
    ylabel('U_{conc} [Volt]')
    legend('Current density=' + string(jrange(counter_current(1))) + ' A/m²','Current density=' + string(jrange(counter_current(2))) + ' A/m²','Current density=' + string(jrange(counter_current(3))) + ' A/m²','Current density=' + string(jrange(counter_current(4))) + ' A/m²','Current density=' + string(jrange(counter_current(5))) + ' A/m²', 'Location', 'northwest');
    hold off
    print -depsc U_conc.eps
    
    counter_temp = 1:100:length(Trange_celcius);
    counter_current = 1:length(Trange_celcius);
    figure(7),hold on
    plot(jrange, U_total(counter_current,counter_temp), 'LineWidth', 1)
    axis tight
    xlabel('Current density [A/m²]')
    ylabel('U_{total} [Volt]')
    legend('Temperature =' + string(Trange_celcius(counter_temp(1))) + ' °C','Temperature=' + string(Trange_celcius(counter_temp(2))) + ' °C','Temperature=' + string(Trange_celcius(counter_temp(3))) + ' °C','Temperature=' + string(Trange_celcius(counter_temp(4))) + ' °C','Temperature=' + string(Trange_celcius(counter_temp(5))) + ' °C');
    hold off
    print -depsc U_tot.eps
    
    figure(8),hold on
    plot(jrange, U_rev(counter_current,counter_temp), 'LineWidth', 1)
    axis tight
    xlabel('Current density [A/m²]')
    ylabel('U_{rev} [Volt]')
    legend('Temperature =' + string(Trange_celcius(counter_temp(1))) + ' °C','Temperature=' + string(Trange_celcius(counter_temp(2))) + ' °C','Temperature=' + string(Trange_celcius(counter_temp(3))) + ' °C','Temperature=' + string(Trange_celcius(counter_temp(4))) + ' °C','Temperature=' + string(Trange_celcius(counter_temp(5))) + ' °C');
    hold off
    print -depsc U_rev.eps
    
    figure(9),hold on
    plot(jrange, U_act(counter_current,counter_temp), 'LineWidth', 1)
    axis tight
    xlabel('Current density [A/m²]')
    ylabel('U_{act} [Volt]')
    legend('Temperature =' + string(Trange_celcius(counter_temp(1))) + ' °C','Temperature=' + string(Trange_celcius(counter_temp(2))) + ' °C','Temperature=' + string(Trange_celcius(counter_temp(3))) + ' °C','Temperature=' + string(Trange_celcius(counter_temp(4))) + ' °C','Temperature=' + string(Trange_celcius(counter_temp(5))) + ' °C');
    hold off
    print -depsc U_act.eps
    
    figure(10),hold on
    plot(jrange, U_ohm(counter_current,counter_temp), 'LineWidth', 1)
    axis tight
    xlabel('Current density [A/m²]')
    ylabel('U_{ohm} [Volt]')
    legend('Temperature =' + string(Trange_celcius(counter_temp(1))) + ' °C','Temperature=' + string(Trange_celcius(counter_temp(2))) + ' °C','Temperature=' + string(Trange_celcius(counter_temp(3))) + ' °C','Temperature=' + string(Trange_celcius(counter_temp(4))) + ' °C','Temperature=' + string(Trange_celcius(counter_temp(5))) + ' °C');
    hold off
    print -depsc U_ohm.eps
    
    
    figure(11), hold on
    plot(jrange, U_conc(counter_current,counter_temp), 'LineWidth', 1)
    axis tight
    xlabel('Current density [A/m²]')
    ylabel('U_{conc} [Volt]')
    legend('Temperature =' + string(Trange_celcius(counter_temp(1))) + ' °C','Temperature=' + string(Trange_celcius(counter_temp(2))) + ' °C','Temperature=' + string(Trange_celcius(counter_temp(3))) + ' °C','Temperature=' + string(Trange_celcius(counter_temp(4))) + ' °C','Temperature=' + string(Trange_celcius(counter_temp(5))) + ' °C');
    hold off
    print -depsc U_conc.eps
    IMatrix = zeros(length(Power));
    TMatrix = zeros(length(Power));
    for i = 1:length(Power)
        for j = 1:length(Power)
            IMatrix(i,j) = 0.21*jrange(i);
            TMatrix(i,j) = Trange(j);
        end
    end

    figure(12), hold on
    Z = 5776*Power-1.2955*IMatrix*5776;
    h = surf(Trange_celcius, jrange, Z);
    set(h,'LineStyle','none')
    title("Power*n_c - U_{tn}*I*n_c - (T - T_a)/R_t")
    view(30,40)
    grid on
    xh = get(gca,'XLabel'); % Handle of the x label
    set(xh, 'Units', 'Normalized')
    pos = get(xh, 'Position');
    set(xh, 'Position',pos.*[1,-0.0,1],'Rotation',-10)
    yh = get(gca,'YLabel'); % Handle of the y label
    set(yh, 'Units', 'Normalized')
    pos = get(yh, 'Position');
    set(yh, 'Position',pos.*[1.12,-0.7,1],'Rotation',40)
   % Create a color map based on Z values
    colormap([1 0 0; 0 1 0])
    caxis([-1 1])
    Z_color = Z >= 0; % Assign a color to each point based on the Z value
    set(h, Z_color); % Assign the color map to the surface plot
end
