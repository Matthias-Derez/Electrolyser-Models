function [U_total,U_rev, U_act, U_ohm, U_conc, Power] = calc_overpotentials_AWE(jmin,jmax,Tmin,Tmax, Area, Tstep, figures)
% Based on model described by Irem Firtina-Ertis "Thermodynamic and
% electrochemical assessment of an alkaline electrolyzer (AE) at different
% operating parameters"  (paper 95)

Trange = Tmin:Tstep:Tmax; %K
jrange = jmin:(jmax-jmin)/((Tmax-Tmin)/Tstep):jmax; %A/m²

R = 8.314; %J/(mol K)
F = 96485.3; %sA/mol
P = 1; %bar

% Gegevens (95)
w = 30; %wt% KOH in electrolyte
d_electrode = 0.003; %m
d_electrolyte = 0.0025; %m
d_membrane = 0.0003; %m
%A = 0.05; %m² (veel groter dan gebruikt voor SOEC)
A = Area; %m² Voorlopig zelfde waarde gebruikt als voor SOEC, niet waarde van paper (95)
n_c = 49; % number of cells
sigma_membrane = 18.75; %Siemens/m = 1/(Ohm*m)
sigma_nickel = 11808768; %Siemens/m

%gegevens (111)
%m = w*10*10/7/56.1056; %mol/kg, Molality of 30wt% KOH solution in water
m = (w/100*1/(56.1056e-3))/(1-w/100*1); %molKOH/kgSOLVENT, Molality of 30wt% KOH solution in water
% Gegevens (101)
gamma_Man = 1.25;
gamma_Mcat = 1.05;
Tref = 298; %K
% j_0aref = 10^-7;%A/m²
% j_0aref = 10^-7; %A/m² OPM dit was oorspronkelijk 10^-7, maar als we hier 10^-9 van maken ligt het dichter bij ulleberg,
%                     % daarnaast valt dit nog steeds binnen de range gegeven
%                     % in (105) die eigenlijk voor PEM is, maar  (98) toont
%                     % dat voor de AWE de exchange current density eigenlijk
%                     % nog lager is dan voor PEM
% j_0cref = 10^1; %A/m²

j_0aref = 10^-5; %A/m² (101)
j_0cref = 10^1; %A/m²
% j_0aref = 0.048; (98)
% j_0cref = 1.9;
% Niet de transfer coefficients van (101) gebruikt, data is daar gefit,
% formules van (95) komen beter overeen met de rest van de literatuur voor
% charge transfer coefficients
% alfa_an = 1.65; 
% alfa_cat = 0.73;
DeltaG = 80510; %J/mol


pH2 = 1; %bar
pO2 = 1; %bar



U_total = zeros(length(jrange), length(Trange));
U_rev = zeros(length(jrange), length(Trange));
U_act = zeros(length(jrange), length(Trange));
U_ohm = zeros(length(jrange), length(Trange));
U_conc = zeros(length(jrange), length(Trange));
Power = zeros(length(jrange), length(Trange));
index_current = 0;
j_0c = zeros(length(Trange),1);
j_0a = zeros(length(Trange),1);
U_actc = zeros(length(jrange), length(Trange));
U_acta = zeros(length(jrange), length(Trange));
U_act = zeros(length(jrange), length(Trange));
for j = jrange
    index_current = index_current + 1;
    indextemp = 0;
    for T = Trange
        indextemp = indextemp + 1;
        % Paper (95)
        rho_T = 868.1+6.81*w+1.137*T+0.04391*w^2 + 0.0138*w*T- 0.002487*T^2 + 0.0001279*w^3 -2.547 * 10^-5* w^2*T-2.644 *10^-5*w*T^2;
        M = w*rho_T/(56.105 * 100);
        sigma_electrolyte = (-2.041*M -0.0028*M.^2 + 0.005332*M.*T + 207.2*M/T + 0.001043*M^3 -3* 10^-7*M^2*T^2 )* 100;
        U_ohm(index_current,indextemp) = j*(2*d_electrode/sigma_nickel + d_membrane/sigma_membrane + d_electrolyte/sigma_electrolyte);
        alfa_cat = 0.1175 + 0.00095*T;
        alfa_an = 0.0675 + 0.00095*T;
        
        % Wat er in (95) staat
        U_conc(index_current,indextemp) = -R*T/(2*F)*log(1-j/(jmax*1.05))*4; %opm, hier (jmax+1) om te voorkomen dat potentiaal naar oneindig gaat op j = jmax
        % hier extra maal factor 4 wat niet in (95) staat, is wel gedaan in
        % (114) en (113) + (95) verwijst naar (113) 


        % (95) heeft dit van (113) (een bron voor PEM)


       

        %(101)
%         alfa_an = 1.65; 
%         alfa_cat = 0.73;
%         alfa_an = 1.85; 
%         alfa_cat = 0.85;
        
        % Paper (111)
        %aw(indextemp) = 10^(-0.02255*m+0.001434*m^2+(1.38*m-0.9254*m^2)/T); %water activity
        aw(indextemp) = 10^(-0.02255*M+0.001434*M^2+(1.38*M-0.9254*M^2)/T); %water activity


        % combinatie van (101) en (95)
        
        U_rev(index_current,indextemp) = 1.229 - 0.9* 10^-3 * (T-298) +  R*T/(2*F) * log(pH2*sqrt(pO2)/aw(indextemp));
        
        % Paper (101)
        
        j_0c(indextemp) = gamma_Mcat*exp(-DeltaG/R*(1/T-1/Tref))*j_0cref; % getest en op deze manier verkrijg je logische waarden als je vergelijkt met gegeven range in (105)
        j_0a(indextemp) = gamma_Man*exp(-DeltaG/R*(1/T-1/Tref))*j_0aref;
        U_actc(index_current,indextemp) = R*T/(alfa_cat*2*F)*log(j/j_0c(indextemp));
        U_acta(index_current,indextemp) = R*T/(alfa_an*2*F)*log(j/j_0a(indextemp));
        U_act(index_current,indextemp) = U_acta(index_current,indextemp) + U_actc(index_current,indextemp);
        


        U_total(index_current,indextemp) = U_rev(index_current,indextemp) + U_act(index_current,indextemp) + U_ohm(index_current,indextemp) + U_conc(index_current,indextemp);
        Power(index_current,indextemp) = j*A*U_total(index_current,indextemp);

    end
end

if figures
    Tmin_celcius = Tmin-273; %°C
    Tmax_celcius = Tmax-273; %°C
    Trange_celcius = (Tmin_celcius):Tstep:(Tmax_celcius); %K
    jrange = jmin:(jmax-jmin)/((Tmax_celcius-Tmin_celcius)/Tstep):jmax; %A/m²
    
    figure(16)
    h = surf(Trange_celcius, jrange, Power) %Beetje vreemd surf functie plot jrange op y as en Trange op x as terwijl Power(index_current, index_temp)
    xlabel("Temperature [°C]", FontSize=12)
    ylabel("Current density [A/m²]",FontSize=12)
    zlabel("Power [W]",FontSize=10)
    set(h,'LineStyle','none')
    title("Power consumption per cell of an AWE",FontSize=14)
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
    print -depsc Power_AWE.eps
    
    counter_current = 1:25:length(Trange_celcius);
    counter_temp = 1:length(Trange_celcius);
    
    figure(17),hold on
    plot(Trange_celcius, U_total(counter_current,counter_temp), 'LineWidth', 1)
    axis tight
    xlabel('Temperature  [°C]')
    ylabel('U_{total} [Volt]')
    legend('Current density=' + string(jrange(counter_current(1))) + ' A/m²','Current density=' + string(jrange(counter_current(2))) + ' A/m²','Current density=' + string(jrange(counter_current(3))) + ' A/m²','Current density=' + string(jrange(counter_current(4))) + ' A/m²','Current density=' + string(jrange(counter_current(5))) + ' A/m²');
    hold off
    print -depsc U_tot_AWE.eps
    
    figure(18),hold on
    plot(Trange_celcius, U_rev(counter_current,counter_temp), 'LineWidth', 1)
    axis tight
    xlabel('Temperature  [°C]')
    ylabel('U_{rev} [Volt]')
    legend('Current density=' + string(jrange(counter_current(1))) + ' A/m²','Current density=' + string(jrange(counter_current(2))) + ' A/m²','Current density=' + string(jrange(counter_current(3))) + ' A/m²','Current density=' + string(jrange(counter_current(4))) + ' A/m²','Current density=' + string(jrange(counter_current(5))) + ' A/m²');
    hold off
    print -depsc U_rev_AWE.eps
    
    figure(19),hold on
    plot(Trange_celcius, U_act(counter_current,counter_temp), 'LineWidth', 1)
    axis tight
    xlabel('Temperature  [°C]')
    ylabel('U_{act} [Volt]')
    legend('Current density=' + string(jrange(counter_current(1))) + ' A/m²','Current density=' + string(jrange(counter_current(2))) + ' A/m²','Current density=' + string(jrange(counter_current(3))) + ' A/m²','Current density=' + string(jrange(counter_current(4))) + ' A/m²','Current density=' + string(jrange(counter_current(5))) + ' A/m²');
    hold off
    print -depsc U_act_AWE.eps
    
    figure(20),hold on
    plot(Trange_celcius, U_ohm(counter_current,counter_temp), 'LineWidth', 1)
    axis tight
    xlabel('Temperature  [°C]')
    ylabel('U_{ohm} [Volt]')
    legend('Current density=' + string(jrange(counter_current(1))) + ' A/m²','Current density=' + string(jrange(counter_current(2))) + ' A/m²','Current density=' + string(jrange(counter_current(3))) + ' A/m²','Current density=' + string(jrange(counter_current(4))) + ' A/m²','Current density=' + string(jrange(counter_current(5))) + ' A/m²');
    hold off
    print -depsc U_ohm_AWE.eps
    
    
    figure(21), hold on
    plot(Trange_celcius, U_conc(counter_current,counter_temp), 'LineWidth', 1)
    axis tight
    xlabel('Temperature  [°C]')
    ylabel('U_{conc} [Volt]')
    legend('Current density=' + string(jrange(counter_current(1))) + ' A/m²','Current density=' + string(jrange(counter_current(2))) + ' A/m²','Current density=' + string(jrange(counter_current(3))) + ' A/m²','Current density=' + string(jrange(counter_current(4))) + ' A/m²','Current density=' + string(jrange(counter_current(5))) + ' A/m²', 'Location', 'northwest');
    hold off
    print -depsc U_conc_AWE.eps
    
    counter_temp = 1:25:length(Trange_celcius);
    counter_current = 1:length(Trange_celcius);
    figure(22),hold on
    plot(jrange, U_total(counter_current,counter_temp), 'LineWidth', 1)
    axis tight
    xlabel('Current density [A/m²]')
    ylabel('U_{total} [Volt]')
    legend('Temperature =' + string(Trange_celcius(counter_temp(1))) + ' °C','Temperature=' + string(Trange_celcius(counter_temp(2))) + ' °C','Temperature=' + string(Trange_celcius(counter_temp(3))) + ' °C','Temperature=' + string(Trange_celcius(counter_temp(4))) + ' °C','Temperature=' + string(Trange_celcius(counter_temp(5))) + ' °C');
    hold off
    print -depsc U_tot.eps
    
    figure(23),hold on
    plot(jrange, U_rev(counter_current,counter_temp), 'LineWidth', 1)
    axis tight
    xlabel('Current density [A/m²]')
    ylabel('U_{rev} [Volt]')
    legend('Temperature =' + string(Trange_celcius(counter_temp(1))) + ' °C','Temperature=' + string(Trange_celcius(counter_temp(2))) + ' °C','Temperature=' + string(Trange_celcius(counter_temp(3))) + ' °C','Temperature=' + string(Trange_celcius(counter_temp(4))) + ' °C','Temperature=' + string(Trange_celcius(counter_temp(5))) + ' °C');
    hold off
    print -depsc U_rev.eps
    
    figure(24),hold on
    plot(jrange, U_act(counter_current,counter_temp), 'LineWidth', 1)
    axis tight
    xlabel('Current density [A/m²]')
    ylabel('U_{act} [Volt]')
    legend('Temperature =' + string(Trange_celcius(counter_temp(1))) + ' °C','Temperature=' + string(Trange_celcius(counter_temp(2))) + ' °C','Temperature=' + string(Trange_celcius(counter_temp(3))) + ' °C','Temperature=' + string(Trange_celcius(counter_temp(4))) + ' °C','Temperature=' + string(Trange_celcius(counter_temp(5))) + ' °C');
    hold off
    print -depsc U_act.eps
    
    figure(25),hold on
    plot(jrange, U_ohm(counter_current,counter_temp), 'LineWidth', 1)
    axis tight
    xlabel('Current density [A/m²]')
    ylabel('U_{ohm} [Volt]')
    legend('Temperature =' + string(Trange_celcius(counter_temp(1))) + ' °C','Temperature=' + string(Trange_celcius(counter_temp(2))) + ' °C','Temperature=' + string(Trange_celcius(counter_temp(3))) + ' °C','Temperature=' + string(Trange_celcius(counter_temp(4))) + ' °C','Temperature=' + string(Trange_celcius(counter_temp(5))) + ' °C');
    hold off
    print -depsc U_ohm.eps
    
    
    figure(26), hold on
    plot(jrange, U_conc(counter_current,counter_temp), 'LineWidth', 1)
    axis tight
    xlabel('Current density [A/m²]')
    ylabel('U_{conc} [Volt]')
    legend('Temperature =' + string(Trange_celcius(counter_temp(1))) + ' °C','Temperature=' + string(Trange_celcius(counter_temp(2))) + ' °C','Temperature=' + string(Trange_celcius(counter_temp(3))) + ' °C','Temperature=' + string(Trange_celcius(counter_temp(4))) + ' °C','Temperature=' + string(Trange_celcius(counter_temp(5))) + ' °C');
    hold off
    print -depsc U_conc.eps

%     IMatrix = zeros(length(Power));
%     TMatrix = zeros(length(Power));
%     for i = 1:length(Power)
%         for j = 1:length(Power)
%             IMatrix(i,j) = 0.21*jrange(i);
%             TMatrix(i,j) = Trange(j);
%         end
%     end
%     figure(27), hold on
%     Z = 6669*Power-1.4813*IMatrix*6669-(TMatrix-293)/0.0001133;
%     h = surf(Trange_celcius, jrange, Z);
%     set(h,'LineStyle','none')
%     title("Power*n_c - U_{tn}*I*n_c - (T - T_a)/R_t")
%     view(30,40)
%     grid on
%     xh = get(gca,'XLabel'); % Handle of the x label
%     set(xh, 'Units', 'Normalized')
%     pos = get(xh, 'Position');
%     set(xh, 'Position',pos.*[1,-0.0,1],'Rotation',-10)
%     yh = get(gca,'YLabel'); % Handle of the y label
%     set(yh, 'Units', 'Normalized')
%     pos = get(yh, 'Position');
%     set(yh, 'Position',pos.*[1.12,-0.7,1],'Rotation',40)
%    % Create a color map based on Z values
%     colormap([1 0 0; 0 1 0])
%     caxis([-1 1])
%     Z_color = Z >= 0; % Assign a color to each point based on the Z value
%     set(h, Z_color); % Assign the color map to the surface plot
end
