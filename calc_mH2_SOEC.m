function [mH2, eff_total] = calc_mH2_SOEC(jmin, jmax, Tmin, Tmax, Power, Area, Tstep, figures)
% SOEC: Calculating mass rate of produced hydrogen

Tmin_celcius = Tmin - 273;
Tmax_celcius = Tmax - 273;
Trange_celcius = (Tmin_celcius):Tstep:(Tmax_celcius); %K
jrange = jmin:(jmax-jmin)/((Tmax_celcius-Tmin_celcius)/Tstep):jmax; %A/m²

% nc = 3800; % Number of cells
% A = 0.37; % m²
A = Area; %m² voorlopig deze waarde om gelijkaardig iets te hebben als bij SOEC
MH2 = 2.016*10^(-3); % Molar mass hydrogen kg/mol
F = 96485.3; % Faraday constant sA/mol
HHV_H2 = 141.88e6; %J/kg Higher heating value hydrogen
LHV_H2 = 119.96e6; %J/kg Lower heating value hydrogen

mH2 = zeros(length(Trange_celcius), length(jrange));
eff_total = zeros(length(Trange_celcius), length(jrange));
indexcurrent = 0;
for j = jrange
    indexcurrent = indexcurrent + 1;
    indextemp = 0;
    for T = Trange_celcius
        indextemp = indextemp + 1;
        mH2(indexcurrent,indextemp) = MH2*j*A/(2*F); %kg/s dit is per cell!!!
        eff_total(indexcurrent,indextemp) = mH2(indexcurrent,indextemp)/Power(indexcurrent,indextemp)*LHV_H2;
    end 
end

if figures
    figure(12)
    h = surf(Trange_celcius, jrange, mH2);%Ook hier weer jrange op x-as en Trange_celcius op y-as terwijl het wel mH2_(indexTrange, indexjrange) is
    xlabel("Temperature [°C]", FontSize=12)
    ylabel("Current density [A/m²]",FontSize=12)
    zlabel("Hydrogen production rate [kg/s]",FontSize=12)
    set(h,'LineStyle','none')
    title("Hydrogen production rate for a SOEC per cell",FontSize=14)
    view(30,40)
    grid on
    xh = get(gca,'XLabel'); % Handle of the x label
    set(xh, 'Units', 'Normalized')
    pos = get(xh, 'Position');
    set(xh, 'Position',pos.*[0.8,-0.0,1],'Rotation',-10)
    yh = get(gca,'YLabel'); % Handle of the y label
    set(yh, 'Units', 'Normalized')
    pos = get(yh, 'Position');
    set(yh, 'Position',pos.*[1,-0.7,1],'Rotation',40)
    print -depsc prod_rate_H2_LT.eps
    
    
    figure(13)
    h = surf(Trange_celcius, jrange, eff_total); %Ook hier weer jrangeSOEC op y-as en TrangeSOEC_celcius op x-as terwijl het wel mH2_SOEC(indexTrange, indexjrange) is
    xlabel("Temperature [°C]", FontSize=12)
    ylabel("Current density [A/m²]",FontSize=12)
    zlabel("Efficiency [-]",FontSize=12)
    set(h,'LineStyle','none')
    title("Efficiency SOEC",FontSize=16)
    view(30,40)
    grid on
    xh = get(gca,'XLabel'); % Handle of the x label
    set(xh, 'Units', 'Normalized')
    pos = get(xh, 'Position');
    set(xh, 'Position',pos.*[0.9,-0.0,1],'Rotation',-10)
    yh = get(gca,'YLabel'); % Handle of the y label
    set(yh, 'Units', 'Normalized')
    pos = get(yh, 'Position');
    set(yh, 'Position',pos.*[1.1,-0.7,1],'Rotation',40)
    print -depsc efficiency_SOEC.eps

    counter_temp = length(Trange_celcius):-100:1;
    counter_current = 1:length(Trange_celcius);
    figure(14),hold on
    plot(jrange, eff_total(counter_current,counter_temp), 'LineWidth', 1)
    axis tight
    ylim([0 1.1])
    xlabel('Current density [A/m²]', FontSize=12)
    ylabel('Efficiency [-]', FontSize=12)
    legend('Temperature = ' + string(Trange_celcius(counter_temp(1))) + '°C','Temperature = ' + string(Trange_celcius(counter_temp(2))) + '°C','Temperature = ' + string(Trange_celcius(counter_temp(3))) + '°C','Temperature = ' + string(Trange_celcius(counter_temp(4))) + '°C','Temperature = ' + string(Trange_celcius(counter_temp(5))) + '°C', Location="southwest");
    title("Efficiency of an SOE", FontSize=14)
    hold off


%     figure(13)
%     h = surf(Trange_celcius, jrange, eff_total); %Ook hier weer jrangeSOEC op y-as en TrangeSOEC_celcius op x-as terwijl het wel mH2_SOEC(indexTrange, indexjrange) is
%     ylabel("Temperature [°C]", FontSize=10)
%     xlabel("Current density [A/m²]",FontSize=10)
%     zlabel("Efficiency [-]",FontSize=10)
%     set(h,'LineStyle','none')
%     title("Efficiency SOEC")
%     view(30,40)
%     grid on
%     xh = get(gca,'XLabel'); % Handle of the x label
%     set(xh, 'Units', 'Normalized')
%     pos = get(xh, 'Position');
%     set(xh, 'Position',pos.*[0.5,-0.2,1])
%     yh = get(gca,'YLabel'); % Handle of the y label
%     set(yh, 'Units', 'Normalized')
%     pos = get(yh, 'Position');
%     set(yh, 'Position',pos.*[1,-0.1,1])
%     print -depsc efficiency_SOEC.eps
end
