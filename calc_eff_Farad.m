function [eff_Farad] = calc_eff_Farad(jmin,jmax,Tmin,Tmax, Tstep, figures)
% Faraday efficiency for AWE: 5 parameter model
% Model from (124), parameters from (97)


B1 = 4.50424*10^(-5); %-
B2 = 1.02116; %-
B3 = -247.26; %A/m²
B4 = 2.06972; %A/(m²°C)
B5 = -0.03571; %A/(m²°C²)
Tmin_celcius = Tmin - 273; %°C
Tmax_celcius = Tmax -273; %°C
Trange_celcius = (Tmin_celcius):Tstep:(Tmax_celcius); %K
jrange = jmin:(jmax-jmin)/((Tmax_celcius-Tmin_celcius)/Tstep):jmax; %A/m²

eff_Farad = zeros(length(Trange_celcius), length(jrange));
indexcurrent = 0;
for j = jrange
    indexcurrent = indexcurrent +1;
    indextemp = 0;
    for T = Trange_celcius
        indextemp = indextemp + 1;
        eff_Farad(indexcurrent,indextemp) = B1 + B2*exp((B3+B4*T+B5*T^2)/j); 
    end 
end

if figures
    figure(30)
    h = surf(Trange_celcius, jrange, eff_Farad) %Ook hier weer jrange2 op x-as en Trange2 op y-as terwijl het wel eff_Farad(indexTrange, indexjrange) is
    xlabel("Temperature [°C]", FontSize=10)
    ylabel("Current density [A/m²]",FontSize=10)
    zlabel("Faraday efficiency [-]",FontSize=10)
    set(h,'LineStyle','none')
    title("Faraday efficiency for an AWE")
    view(30,40)
    grid on
    xh = get(gca,'XLabel'); % Handle of the x label
    set(xh, 'Units', 'Normalized')
    pos = get(xh, 'Position');
    set(xh, 'Position',pos.*[0.8,-0.0,1],'Rotation',-10)
    yh = get(gca,'YLabel'); % Handle of the y label
    set(yh, 'Units', 'Normalized')
    pos = get(yh, 'Position');
    set(yh, 'Position',pos.*[1.1,-0.7,1],'Rotation',40)
    print -depsc eff_Farad_LT.eps
end

