function [U_total, U_rev, Power_density, Power] = calc_overpotentials_Ulleberg(jmin,jmax,Tmin,Tmax, PowerAWE, Area, Tstep, figures)
% Based on model described by Ulleberg "Modeling of advanced alkaline
% electrolyzers: a system simulation approach" (paper 67)

Trange = Tmin:Tstep:Tmax; %K
jrange = jmin:(jmax-jmin)/((Tmax-Tmin)/Tstep):jmax; %A/m²
A = Area; %Voorlopig zelfde genomen als SOEC

% Ulleberg
r1 = 8.05e-5; %Ohm*m²
r2 = -2.5e-7; %Ohm*m²/°C
s = 0.185; %V
% t1 = -1.002; %m²/A waarde die in (67) staat
t1= -0.1002; % wat jij & alexander denken dat ulleberg verkeerd heeft gedaan als je naar die grafiek kijkt
t2 = 8.424; %m²°C/A
t3 = 247.3; %m²°C²/A
% Mooi verwoord: values of the constants r1, r2, t1, t2, and t3 in the Ulleberg model can vary depending on the specific type of electrolyzer and operating conditions.
%  it's important to note that these values are not universal and can vary depending on the specific electrolyzer and operating conditions.
% Zheng
% r1 = 8.05e-5; %Ohm*m²
% r2 = -2.5e-7; %Ohm*m²/°C
% s = 0.08; %V
% t1 = 1.002; %m²/A
% t2 = 8.424; %m²°C/A
% t3 = 247.3e3; %m²°C²/A
% paper (112)
% r1 = 3.53855e-4; %Ohm*m²
% r2 = -3.0215e-6; %Ohm*m²/°C
% s = 2.2396e-1; %V
% t1 = 5.13093; %m²/A
% t2 = -2.40447e2; %m²°C/A
% t3 = 5.99576e3; %m²°C²/A


   
U_total = zeros(length(jrange), length(Trange));
U_rev = zeros(length(jrange), length(Trange));
Power = zeros(length(jrange), length(Trange));
Power_density = zeros(length(jrange), length(Trange));
index_current = 0;
for j = jrange
    index_current = index_current + 1;
    indextemp = 0;
    for T = Trange
        indextemp = indextemp + 1;
        % Zheng
%         U_rev(index_current, indextemp) = 1.5184 - 1.5421e-3*T + 9.523e-5*T*log(T) + 9.84e-8*T^2; 
%         U_total(index_current, indextemp) = U_rev(index_current, indextemp) + (r1+r2*(T-273))*j + s*log((t1 + t2/(T-273) + t3/(T-273)^2)*j+1);
%         Power_density(index_current, indextemp) = U_total(index_current, indextemp)*j; %W/m²
%         Power(index_current, indextemp) = U_total(index_current, indextemp)*j*A;
        % Ulleberg
        U_rev(index_current, indextemp) = 1.5184 - 1.5421e-3*T + 9.523e-5*T*log(T) + 9.84e-8*T^2; 
        U_total(index_current, indextemp) = U_rev(index_current, indextemp) + (r1+r2*(T-273))*j + s*log((t1 + t2/(T-273) + t3/(T-273)^2)*j+1);
        Power_density(index_current, indextemp) = U_total(index_current, indextemp)*j; %W/m²
        Power(index_current, indextemp) = U_total(index_current, indextemp)*j*A;
    
    end
end
if figures

    Tmin_celcius = Tmin - 273; %°C
    Tmax_celcius = Tmax - 273; %°C
    Trange_celcius = (Tmin_celcius):Tstep:(Tmax_celcius); %K
    jrange = jmin:(jmax-jmin)/((Tmax_celcius-Tmin_celcius)/Tstep):jmax; %A/m²
    
    figure(27)
    h = surf(Trange_celcius, jrange, Power_density) %Beetje vreemd surf functie plot jrange op y as en Trange op x as terwijl Power(index_current, index_temp)
    xlabel("Temperature [°C]", FontSize=10)
    ylabel("Current density [A/m²]",FontSize=10)
    zlabel("Power density [W/m²]",FontSize=10)
    set(h,'LineStyle','none')
    title("Power density per cell of an AWE electrolyser using ulleberg")
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
    
    figure(28)
    h = surf(Trange_celcius, jrange, Power) %Beetje vreemd surf functie plot jrange op y as en Trange op x as terwijl Power(index_current, index_temp)
    xlabel("Temperature [°C]", FontSize=10)
    ylabel("Current density [A/m²]",FontSize=10)
    zlabel("Power [W]",FontSize=10)
    set(h,'LineStyle','none')
    title("Power per cell of an AWE electrolyser using ulleberg")
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
    
    difference = (PowerAWE-Power)./PowerAWE;
    
    figure(29)
    h = surf(Trange_celcius, jrange, difference) %Beetje vreemd surf functie plot jrange op y as en Trange op x as terwijl Power(index_current, index_temp)
    xlabel("Temperature [°C]", FontSize=10)
    ylabel("Current density [A/m²]",FontSize=10)
    zlabel("Relative difference [-]",FontSize=10)
    set(h,'LineStyle','none')
    title({"Relative difference in consumed power","between theoretical and Ulleberg model"})
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
end
