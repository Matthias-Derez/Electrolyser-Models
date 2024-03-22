function [Power] = calc_overpotentials_ASR(jmin,jmax,Tmin,Tmax, PowerSOEC, Area, Tstep, figures)
% Based on model described by (74) (ASR also used in 72)

Trange = Tmin:Tstep:Tmax; %K
jrange = jmin:(jmax-jmin)/((Tmax-Tmin)/Tstep):jmax; %A/m²

A = Area; %Voorlopig zelfde genomen als SOEC

Power = zeros(length(jrange), length(Trange));
for i = 1:length(Trange)
    for j = 1:length(jrange)
        Power(j,i) = (3-0.463+3.973*10^(-5)*exp(10300/Trange(i)))*10^-4*jrange(j)*jrange(j)*A;
    end 
end
if figures
    figure(41)
    h = surf(Trange-273, jrange, Power) %Beetje vreemd surf functie plot jrange op y as en Trange op x as terwijl Power(index_current, index_temp)
    xlabel("Temperature [°C]", FontSize=10)
    ylabel("Current density [A/m²]",FontSize=10)
    zlabel("Power [W]",FontSize=10)
    set(h,'LineStyle','none')
    title("Power consumption per cell of a SOEC electrolyser using the ASR empirical model")
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
    print -depsc Power.eps

    figure(42)
    h = surf(Trange-273, jrange, Power-PowerSOEC) %Beetje vreemd surf functie plot jrange op y as en Trange op x as terwijl Power(index_current, index_temp)
    xlabel("Temperature [°C]", FontSize=10)
    ylabel("Current density [A/m²]",FontSize=10)
    zlabel("Power [W]",FontSize=10)
    set(h,'LineStyle','none')
    title("Power ASR - Power via vgl")
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
    print -depsc Power.eps
    figure(42)
    h = surf(Trange-273, jrange, (PowerSOEC-Power)./PowerSOEC) %Beetje vreemd surf functie plot jrange op y as en Trange op x as terwijl Power(index_current, index_temp)
    xlabel("Temperature [°C]", FontSize=10)
    ylabel("Current density [A/m²]",FontSize=10)
    zlabel("Power [W]",FontSize=10)
    set(h,'LineStyle','none')
    title("Relative difference")
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
    print -depsc Power.eps
end