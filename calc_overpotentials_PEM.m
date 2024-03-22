function [U_total,U_rev, U_act, U_ohm, U_conc, Power] = calc_overpotentials_PEM(jmin,jmax,Tmin,Tmax, Area, Tstep, figures)
% Based on different papers


Trange = (Tmin):Tstep:(Tmax); %K
jrange = jmin:(jmax-jmin)/((Tmax-Tmin)/Tstep):jmax; %A/m²

R = 8.314; %J/(mol K)
F = 96485.3; %sA/mol
P = 1; %bar
pH2 = 1; %bar
pO2 = 1; %bar

A = Area; %m² voorlopig zelfde waarde gebruiken als voor SOEC

% Gegevens (42)
aw = 1; %activity between electrode and membrane is 1 for liquid water
% A = 0.016; %m²
d_membrane = 0.0254e-2; %m
d_electrode = 0.008e-2; %m
rho_electrode = 10.6e-8; %Ohm*m
epsilon = 0.3;
% lambda = 21;
rho_eff = rho_electrode/(1-epsilon)^1.5;
sigma_electrode = 1/rho_eff;
alfa = 1/2;
n = 2;



U_total = zeros(length(jrange), length(Trange));
U_rev = zeros(length(jrange), length(Trange));
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
%         Paper (66)
        lambda = 0.08533*T-6.77632;
%         lambda = 21;


        % Paper (42)
        U_rev(index_current,indextemp) = 1.229 - 0.9* 10^-3 * (T-298) +  R*T/(2*F) * log(pH2*sqrt(pO2)/aw); %aw = 1
        sigma_membrane = (0.005139*lambda - 0.00326)*exp(1268*(1/303-1/T))*100; %1/(Ohm*m) Maal 100 tov (42), heel belangrijk, alles daar in eenheden van cm dus sigma daar 1/(Ohm*cm) maar ik heb hier alles omgezet naar m!!
        U_ohm(index_current, indextemp) = j*(d_electrode/sigma_electrode+d_membrane/sigma_membrane); %A uit beide termen gehaald en via I naar j


        % Paper (113)
        j_0a(index_current, indextemp) = 1.08e-17*exp(0.086*T); %A/m²
        j_0c(index_current, indextemp) =  j_0a(index_current, indextemp)*10^4; %maal 10^4 zodat grootste waarde orde 10^3 heeft, 
                                                                                   % wat een acceptabele grootorde is als je kijkt naar (42) of (105) en de range die daar in tabel 6 wordt beschreven.
         
        U_act(index_current, indextemp) = R*T/(alfa*n*F)*log(j/j_0c(index_current, indextemp)) + R*T/(alfa*n*F)*log(j/j_0a(index_current, indextemp));
        U_conc(index_current, indextemp) = -4*R*T/(n*F)*log(1-j/(jmax*1.05)); 


        U_total(index_current, indextemp) = U_rev(index_current, indextemp) + U_conc(index_current, indextemp) + U_act(index_current, indextemp) + U_ohm(index_current, indextemp);
        Power(index_current,indextemp) = j*A*U_total(index_current,indextemp);


    end
end

if figures    
    Tmin_celcius = Tmin - 273; %°C
    Tmax_celcius = Tmax - 273; %°C
    Trange_celcius = (Tmin_celcius):Tstep:(Tmax_celcius); %K
    jrange = jmin:(jmax-jmin)/((Tmax_celcius-Tmin_celcius)/Tstep):jmax; %A/m²
    
    
    figure(35)
    h = surf(Trange_celcius, jrange, Power) %Beetje vreemd surf functie plot jrange op y as en Trange op x as terwijl Power(index_current, index_temp)
    xlabel("Temperature [°C]", FontSize=12)
    ylabel("Current density [A/m²]",FontSize=12)
    zlabel("Power [W]",FontSize=12)
    set(h,'LineStyle','none')
    title("Power consumption per cell of a PEM electrolyser",'FontSize', 14)
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
    set(yh, 'Position',pos.*[1.1,-0.7,1],'Rotation',38)
    print -depsc Power_PEM.eps
    
    % difference = Power-PowerUlle;
    % 
    % figure(36)
    % h = surf(Trange_celcius, jrange, difference) %Beetje vreemd surf functie plot jrange op y as en Trange op x as terwijl Power(index_current, index_temp)
    % xlabel("Temperature [°C]", FontSize=10)
    % ylabel("Current density [A/m²]",FontSize=10)
    % zlabel("Power [W]",FontSize=10)
    % set(h,'LineStyle','none')
    % title("Power - PowerUlle")
    % view(30,40)
    % grid on
    % xh = get(gca,'XLabel'); % Handle of the x label
    % set(xh, 'Units', 'Normalized')
    % pos = get(xh, 'Position');
    % set(xh, 'Position',pos.*[1,-0.0,1],'Rotation',-10)
    % yh = get(gca,'YLabel'); % Handle of the y label
    % set(yh, 'Units', 'Normalized')
    % pos = get(yh, 'Position');
    % set(yh, 'Position',pos.*[1.12,-0.7,1],'Rotation',40)
    % print -depsc Power_PEM.eps

%     IMatrix = zeros(length(Power));
%     TMatrix = zeros(length(Power));
%     for i = 1:length(Power)
%         for j = 1:length(Power)
%             IMatrix(i,j) = 0.21*jrange(i);
%             TMatrix(i,j) = Trange(j);
%         end
%     end
%     figure(27), hold on
%     Z = 1176*Power-1.4813*IMatrix*1176-(TMatrix-293)/0.00010666;
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
