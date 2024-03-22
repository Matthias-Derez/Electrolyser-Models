%% SOEC: Calculating overpotentials
jminSOEC = 2000; %A/m²
jmaxSOEC = 10000; %A/m²
TminSOEC = 800+273; %K
TmaxSOEC = 1000+273; %K
[U_totalSOEC,U_revSOEC, U_actSOEC, U_ohmSOEC, U_concSOEC,PowerSOEC, j_0aSOEC, j_0cSOEC] = calc_overpotentials_SOEC(jminSOEC,jmaxSOEC,TminSOEC,TmaxSOEC);
TminSOEC_celcius = TminSOEC -273;
TmaxSOEC_celcius = TmaxSOEC - 273;

TrangeSOEC_celcius = TminSOEC_celcius:1:TmaxSOEC_celcius; %°C
jrangeSOEC = jminSOEC:(jmaxSOEC-jminSOEC)/(TmaxSOEC-TminSOEC):jmaxSOEC; %A/m²
% Power = calc_power(jrange,Trange);
figure(1)
h = surf(TrangeSOEC_celcius, jrangeSOEC, PowerSOEC) %Beetje vreemd surf functie plot jrange op y as en Trange op x as terwijl Power(index_current, index_temp)
xlabel("Temperature [°C]", FontSize=10)
ylabel("Current density [A/m²]",FontSize=10)
zlabel("Power [W]",FontSize=10)
set(h,'LineStyle','none')
title("Power consumption per cell of a SOEC electrolyser")
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

counter_current = 1:50:201;
counter_temp = 1:201;

figure(2),hold on
plot(TrangeSOEC_celcius, U_totalSOEC(counter_current,counter_temp), 'LineWidth', 1)
axis tight
xlabel('Temperature  [°C]')
ylabel('U_{total} [Volt]')
legend('Current density=' + string(jrangeSOEC(counter_current(1))) + ' A/m²','Current density=' + string(jrangeSOEC(counter_current(2))) + ' A/m²','Current density=' + string(jrangeSOEC(counter_current(3))) + ' A/m²','Current density=' + string(jrangeSOEC(counter_current(4))) + ' A/m²','Current density=' + string(jrangeSOEC(counter_current(5))) + ' A/m²');
hold off
print -depsc U_tot.eps

figure(3),hold on
plot(TrangeSOEC_celcius, U_revSOEC(counter_current,counter_temp), 'LineWidth', 1)
axis tight
xlabel('Temperature  [°C]')
ylabel('U_{rev} [Volt]')
legend('Current density=' + string(jrangeSOEC(counter_current(1))) + ' A/m²','Current density=' + string(jrangeSOEC(counter_current(2))) + ' A/m²','Current density=' + string(jrangeSOEC(counter_current(3))) + ' A/m²','Current density=' + string(jrangeSOEC(counter_current(4))) + ' A/m²','Current density=' + string(jrangeSOEC(counter_current(5))) + ' A/m²');
hold off
print -depsc U_rev.eps

figure(4),hold on
plot(TrangeSOEC_celcius, U_actSOEC(counter_current,counter_temp), 'LineWidth', 1)
axis tight
xlabel('Temperature  [°C]')
ylabel('U_{act} [Volt]')
legend('Current density=' + string(jrangeSOEC(counter_current(1))) + ' A/m²','Current density=' + string(jrangeSOEC(counter_current(2))) + ' A/m²','Current density=' + string(jrangeSOEC(counter_current(3))) + ' A/m²','Current density=' + string(jrangeSOEC(counter_current(4))) + ' A/m²','Current density=' + string(jrangeSOEC(counter_current(5))) + ' A/m²');
hold off
print -depsc U_act.eps

figure(5),hold on
plot(TrangeSOEC_celcius, U_ohmSOEC(counter_current,counter_temp), 'LineWidth', 1)
axis tight
xlabel('Temperature  [°C]')
ylabel('U_{ohm} [Volt]')
legend('Current density=' + string(jrangeSOEC(counter_current(1))) + ' A/m²','Current density=' + string(jrangeSOEC(counter_current(2))) + ' A/m²','Current density=' + string(jrangeSOEC(counter_current(3))) + ' A/m²','Current density=' + string(jrangeSOEC(counter_current(4))) + ' A/m²','Current density=' + string(jrangeSOEC(counter_current(5))) + ' A/m²');
hold off
print -depsc U_ohm.eps


figure(6), hold on
plot(TrangeSOEC_celcius, U_concSOEC(counter_current,counter_temp), 'LineWidth', 1)
axis tight
xlabel('Temperature  [°C]')
ylabel('U_{conc} [Volt]')
legend('Current density=' + string(jrangeSOEC(counter_current(1))) + ' A/m²','Current density=' + string(jrangeSOEC(counter_current(2))) + ' A/m²','Current density=' + string(jrangeSOEC(counter_current(3))) + ' A/m²','Current density=' + string(jrangeSOEC(counter_current(4))) + ' A/m²','Current density=' + string(jrangeSOEC(counter_current(5))) + ' A/m²', 'Location', 'northwest');
hold off
print -depsc U_conc.eps

counter_temp = 1:50:201;
counter_current = 1:201;
figure(7),hold on
plot(jrangeSOEC, U_totalSOEC(counter_current,counter_temp), 'LineWidth', 1)
axis tight
xlabel('Current density [A/m²]')
ylabel('U_{total} [Volt]')
legend('Temperature =' + string(TrangeSOEC_celcius(counter_temp(1))) + ' °C','Temperature=' + string(TrangeSOEC_celcius(counter_temp(2))) + ' °C','Temperature=' + string(TrangeSOEC_celcius(counter_temp(3))) + ' °C','Temperature=' + string(TrangeSOEC_celcius(counter_temp(4))) + ' °C','Temperature=' + string(TrangeSOEC_celcius(counter_temp(5))) + ' °C');
hold off
print -depsc U_tot.eps

figure(8),hold on
plot(jrangeSOEC, U_revSOEC(counter_current,counter_temp), 'LineWidth', 1)
axis tight
xlabel('Current density [A/m²]')
ylabel('U_{rev} [Volt]')
legend('Temperature =' + string(TrangeSOEC_celcius(counter_temp(1))) + ' °C','Temperature=' + string(TrangeSOEC_celcius(counter_temp(2))) + ' °C','Temperature=' + string(TrangeSOEC_celcius(counter_temp(3))) + ' °C','Temperature=' + string(TrangeSOEC_celcius(counter_temp(4))) + ' °C','Temperature=' + string(TrangeSOEC_celcius(counter_temp(5))) + ' °C');
hold off
print -depsc U_rev.eps

figure(9),hold on
plot(jrangeSOEC, U_actSOEC(counter_current,counter_temp), 'LineWidth', 1)
axis tight
xlabel('Current density [A/m²]')
ylabel('U_{act} [Volt]')
legend('Temperature =' + string(TrangeSOEC_celcius(counter_temp(1))) + ' °C','Temperature=' + string(TrangeSOEC_celcius(counter_temp(2))) + ' °C','Temperature=' + string(TrangeSOEC_celcius(counter_temp(3))) + ' °C','Temperature=' + string(TrangeSOEC_celcius(counter_temp(4))) + ' °C','Temperature=' + string(TrangeSOEC_celcius(counter_temp(5))) + ' °C');
hold off
print -depsc U_act.eps

figure(10),hold on
plot(jrangeSOEC, U_ohmSOEC(counter_current,counter_temp), 'LineWidth', 1)
axis tight
xlabel('Current density [A/m²]')
ylabel('U_{ohm} [Volt]')
legend('Temperature =' + string(TrangeSOEC_celcius(counter_temp(1))) + ' °C','Temperature=' + string(TrangeSOEC_celcius(counter_temp(2))) + ' °C','Temperature=' + string(TrangeSOEC_celcius(counter_temp(3))) + ' °C','Temperature=' + string(TrangeSOEC_celcius(counter_temp(4))) + ' °C','Temperature=' + string(TrangeSOEC_celcius(counter_temp(5))) + ' °C');
hold off
print -depsc U_ohm.eps


figure(11), hold on
plot(jrangeSOEC, U_concSOEC(counter_current,counter_temp), 'LineWidth', 1)
axis tight
xlabel('Current density [A/m²]')
ylabel('U_{conc} [Volt]')
legend('Temperature =' + string(TrangeSOEC_celcius(counter_temp(1))) + ' °C','Temperature=' + string(TrangeSOEC_celcius(counter_temp(2))) + ' °C','Temperature=' + string(TrangeSOEC_celcius(counter_temp(3))) + ' °C','Temperature=' + string(TrangeSOEC_celcius(counter_temp(4))) + ' °C','Temperature=' + string(TrangeSOEC_celcius(counter_temp(5))) + ' °C');
hold off
print -depsc U_conc.eps

