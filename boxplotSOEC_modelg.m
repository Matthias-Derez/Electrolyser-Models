clearvars
format long
%% Model 1g
objectives_SOEC_model1g = load("objectives_SOEC_model1g.mat");

objectives_SOEC_model1g_15 = load("objectives_SOEC_model1g_15.mat");

obj1g_1_5 = [objectives_SOEC_model1g_15.objective_values(1:158,1);objectives_SOEC_model1g_15.objective_values(160:365,1)];
obj1g_2_5 = [objectives_SOEC_model1g.objective_values(1:158,1);objectives_SOEC_model1g.objective_values(160:365,1)];
obj1g_3_5 = [objectives_SOEC_model1g.objective_values(1:158,2);objectives_SOEC_model1g.objective_values(160:365,2)];
obj1g_4_5 = [objectives_SOEC_model1g.objective_values(1:158,3);objectives_SOEC_model1g.objective_values(160:365,3)];
obj1g_5_5 = [objectives_SOEC_model1g.objective_values(1:158,4);objectives_SOEC_model1g.objective_values(160:365,4)];


figure
boxplot(obj1g_1_5);  % Create the boxplot
title({"Revenues per day for an SOEC model1g in 2019"},{"with hydrogen price = 1.5€/kg"})
hold on
plot(mean(obj1g_1_5), 'xg')
hold off
print -depsc boxplot_SOEC5e_1.5.eps


figure
boxplot(obj1g_2_5);  % Create the boxplot
title({"Revenues per day for an SOEC model1g in 2019"},{"with hydrogen price = 2.5€/kg"})
hold on
plot(mean(obj1g_2_5), 'xg')
hold off
print -depsc boxplot_SOEC5e_2.5.eps

figure
boxplot(obj1g_3_5);  % Create the boxplot
title({"Revenues per day for an SOEC model1g in 2019"},{"with hydrogen price = 3.5€/kg"})
hold on
plot(mean(obj1g_3_5), 'xg')
hold off
print -depsc boxplot_SOEC1g_3.5.eps

figure
boxplot(obj1g_4_5);  % Create the boxplot
title({"Revenues per day for an SOEC model1g in 2019"},{"with hydrogen price = 4.5€/kg"})
hold on
plot(mean(obj1g_4_5), 'xg')
hold off
print -depsc boxplot_SOEC1g_4.5.eps

figure
boxplot(obj1g_5_5);  % Create the boxplot
title({"Revenues per day for an SOEC model5e in 2019"},{"with hydrogen price = 5.5€/kg"})
hold on
plot(mean(obj1g_5_5), 'xg')
hold off
print -depsc boxplot_SOEC1g_5.5.eps

mean_rev_1g_1_5 = mean(obj1g_1_5)
mean_rev_1g_2_5 = mean(obj1g_2_5)
mean_rev_1g_3_5 = mean(obj1g_3_5)
mean_rev_1g_4_5 = mean(obj1g_4_5)
mean_rev_1g_5_5 = mean(obj1g_5_5)






HHV_H2 = 141.88e6; %J/kg Higher heating value hydrogen
LHV_H2 = 119.96e6; %J/kg Lower heating value hydrogen

indices = find(abs(objectives_SOEC_model1g_15.objective_values(:,5)) > 10^-6);
mH2_1g_1_5 = objectives_SOEC_model1g_15.objective_values(indices,5);
p_E_1g_1_5 = objectives_SOEC_model1g_15.objective_values(indices,9);

indices = find(abs(objectives_SOEC_model1g.objective_values(:,5)) > 10^-6);
mH2_1g_2_5 = objectives_SOEC_model1g.objective_values(indices,5);
p_E_1g_2_5 = objectives_SOEC_model1g.objective_values(indices,9);

indices = find(abs(objectives_SOEC_model1g.objective_values(:,6)) > 10^-6);
mH2_1g_3_5 = objectives_SOEC_model1g.objective_values(indices,6);
p_E_1g_3_5 = objectives_SOEC_model1g.objective_values(indices,10);

indices = find(abs(objectives_SOEC_model1g.objective_values(:,7)) > 10^-6);
mH2_1g_4_5 = objectives_SOEC_model1g.objective_values(indices,7);
p_E_1g_4_5 = objectives_SOEC_model1g.objective_values(indices,11);


indices = find(abs(objectives_SOEC_model1g.objective_values(:,8)) > 10^-6);
mH2_1g_5_5 = objectives_SOEC_model1g.objective_values(indices,8);
p_E_1g_5_5 = objectives_SOEC_model1g.objective_values(indices,12);


eff_1g_1_5 = zeros(length(mH2_1g_1_5),1);
eff_1g_2_5 = zeros(length(mH2_1g_2_5),1);
eff_1g_3_5 = zeros(length(mH2_1g_3_5),1);
eff_1g_4_5 = zeros(length(mH2_1g_4_5),1);
eff_1g_5_5 = zeros(length(mH2_1g_5_5),1);
for i = 1:length(mH2_1g_1_5)
    eff_1g_1_5(i) =  mH2_1g_1_5(i)*LHV_H2/p_E_1g_1_5(i);
end
for i = 1:length(mH2_1g_2_5)
    eff_1g_2_5(i) =  mH2_1g_2_5(i)*LHV_H2/p_E_1g_2_5(i);
end
for i = 1:length(mH2_1g_3_5)
    eff_1g_3_5(i) =  mH2_1g_3_5(i)*LHV_H2/p_E_1g_3_5(i);
end
for i = 1:length(mH2_1g_4_5)
    eff_1g_4_5(i) =  mH2_1g_4_5(i)*LHV_H2/p_E_1g_4_5(i);
end
for i = 1:length(mH2_1g_5_5)
    eff_1g_5_5(i) =  mH2_1g_5_5(i)*LHV_H2/p_E_1g_5_5(i);
end

% for i = 1:8760
% %     eff_1g_2_5(i) =  mH2_1g_2_5(i)*LHV_H2/p_E_1g_2_5(i);
% %     eff_1g_3_5(i) =  mH2_1g_3_5(i)*LHV_H2/p_E_1g_3_5(i);
%     eff_1g_4_5(i) =  mH2_1g_4_5(i)*LHV_H2/p_E_1g_4_5(i);
%     eff_1g_5_5(i) =  mH2_1g_5_5(i)*LHV_H2/p_E_1g_5_5(i);
% end
eff_1g_1_5 = eff_1g_1_5(~isnan(eff_1g_1_5));
eff_1g_2_5 = eff_1g_2_5(~isnan(eff_1g_2_5));
eff_1g_3_5 = eff_1g_3_5(~isnan(eff_1g_3_5));
eff_1g_4_5 = eff_1g_4_5(~isnan(eff_1g_4_5));
eff_1g_5_5 = eff_1g_5_5(~isnan(eff_1g_5_5));
lower_limit = -1e-6;
upper_limit = 10;

mean_eff_1g_1_5 = mean(eff_1g_1_5(eff_1g_1_5 > lower_limit & eff_1g_1_5 < upper_limit));
mean_eff_1g_2_5 = mean(eff_1g_2_5(eff_1g_2_5 > lower_limit & eff_1g_2_5 < upper_limit));
mean_eff_1g_3_5 = mean(eff_1g_3_5(eff_1g_3_5 > lower_limit & eff_1g_3_5 < upper_limit));
mean_eff_1g_4_5 = mean(eff_1g_4_5(eff_1g_4_5 > lower_limit & eff_1g_4_5 < upper_limit));
mean_eff_1g_5_5 = mean(eff_1g_5_5(eff_1g_5_5 > lower_limit & eff_1g_5_5 < upper_limit));
%% Model 2g
objectives_SOEC_model2g = load("objectives_SOEC_model2g.mat");
objectives_SOEC_model2g_15 = load("objectives_SOEC_model2g_15.mat");

obj2g_1_5 = [objectives_SOEC_model2g_15.objective_values(1:158,1);objectives_SOEC_model2g_15.objective_values(160:365,1)];
obj2g_2_5 = [objectives_SOEC_model2g.objective_values(1:158,1);objectives_SOEC_model2g.objective_values(160:365,1)];
obj2g_3_5 = [objectives_SOEC_model2g.objective_values(1:158,2);objectives_SOEC_model2g.objective_values(160:365,2)];
obj2g_4_5 = [objectives_SOEC_model2g.objective_values(1:158,3);objectives_SOEC_model2g.objective_values(160:365,3)];
obj2g_5_5 = [objectives_SOEC_model2g.objective_values(1:158,4);objectives_SOEC_model2g.objective_values(160:365,4)];

figure
boxplot(obj2g_1_5);  % Create the boxplot
title({"Revenues per day for an SOEC model2g in 2019"},{"with hydrogen price = 1.5€/kg"})
hold on
plot(mean(obj2g_1_5), 'xg')
hold off
print -depsc boxplot_SOEC5e_1.5.eps

figure
boxplot(obj2g_2_5);  % Create the boxplot
title({"Revenues per day for an SOEC model2g in 2019"},{"with hydrogen price = 2.5€/kg"})
hold on
plot(mean(obj2g_2_5), 'xg')
hold off
print -depsc boxplot_SOEC5e_2.5.eps

figure
boxplot(obj2g_3_5);  % Create the boxplot
title({"Revenues per day for an SOEC model2g in 2019"},{"with hydrogen price = 3.5€/kg"})
hold on
plot(mean(obj2g_3_5), 'xg')
hold off
print -depsc boxplot_SOEC2g_3.5.eps

figure
boxplot(obj2g_4_5);  % Create the boxplot
title({"Revenues per day for an SOEC model2g in 2019"},{"with hydrogen price = 4.5€/kg"})
hold on
plot(mean(obj2g_4_5), 'xg')
hold off
print -depsc boxplot_SOEC2g_4.5.eps

figure
boxplot(obj2g_5_5);  % Create the boxplot
title({"Revenues per day for an SOEC model5e in 2019"},{"with hydrogen price = 5.5€/kg"})
hold on
plot(mean(obj2g_5_5), 'xg')
hold off
print -depsc boxplot_SOEC2g_5.5.eps

mean_rev_2g_1_5 = mean(obj2g_1_5)
mean_rev_2g_2_5 = mean(obj2g_2_5)
mean_rev_2g_3_5 = mean(obj2g_3_5)
mean_rev_2g_4_5 = mean(obj2g_4_5)
mean_rev_2g_5_5 = mean(obj2g_5_5)





HHV_H2 = 141.88e6; %J/kg Higher heating value hydrogen
LHV_H2 = 119.96e6; %J/kg Lower heating value hydrogen

indices = find(abs(objectives_SOEC_model2g_15.objective_values(:,5)) > 10^-6);
mH2_2g_1_5 = objectives_SOEC_model2g_15.objective_values(indices,5);
p_E_2g_1_5 = objectives_SOEC_model2g_15.objective_values(indices,9);

indices = find(abs(objectives_SOEC_model2g.objective_values(:,5)) > 10^-6);
mH2_2g_2_5 = objectives_SOEC_model2g.objective_values(indices,5);
p_E_2g_2_5 = objectives_SOEC_model2g.objective_values(indices,9);

indices = find(abs(objectives_SOEC_model2g.objective_values(:,6)) > 10^-6);
mH2_2g_3_5 = objectives_SOEC_model2g.objective_values(indices,6);
p_E_2g_3_5 = objectives_SOEC_model2g.objective_values(indices,10);

indices = find(abs(objectives_SOEC_model2g.objective_values(:,7)) > 10^-6);
mH2_2g_4_5 = objectives_SOEC_model2g.objective_values(indices,7);
p_E_2g_4_5 = objectives_SOEC_model2g.objective_values(indices,11);


indices = find(abs(objectives_SOEC_model2g.objective_values(:,8)) > 10^-6);
mH2_2g_5_5 = objectives_SOEC_model2g.objective_values(indices,8);
p_E_2g_5_5 = objectives_SOEC_model2g.objective_values(indices,12);


eff_2g_1_5 = zeros(length(mH2_2g_1_5),1);
eff_2g_2_5 = zeros(length(mH2_2g_2_5),1);
eff_2g_3_5 = zeros(length(mH2_2g_3_5),1);
eff_2g_4_5 = zeros(length(mH2_2g_4_5),1);
eff_2g_5_5 = zeros(length(mH2_2g_5_5),1);
for i = 1:length(mH2_2g_1_5)
    eff_2g_1_5(i) =  mH2_2g_1_5(i)*LHV_H2/p_E_2g_1_5(i);
end
for i = 1:length(mH2_2g_2_5)
    eff_2g_2_5(i) =  mH2_2g_2_5(i)*LHV_H2/p_E_2g_2_5(i);
end
for i = 1:length(mH2_2g_3_5)
    eff_2g_3_5(i) =  mH2_2g_3_5(i)*LHV_H2/p_E_2g_3_5(i);
end
for i = 1:length(mH2_2g_4_5)
    eff_2g_4_5(i) =  mH2_2g_4_5(i)*LHV_H2/p_E_2g_4_5(i);
end
for i = 1:length(mH2_2g_5_5)
    eff_2g_5_5(i) =  mH2_2g_5_5(i)*LHV_H2/p_E_2g_5_5(i);
end

% for i = 1:8760
% %     eff_2g_2_5(i) =  mH2_2g_2_5(i)*LHV_H2/p_E_2g_2_5(i);
% %     eff_2g_3_5(i) =  mH2_2g_3_5(i)*LHV_H2/p_E_2g_3_5(i);
%     eff_2g_4_5(i) =  mH2_2g_4_5(i)*LHV_H2/p_E_2g_4_5(i);
%     eff_2g_5_5(i) =  mH2_2g_5_5(i)*LHV_H2/p_E_2g_5_5(i);
% end
eff_2g_1_5 = eff_2g_1_5(~isnan(eff_2g_1_5));
eff_2g_2_5 = eff_2g_2_5(~isnan(eff_2g_2_5));
eff_2g_3_5 = eff_2g_3_5(~isnan(eff_2g_3_5));
eff_2g_4_5 = eff_2g_4_5(~isnan(eff_2g_4_5));
eff_2g_5_5 = eff_2g_5_5(~isnan(eff_2g_5_5));
lower_limit = -1e-6;
upper_limit = 10;

mean_eff_2g_1_5 = mean(eff_2g_1_5(eff_2g_1_5 > lower_limit & eff_2g_1_5 < upper_limit));
mean_eff_2g_2_5 = mean(eff_2g_2_5(eff_2g_2_5 > lower_limit & eff_2g_2_5 < upper_limit));
mean_eff_2g_3_5 = mean(eff_2g_3_5(eff_2g_3_5 > lower_limit & eff_2g_3_5 < upper_limit));
mean_eff_2g_4_5 = mean(eff_2g_4_5(eff_2g_4_5 > lower_limit & eff_2g_4_5 < upper_limit));
mean_eff_2g_5_5 = mean(eff_2g_5_5(eff_2g_5_5 > lower_limit & eff_2g_5_5 < upper_limit));
revenue_week20_perday = sum(obj2g_2_5(141:147))
revenue_week20_perweek = 72713.63 % nog niet klaar met runnen incumbent al lang 72713.63, best bound 72741.3
rel_diff = (revenue_week20_perweek-revenue_week20_perday)/revenue_week20_perday

revenue_week20_perday = sum(obj2g_5_5(141:147))
revenue_week20_perweek = 302634.77029693074 % wel heel rap klaar met runnen 
rel_diff = (revenue_week20_perweek-revenue_week20_perday)/revenue_week20_perday


%% Model 5g
objectives_SOEC_model5g = load("objectives_SOEC_model5g.mat");
objectives_SOEC_model5g_15 = load("objectives_SOEC_model5g_15.mat");

obj5g_1_5 = [objectives_SOEC_model5g_15.objective_values(1:158,1);objectives_SOEC_model5g_15.objective_values(160:365,1)];

obj5g_2_5 = [objectives_SOEC_model5g.objective_values(1:158,1);objectives_SOEC_model5g.objective_values(160:365,1)];
obj5g_3_5 = [objectives_SOEC_model5g.objective_values(1:158,2);objectives_SOEC_model5g.objective_values(160:365,2)];
obj5g_4_5 = [objectives_SOEC_model5g.objective_values(1:158,3);objectives_SOEC_model5g.objective_values(160:365,3)];
obj5g_5_5 = [objectives_SOEC_model5g.objective_values(1:158,4);objectives_SOEC_model5g.objective_values(160:365,4)];


figure
boxplot(obj5g_1_5);  % Create the boxplot
title({"Revenues per day for an SOEC model5g in 2019"},{"with hydrogen price = 1.5€/kg"})
hold on
plot(mean(obj5g_1_5), 'xg')
hold off
print -depsc boxplot_SOEC5e_1.5.eps


figure
boxplot(obj5g_2_5);  % Create the boxplot
title({"Revenues per day for an SOEC model5g in 2019"},{"with hydrogen price = 2.5€/kg"})
hold on
plot(mean(obj5g_2_5), 'xg')
hold off
print -depsc boxplot_SOEC5e_2.5.eps

figure
boxplot(obj5g_3_5);  % Create the boxplot
title({"Revenues per day for an SOEC model5g in 2019"},{"with hydrogen price = 3.5€/kg"})
hold on
plot(mean(obj5g_3_5), 'xg')
hold off
print -depsc boxplot_SOEC5g_3.5.eps

figure
boxplot(obj5g_4_5);  % Create the boxplot
title({"Revenues per day for an SOEC model5g in 2019"},{"with hydrogen price = 4.5€/kg"})
hold on
plot(mean(obj5g_4_5), 'xg')
hold off
print -depsc boxplot_SOEC5g_4.5.eps

figure
boxplot(obj5g_5_5);  % Create the boxplot
title({"Revenues per day for an SOEC model5e in 2019"},{"with hydrogen price = 5.5€/kg"})
hold on
plot(mean(obj5g_5_5), 'xg')
hold off
print -depsc boxplot_SOEC5g_5.5.eps

mean_rev_5g_1_5 = mean(obj5g_1_5)
mean_rev_5g_2_5 = mean(obj5g_2_5)
mean_rev_5g_3_5 = mean(obj5g_3_5)
mean_rev_5g_4_5 = mean(obj5g_4_5)
mean_rev_5g_5_5 = mean(obj5g_5_5)






HHV_H2 = 141.88e6; %J/kg Higher heating value hydrogen
LHV_H2 = 119.96e6; %J/kg Lower heating value hydrogen

indices = find(abs(objectives_SOEC_model5g_15.objective_values(:,5)) > 10^-6);
mH2_5g_1_5 = objectives_SOEC_model5g_15.objective_values(indices,5);
p_E_5g_1_5 = objectives_SOEC_model5g_15.objective_values(indices,9);

indices = find(abs(objectives_SOEC_model5g.objective_values(:,5)) > 10^-6);
mH2_5g_2_5 = objectives_SOEC_model5g.objective_values(indices,5);
p_E_5g_2_5 = objectives_SOEC_model5g.objective_values(indices,9);

indices = find(abs(objectives_SOEC_model5g.objective_values(:,6)) > 10^-6);
mH2_5g_3_5 = objectives_SOEC_model5g.objective_values(indices,6);
p_E_5g_3_5 = objectives_SOEC_model5g.objective_values(indices,10);

indices = find(abs(objectives_SOEC_model5g.objective_values(:,7)) > 10^-6);
mH2_5g_4_5 = objectives_SOEC_model5g.objective_values(indices,7);
p_E_5g_4_5 = objectives_SOEC_model5g.objective_values(indices,11);


indices = find(abs(objectives_SOEC_model5g.objective_values(:,8)) > 10^-6);
mH2_5g_5_5 = objectives_SOEC_model5g.objective_values(indices,8);
p_E_5g_5_5 = objectives_SOEC_model5g.objective_values(indices,12);


eff_5g_1_5 = zeros(length(mH2_5g_1_5),1);
eff_5g_2_5 = zeros(length(mH2_5g_2_5),1);
eff_5g_3_5 = zeros(length(mH2_5g_3_5),1);
eff_5g_4_5 = zeros(length(mH2_5g_4_5),1);
eff_5g_5_5 = zeros(length(mH2_5g_5_5),1);
for i = 1:length(mH2_5g_1_5)
    eff_5g_1_5(i) =  mH2_5g_1_5(i)*LHV_H2/p_E_5g_1_5(i);
end
for i = 1:length(mH2_5g_2_5)
    eff_5g_2_5(i) =  mH2_5g_2_5(i)*LHV_H2/p_E_5g_2_5(i);
end
for i = 1:length(mH2_5g_3_5)
    eff_5g_3_5(i) =  mH2_5g_3_5(i)*LHV_H2/p_E_5g_3_5(i);
end
for i = 1:length(mH2_5g_4_5)
    eff_5g_4_5(i) =  mH2_5g_4_5(i)*LHV_H2/p_E_5g_4_5(i);
end
for i = 1:length(mH2_5g_5_5)
    eff_5g_5_5(i) =  mH2_5g_5_5(i)*LHV_H2/p_E_5g_5_5(i);
end

% for i = 1:8760
% %     eff_5g_2_5(i) =  mH2_5g_2_5(i)*LHV_H2/p_E_5g_2_5(i);
% %     eff_5g_3_5(i) =  mH2_5g_3_5(i)*LHV_H2/p_E_5g_3_5(i);
%     eff_5g_4_5(i) =  mH2_5g_4_5(i)*LHV_H2/p_E_5g_4_5(i);
%     eff_5g_5_5(i) =  mH2_5g_5_5(i)*LHV_H2/p_E_5g_5_5(i);
% end
eff_5g_1_5 = eff_5g_1_5(~isnan(eff_5g_1_5));
eff_5g_2_5 = eff_5g_2_5(~isnan(eff_5g_2_5));
eff_5g_3_5 = eff_5g_3_5(~isnan(eff_5g_3_5));
eff_5g_4_5 = eff_5g_4_5(~isnan(eff_5g_4_5));
eff_5g_5_5 = eff_5g_5_5(~isnan(eff_5g_5_5));
lower_limit = -1e-6;
upper_limit = 10;

mean_eff_5g_1_5 = mean(eff_5g_1_5(eff_5g_1_5 > lower_limit & eff_5g_1_5 < upper_limit));
mean_eff_5g_2_5 = mean(eff_5g_2_5(eff_5g_2_5 > lower_limit & eff_5g_2_5 < upper_limit));
mean_eff_5g_3_5 = mean(eff_5g_3_5(eff_5g_3_5 > lower_limit & eff_5g_3_5 < upper_limit));
mean_eff_5g_4_5 = mean(eff_5g_4_5(eff_5g_4_5 > lower_limit & eff_5g_4_5 < upper_limit));
mean_eff_5g_5_5 = mean(eff_5g_5_5(eff_5g_5_5 > lower_limit & eff_5g_5_5 < upper_limit));
revenue_week20_perday = sum(obj5g_2_5(141:147))
revenue_week20_perweek = 83506.65 % nog niet klaar met runnen incumbent al lang 72713.63, best bound 72741.3
rel_diff = (revenue_week20_perweek-revenue_week20_perday)/revenue_week20_perday

revenue_week20_perday = sum(obj5g_5_5(141:147))
revenue_week20_perweek = 313427.78 % wel heel rap klaar met runnen 
rel_diff = (revenue_week20_perweek-revenue_week20_perday)/revenue_week20_perday


%% Model 6g
objectives_SOEC_model6g = load("objectives_SOEC_model6g.mat");
objectives_SOEC_model6g_15 = load("objectives_SOEC_model6g_15.mat");

obj6g_1_5 = [objectives_SOEC_model6g_15.objective_values(1:158,1);objectives_SOEC_model6g_15.objective_values(160:365,1)];

obj6g_2_5 = [objectives_SOEC_model6g.objective_values(1:158,1);objectives_SOEC_model6g.objective_values(160:365,1)];
obj6g_3_5 = [objectives_SOEC_model6g.objective_values(1:158,2);objectives_SOEC_model6g.objective_values(160:365,2)];
obj6g_4_5 = [objectives_SOEC_model6g.objective_values(1:158,3);objectives_SOEC_model6g.objective_values(160:365,3)];
obj6g_5_5 = [objectives_SOEC_model6g.objective_values(1:158,4);objectives_SOEC_model6g.objective_values(160:365,4)];


figure
boxplot(obj6g_1_5);  % Create the boxplot
title({"Revenues per day for an SOEC model6g in 2019"},{"with hydrogen price = 1.5€/kg"})
hold on
plot(mean(obj6g_1_5), 'xg')
hold off
print -depsc boxplot_SOEC5e_1.5.eps


figure
boxplot(obj6g_2_5);  % Create the boxplot
title({"Revenues per day for an SOEC model6g in 2019"},{"with hydrogen price = 2.5€/kg"})
hold on
plot(mean(obj6g_2_5), 'xg')
hold off
print -depsc boxplot_SOEC5e_2.5.eps

figure
boxplot(obj6g_3_5);  % Create the boxplot
title({"Revenues per day for an SOEC model6g in 2019"},{"with hydrogen price = 3.5€/kg"})
hold on
plot(mean(obj6g_3_5), 'xg')
hold off
print -depsc boxplot_SOEC6g_3.5.eps

figure
boxplot(obj6g_4_5);  % Create the boxplot
title({"Revenues per day for an SOEC model6g in 2019"},{"with hydrogen price = 4.5€/kg"})
hold on
plot(mean(obj6g_4_5), 'xg')
hold off
print -depsc boxplot_SOEC6g_4.5.eps

figure
boxplot(obj6g_5_5);  % Create the boxplot
title({"Revenues per day for an SOEC model5e in 2019"},{"with hydrogen price = 5.5€/kg"})
hold on
plot(mean(obj6g_5_5), 'xg')
hold off
print -depsc boxplot_SOEC6g_5.5.eps

mean_rev_6g_1_5 = mean(obj6g_1_5)
mean_rev_6g_2_5 = mean(obj6g_2_5)
mean_rev_6g_3_5 = mean(obj6g_3_5)
mean_rev_6g_4_5 = mean(obj6g_4_5)
mean_rev_6g_5_5 = mean(obj6g_5_5)






HHV_H2 = 141.88e6; %J/kg Higher heating value hydrogen
LHV_H2 = 119.96e6; %J/kg Lower heating value hydrogen

indices = find(abs(objectives_SOEC_model6g_15.objective_values(:,5)) > 10^-6);
mH2_6g_1_5 = objectives_SOEC_model6g_15.objective_values(indices,5);
p_E_6g_1_5 = objectives_SOEC_model6g_15.objective_values(indices,9);

indices = find(abs(objectives_SOEC_model6g.objective_values(:,5)) > 10^-6);
mH2_6g_2_5 = objectives_SOEC_model6g.objective_values(indices,5);
p_E_6g_2_5 = objectives_SOEC_model6g.objective_values(indices,9);

indices = find(abs(objectives_SOEC_model6g.objective_values(:,6)) > 10^-6);
mH2_6g_3_5 = objectives_SOEC_model6g.objective_values(indices,6);
p_E_6g_3_5 = objectives_SOEC_model6g.objective_values(indices,10);

indices = find(abs(objectives_SOEC_model6g.objective_values(:,7)) > 10^-6);
mH2_6g_4_5 = objectives_SOEC_model6g.objective_values(indices,7);
p_E_6g_4_5 = objectives_SOEC_model6g.objective_values(indices,11);


indices = find(abs(objectives_SOEC_model6g.objective_values(:,8)) > 10^-6);
mH2_6g_5_5 = objectives_SOEC_model6g.objective_values(indices,8);
p_E_6g_5_5 = objectives_SOEC_model6g.objective_values(indices,12);


eff_6g_1_5 = zeros(length(mH2_6g_1_5),1);
eff_6g_2_5 = zeros(length(mH2_6g_2_5),1);
eff_6g_3_5 = zeros(length(mH2_6g_3_5),1);
eff_6g_4_5 = zeros(length(mH2_6g_4_5),1);
eff_6g_5_5 = zeros(length(mH2_6g_5_5),1);
for i = 1:length(mH2_6g_1_5)
    eff_6g_1_5(i) =  mH2_6g_1_5(i)*LHV_H2/p_E_6g_1_5(i);
end
for i = 1:length(mH2_6g_2_5)
    eff_6g_2_5(i) =  mH2_6g_2_5(i)*LHV_H2/p_E_6g_2_5(i);
end
for i = 1:length(mH2_6g_3_5)
    eff_6g_3_5(i) =  mH2_6g_3_5(i)*LHV_H2/p_E_6g_3_5(i);
end
for i = 1:length(mH2_6g_4_5)
    eff_6g_4_5(i) =  mH2_6g_4_5(i)*LHV_H2/p_E_6g_4_5(i);
end
for i = 1:length(mH2_6g_5_5)
    eff_6g_5_5(i) =  mH2_6g_5_5(i)*LHV_H2/p_E_6g_5_5(i);
end

% for i = 1:8760
% %     eff_6g_2_5(i) =  mH2_6g_2_5(i)*LHV_H2/p_E_6g_2_5(i);
% %     eff_6g_3_5(i) =  mH2_6g_3_5(i)*LHV_H2/p_E_6g_3_5(i);
%     eff_6g_4_5(i) =  mH2_6g_4_5(i)*LHV_H2/p_E_6g_4_5(i);
%     eff_6g_5_5(i) =  mH2_6g_5_5(i)*LHV_H2/p_E_6g_5_5(i);
% end
eff_6g_1_5 = eff_6g_1_5(~isnan(eff_6g_1_5));
eff_6g_2_5 = eff_6g_2_5(~isnan(eff_6g_2_5));
eff_6g_3_5 = eff_6g_3_5(~isnan(eff_6g_3_5));
eff_6g_4_5 = eff_6g_4_5(~isnan(eff_6g_4_5));
eff_6g_5_5 = eff_6g_5_5(~isnan(eff_6g_5_5));
lower_limit = -1e-6;
upper_limit = 10;

mean_eff_6g_1_5 = mean(eff_6g_1_5(eff_6g_1_5 > lower_limit & eff_6g_1_5 < upper_limit));
mean_eff_6g_2_5 = mean(eff_6g_2_5(eff_6g_2_5 > lower_limit & eff_6g_2_5 < upper_limit));
mean_eff_6g_3_5 = mean(eff_6g_3_5(eff_6g_3_5 > lower_limit & eff_6g_3_5 < upper_limit));
mean_eff_6g_4_5 = mean(eff_6g_4_5(eff_6g_4_5 > lower_limit & eff_6g_4_5 < upper_limit));
mean_eff_6g_5_5 = mean(eff_6g_5_5(eff_6g_5_5 > lower_limit & eff_6g_5_5 < upper_limit));
revenue_week20_perday = sum(obj6g_2_5(141:147))
revenue_week20_perweek = 87251.88 % nog niet klaar met runnen incumbent al lang 72713.63, best bound 72741.3
rel_diff = (revenue_week20_perweek-revenue_week20_perday)/revenue_week20_perday

revenue_week20_perday = sum(obj6g_5_5(141:147))
revenue_week20_perweek = 317169.88 % wel heel rap klaar met runnen 
rel_diff = (revenue_week20_perweek-revenue_week20_perday)/revenue_week20_perday


%% Model 8g
objectives_SOEC_model8g = load("objectives_SOEC_model8g.mat");
objectives_SOEC_model8g_15 = load("objectives_SOEC_model8g_15.mat");

obj8g_1_5 = [objectives_SOEC_model8g_15.objective_values(1:158,1);objectives_SOEC_model8g_15.objective_values(160:365,1)];

obj8g_2_5 = [objectives_SOEC_model8g.objective_values(1:158,1);objectives_SOEC_model8g.objective_values(160:365,1)];
obj8g_3_5 = [objectives_SOEC_model8g.objective_values(1:158,2);objectives_SOEC_model8g.objective_values(160:365,2)];
obj8g_4_5 = [objectives_SOEC_model8g.objective_values(1:158,3);objectives_SOEC_model8g.objective_values(160:365,3)];
obj8g_5_5 = [objectives_SOEC_model8g.objective_values(1:158,4);objectives_SOEC_model8g.objective_values(160:365,4)];


figure
boxplot(obj8g_1_5);  % Create the boxplot
title({"Revenues per day for an SOEC model8g in 2019"},{"with hydrogen price = 1.5€/kg"})
hold on
plot(mean(obj8g_1_5), 'xg')
hold off
print -depsc boxplot_SOEC5e_1.5.eps


figure
boxplot(obj8g_2_5);  % Create the boxplot
title({"Revenues per day for an SOEC model8g in 2019"},{"with hydrogen price = 2.5€/kg"})
hold on
plot(mean(obj8g_2_5), 'xg')
hold off
print -depsc boxplot_SOEC5e_2.5.eps

figure
boxplot(obj8g_3_5);  % Create the boxplot
title({"Revenues per day for an SOEC model8g in 2019"},{"with hydrogen price = 3.5€/kg"})
hold on
plot(mean(obj8g_3_5), 'xg')
hold off
print -depsc boxplot_SOEC8g_3.5.eps

figure
boxplot(obj8g_4_5);  % Create the boxplot
title({"Revenues per day for an SOEC model8g in 2019"},{"with hydrogen price = 4.5€/kg"})
hold on
plot(mean(obj8g_4_5), 'xg')
hold off
print -depsc boxplot_SOEC8g_4.5.eps

figure
boxplot(obj8g_5_5);  % Create the boxplot
title({"Revenues per day for an SOEC model5e in 2019"},{"with hydrogen price = 5.5€/kg"})
hold on
plot(mean(obj8g_5_5), 'xg')
hold off
print -depsc boxplot_SOEC8g_5.5.eps

mean_rev_8g_1_5 = mean(obj8g_1_5)
mean_rev_8g_2_5 = mean(obj8g_2_5);
mean_rev_8g_3_5 = mean(obj8g_3_5);
mean_rev_8g_4_5 = mean(obj8g_4_5);
mean_rev_8g_5_5 = mean(obj8g_5_5);







HHV_H2 = 141.88e6; %J/kg Higher heating value hydrogen
LHV_H2 = 119.96e6; %J/kg Lower heating value hydrogen
mH2_8g_2_5 = objectives_SOEC_model8g.objective_values(:,5);
mH2_8g_3_5 = objectives_SOEC_model8g.objective_values(:,6);
mH2_8g_4_5 = objectives_SOEC_model8g.objective_values(:,7);
mH2_8g_5_5 = objectives_SOEC_model8g.objective_values(:,8);
p_E_8g_2_5 = objectives_SOEC_model8g.objective_values(:,9);
p_E_8g_3_5 = objectives_SOEC_model8g.objective_values(:,10);
p_E_8g_4_5 = objectives_SOEC_model8g.objective_values(:,11);
p_E_8g_5_5 = objectives_SOEC_model8g.objective_values(:,12);
eff_8g_2_5 = zeros(length(mH2_8g_5_5),1);
eff_8g_3_5 = zeros(length(mH2_8g_5_5),1);
eff_8g_4_5 = zeros(length(mH2_8g_5_5),1);
eff_8g_5_5 = zeros(length(mH2_8g_5_5),1);
for i = 1:8760
    eff_8g_2_5(i) =  mH2_8g_2_5(i)*LHV_H2/p_E_8g_2_5(i);
    eff_8g_3_5(i) =  mH2_8g_3_5(i)*LHV_H2/p_E_8g_3_5(i);
    eff_8g_4_5(i) =  mH2_8g_4_5(i)*LHV_H2/p_E_8g_4_5(i);
    eff_8g_5_5(i) =  mH2_8g_5_5(i)*LHV_H2/p_E_8g_5_5(i);
end
eff_8g_2_5 = eff_8g_2_5(~isnan(eff_8g_2_5));
eff_8g_3_5 = eff_8g_3_5(~isnan(eff_8g_3_5));
eff_8g_4_5 = eff_8g_4_5(~isnan(eff_8g_4_5));
eff_8g_5_5 = eff_8g_5_5(~isnan(eff_8g_5_5));
mean_eff_8g_2_5 =  mean(eff_8g_2_5(eff_8g_2_5 ~=0));
mean_eff_8g_3_5 =  mean(eff_8g_3_5(eff_8g_3_5 ~=0));
mean_eff_8g_4_5= mean(eff_8g_4_5(eff_8g_4_5 ~=0));
mean_eff_8g_5_5 =  mean(eff_8g_5_5(eff_8g_5_5 ~=0));


% %% 2.5
% data = [obj1g_2_5, obj2g_2_5, obj5g_2_5, obj6g_2_5];
% 
% % Create a figure
% figure;
% 
% % Plot the boxplots
% boxplot(data, 'Labels', {'a)', 'b)', 'c)', 'd)'});
% ylabel('Profit per day [€]'); % Add a y-axis label
% 
% % Set the title with the second line in bold
% title(sprintf('Profit per day for an SOE, a hydrogen price of 2.5 €/kg and\n\\bf{a) 1 power segment, b) 4 power segments, c) 4 power segments with heat integration} \n\\bf{from a LWR and d) 4 power segments with heat integration from a HTGR}'));
% 
% % Add mean value as a red diamond
% hold on;
% for i = 1:size(data, 2)
%     x = i; % x-coordinate for the mean value
%     y = mean(data(:, i)); % y-coordinate for the mean value
%     plot(x, y, 'rd', 'MarkerSize', 4, 'MarkerFaceColor', 'r'); % Plot red diamond at (x, y)
% end
% hold off;
% 
% % Set other optional properties as desired (e.g., axis limits, etc.)
% 
% % Save the figure (optional)
% saveas(gcf, 'boxplot_SOE_25.png');
% % Set other optional properties as desired (e.g., axis limits, etc.)
% yticklabels(num2str(get(gca, 'YTick').'));
% 
% %% 3.5
% data = [obj1g_3_5, obj2g_3_5, obj5g_3_5, obj6g_3_5];
% 
% % Create a figure
% figure;
% 
% % Plot the boxplots
% boxplot(data, 'Labels', {'a)', 'b)', 'c)', 'd)'});
% ylabel('Profit per day [€]'); % Add a y-axis label
% 
% % Set the title with the second line in bold
% title(sprintf('Profit per day for an SOE, a hydrogen price of 3.5 €/kg and\n\\bf{a) 1 power segment, b) 4 power segments, c) 4 power segments with heat integration} \n\\bf{from a LWR and d) 4 power segments with heat integration from a HTGR}'));
% 
% % Add mean value as a red diamond
% hold on;
% for i = 1:size(data, 2)
%     x = i; % x-coordinate for the mean value
%     y = mean(data(:, i)); % y-coordinate for the mean value
%     plot(x, y, 'rd', 'MarkerSize', 4, 'MarkerFaceColor', 'r'); % Plot red diamond at (x, y)
% end
% hold off;
% 
% % Set other optional properties as desired (e.g., axis limits, etc.)
% 
% % Save the figure (optional)
% saveas(gcf, 'boxplot_SOE_25.png');
% % Set other optional properties as desired (e.g., axis limits, etc.)
% yticklabels(num2str(get(gca, 'YTick').'));
% 
% %% 4.5
% data = [obj1g_4_5, obj2g_4_5, obj5g_4_5, obj6g_4_5];
% 
% % Create a figure
% figure;
% 
% % Plot the boxplots
% boxplot(data, 'Labels', {'a)', 'b)', 'c)', 'd)'});
% ylabel('Profit per day [€]'); % Add a y-axis label
% 
% % Set the title with the second line in bold
% title(sprintf('Profit per day for an SOE, a hydrogen price of 4.5 €/kg and\n\\bf{a) 1 power segment, b) 4 power segments, c) 4 power segments with heat integration} \n\\bf{from a LWR and d) 4 power segments with heat integration from a HTGR}'));
% 
% % Add mean value as a red diamond
% hold on;
% for i = 1:size(data, 2)
%     x = i; % x-coordinate for the mean value
%     y = mean(data(:, i)); % y-coordinate for the mean value
%     plot(x, y, 'rd', 'MarkerSize', 4, 'MarkerFaceColor', 'r'); % Plot red diamond at (x, y)
% end
% hold off;
% 
% % Set other optional properties as desired (e.g., axis limits, etc.)
% 
% % Save the figure (optional)
% saveas(gcf, 'boxplot_SOE_25.png');
% % Set other optional properties as desired (e.g., axis limits, etc.)
% yticklabels(num2str(get(gca, 'YTick').'));
% 
% %% 5.5
% data = [obj1g_5_5, obj2g_5_5, obj5g_5_5, obj6g_5_5];
% 
% % Create a figure
% figure;
% 
% % Plot the boxplots
% boxplot(data, 'Labels', {'a)', 'b)', 'c)', 'd)'});
% ylabel('Profit per day [€]'); % Add a y-axis label
% 
% % Set the title with the second line in bold
% title(sprintf('Profit per day for an SOE, a hydrogen price of 5.5 €/kg and\n\\bf{a) 1 power segment, b) 4 power segments, c) 4 power segments with heat integration} \n\\bf{from a LWR and d) 4 power segments with heat integration from a HTGR}'));
% 
% % Add mean value as a red diamond
% hold on;
% for i = 1:size(data, 2)
%     x = i; % x-coordinate for the mean value
%     y = mean(data(:, i)); % y-coordinate for the mean value
%     plot(x, y, 'rd', 'MarkerSize', 4, 'MarkerFaceColor', 'r'); % Plot red diamond at (x, y)
% end
% hold off;
% 
% % Set other optional properties as desired (e.g., axis limits, etc.)
% 
% % Save the figure (optional)
% saveas(gcf, 'boxplot_SOE_25.png');
% % Set other optional properties as desired (e.g., axis limits, etc.)
% yticklabels(num2str(get(gca, 'YTick').'));
% 
% 

%% 1.5
data = [obj1g_1_5, obj2g_1_5, obj5g_1_5, obj6g_1_5];

% Create a figure
figure;

% Plot the boxplots
boxplot(data, 'Labels', {'1)', '2)', '3)', '4)'});
ylabel('Profit per day [€]', FontSize=12); % Add a y-axis label

% Set the title with the second line in bold
title(sprintf('Profit per day for an SOE with\n a hydrogen price of 1.5 €/kg'), FontSize=14);

% Add mean value as a red diamond
hold on;
for i = 1:size(data, 2)
    x = i; % x-coordinate for the mean value
    y = mean(data(:, i)); % y-coordinate for the mean value
    plot(x, y, 'rd', 'MarkerSize', 4, 'MarkerFaceColor', 'r'); % Plot red diamond at (x, y)
end
hold off;
ylim([-3000,60000])
% Set other optional properties as desired (e.g., axis limits, etc.)

% Save the figure (optional)
saveas(gcf, 'boxplot_SOE_15.png');
% Set other optional properties as desired (e.g., axis limits, etc.)
yticklabels(num2str(get(gca, 'YTick').'));
print -depsc boxplot_combined_SOEC_15.eps
%% 2.5
data = [obj1g_2_5, obj2g_2_5, obj5g_2_5, obj6g_2_5];

% Create a figure
figure;

% Plot the boxplots
boxplot(data, 'Labels', {'1)', '2)', '3)', '4)'});
ylabel('Profit per day [€]', FontSize=12); % Add a y-axis label

% Set the title with the second line in bold
title(sprintf('Profit per day for an SOE with\n a hydrogen price of 2.5 €/kg'), FontSize=14);

% Add mean value as a red diamond
hold on;
for i = 1:size(data, 2)
    x = i; % x-coordinate for the mean value
    y = mean(data(:, i)); % y-coordinate for the mean value
    plot(x, y, 'rd', 'MarkerSize', 4, 'MarkerFaceColor', 'r'); % Plot red diamond at (x, y)
end
hold off;
ylim([-3000,30000])
% Set other optional properties as desired (e.g., axis limits, etc.)

% Save the figure (optional)
saveas(gcf, 'boxplot_SOE_25.png');
% Set other optional properties as desired (e.g., axis limits, etc.)
yticklabels(num2str(get(gca, 'YTick').'));
print -depsc boxplot_combined_SOEC_25.eps
%% 3.5
data = [obj1g_3_5, obj2g_3_5, obj5g_3_5, obj6g_3_5];

% Create a figure
figure;

% Plot the boxplots
boxplot(data, 'Labels', {'1)', '2)', '3)', '4)'});
ylabel('Profit per day [€]', FontSize=12); % Add a y-axis label

% Set the title with the second line in bold
title(sprintf('Profit per day for an SOE with\n a hydrogen price of 3.5 €/kg'), FontSize=14);

% Add mean value as a red diamond
hold on;
for i = 1:size(data, 2)
    x = i; % x-coordinate for the mean value
    y = mean(data(:, i)); % y-coordinate for the mean value
    plot(x, y, 'rd', 'MarkerSize', 4, 'MarkerFaceColor', 'r'); % Plot red diamond at (x, y)
end
hold off;
ylim([-3000,60000])
% Set other optional properties as desired (e.g., axis limits, etc.)

% Save the figure (optional)
saveas(gcf, 'boxplot_SOE_25.png');
% Set other optional properties as desired (e.g., axis limits, etc.)
yticklabels(num2str(get(gca, 'YTick').'));
print -depsc boxplot_combined_SOEC_35.eps

%% 4.5
data = [obj1g_4_5, obj2g_4_5, obj5g_4_5, obj6g_4_5];

% Create a figure
figure;

% Plot the boxplots
boxplot(data, 'Labels', {'1)', '2)', '3)', '4)'});
ylabel('Profit per day [€]', FontSize=12); % Add a y-axis label

% Set the title with the second line in bold
title(sprintf('Profit per day for an SOE with\n a hydrogen price of 4.5 €/kg'), FontSize=14);

% Add mean value as a red diamond
hold on;
for i = 1:size(data, 2)
    x = i; % x-coordinate for the mean value
    y = mean(data(:, i)); % y-coordinate for the mean value
    plot(x, y, 'rd', 'MarkerSize', 4, 'MarkerFaceColor', 'r'); % Plot red diamond at (x, y)
end
hold off;
ylim([-3000,60000])
% Set other optional properties as desired (e.g., axis limits, etc.)

% Save the figure (optional)
saveas(gcf, 'boxplot_SOE_25.png');
% Set other optional properties as desired (e.g., axis limits, etc.)
yticklabels(num2str(get(gca, 'YTick').'));
print -depsc boxplot_combined_SOEC_45.eps

%% 5.5
data = [obj1g_5_5, obj2g_5_5, obj5g_5_5, obj6g_5_5];

% Create a figure
figure;

% Plot the boxplots
boxplot(data, 'Labels', {'1)', '2)', '3)', '4)'});
ylabel('Profit per day [€]', FontSize=12); % Add a y-axis label

% Set the title with the second line in bold
title(sprintf('Profit per day for an SOE with\n a hydrogen price of 5.5 €/kg'), FontSize=14);

% Add mean value as a red diamond
hold on;
for i = 1:size(data, 2)
    x = i; % x-coordinate for the mean value
    y = mean(data(:, i)); % y-coordinate for the mean value
    plot(x, y, 'rd', 'MarkerSize', 4, 'MarkerFaceColor', 'r'); % Plot red diamond at (x, y)
end
hold off;
ylim([-3000,60000])
% Set other optional properties as desired (e.g., axis limits, etc.)

% Save the figure (optional)
saveas(gcf, 'boxplot_SOE_25.png');
% Set other optional properties as desired (e.g., axis limits, etc.)
yticklabels(num2str(get(gca, 'YTick').'));
print -depsc boxplot_combined_SOEC_55.eps


%% all gemiddeldes in 1 figuur, per model een lijntje in functie van waterstofprijs
hprice = [1.5 2.5 3.5 4.5 5.5];
figure
hold on
plot(hprice, [mean_rev_1g_1_5 mean_rev_1g_2_5 mean_rev_1g_3_5 mean_rev_1g_4_5 mean_rev_1g_5_5])
plot(hprice, [mean_rev_2g_1_5 mean_rev_2g_2_5 mean_rev_2g_3_5 mean_rev_2g_4_5 mean_rev_2g_5_5])
plot(hprice, [mean_rev_5g_1_5 mean_rev_5g_2_5 mean_rev_5g_3_5 mean_rev_5g_4_5 mean_rev_5g_5_5])
plot(hprice, [mean_rev_6g_1_5 mean_rev_6g_2_5 mean_rev_6g_3_5 mean_rev_6g_4_5 mean_rev_6g_5_5])
legend('1 power segment','4 power segments','4 power segments with heat from LWR','4 power segments with heat from HTGR', Location='southeast')
ylabel('Average profit per day [€]')
xlabel('Hydrogen price [€/kg]')
title(sprintf('Average profit per day for an SOE for different models\n as a function of hydrogen price'))
print -depsc averageprofit_SOEC.eps



%% figuur gemiddelde verbeteringen van meer segmenten ZONDER WARMTE
hprice = [2.5 3.5 4.5 5.5];
figure
hold on
plot(hprice, [(mean_rev_2g_2_5-mean_rev_1g_2_5)/mean_rev_1g_2_5 (mean_rev_2g_3_5-mean_rev_1g_3_5)/mean_rev_1g_3_5 (mean_rev_2g_4_5-mean_rev_1g_4_5)/mean_rev_1g_4_5 (mean_rev_2g_5_5-mean_rev_1g_5_5)/mean_rev_1g_5_5]*100)
ylabel('Relative improvement [%]', FontSize=12)
xlabel('Hydrogen price [€/kg]', FontSize=12)
title(sprintf('Relative improvement in representation of\n mean profit per day by increasing the number\n of power segments from 1 to 4'), FontSize=14)
print -depsc improvement_segments_SOEC.eps

%% figuur gemiddelde verbeteringen van meer segmenten MET WARMTE VAN HTGR
hprice = [1.5 2.5 3.5 4.5 5.5];
figure
hold on
plot(hprice, [(mean_rev_6g_1_5-mean_rev_8g_1_5)/mean_rev_8g_1_5 (mean_rev_6g_2_5-mean_rev_8g_2_5)/mean_rev_8g_2_5 (mean_rev_6g_3_5-mean_rev_8g_3_5)/mean_rev_8g_3_5 (mean_rev_6g_4_5-mean_rev_8g_4_5)/mean_rev_8g_4_5 (mean_rev_6g_5_5-mean_rev_8g_5_5)/mean_rev_8g_5_5]*100)
ylabel('Relative improvement in mean profit per day [%]')
xlabel('Hydrogen price [€/kg]')
title(sprintf('Relative improvement in mean profit per day\n by increasing the number of segments from 1 to 4'))


%% figuur gemiddelde verbeteringen van toevoegen warmte
hprice = [2.5 3.5 4.5 5.5];
figure
hold on
plot(hprice, [(mean_rev_5g_2_5-mean_rev_2g_2_5)/mean_rev_2g_2_5 (mean_rev_5g_3_5-mean_rev_2g_3_5)/mean_rev_2g_3_5 (mean_rev_5g_4_5-mean_rev_2g_4_5)/mean_rev_2g_4_5 (mean_rev_5g_5_5-mean_rev_2g_5_5)/mean_rev_2g_5_5]*100)
plot(hprice, [(mean_rev_6g_2_5-mean_rev_2g_2_5)/mean_rev_2g_2_5 (mean_rev_6g_3_5-mean_rev_2g_3_5)/mean_rev_2g_3_5 (mean_rev_6g_4_5-mean_rev_2g_4_5)/mean_rev_2g_4_5 (mean_rev_6g_5_5-mean_rev_2g_5_5)/mean_rev_2g_5_5]*100)
legend('4 power segments with heat from LWR','4 power segments with heat from HTGR', Location='northeast', FontSize=13)
ylabel('Relative improvement [%]', FontSize=12)
xlabel('Hydrogen price [€/kg]', FontSize=12)
title(sprintf('Relative improvement in mean profit per day\n due to heat integration'), FontSize=14)
print -depsc improvement_heat_SOEC.eps

%% Gemiddelde lengte in standby model 1g 1.5€/kg

% MOET BEGINNEN EN EINDIGEN MET EEN 1
% Example vector
vector = [objectives_SOEC_model1g_15.objective_values(1:158,5);objectives_SOEC_model1g_15.objective_values(160:365,5)]; % Your vector of length 8760

% Define the range
lower_limit = 0 - 10^-5;
upper_limit = 0 + 10^-5;

% Find the indices where the values are outside the range
nonzero_indices = find(vector < lower_limit | vector > upper_limit);
% Compute the differences between consecutive indices
diff_indices = diff(nonzero_indices);
diff_indices = diff_indices(diff_indices~=1)-1;
% Find the average length of a series of zeros
av_SB_1g_15 = mean(diff_indices);
%% Gemiddelde lengte in standby model 2g 1.5€/kg

% MOET BEGINNEN EN EINDIGEN MET EEN 1
% Example vector
vector = [objectives_SOEC_model2g_15.objective_values(1:158,5);objectives_SOEC_model2g_15.objective_values(160:365,5)]; % Your vector of length 8760

% Define the range
lower_limit = 0 - 10^-5;
upper_limit = 0 + 10^-5;

% Find the indices where the values are outside the range
nonzero_indices = find(vector < lower_limit | vector > upper_limit);
% Compute the differences between consecutive indices
diff_indices = diff(nonzero_indices);
diff_indices = diff_indices(diff_indices~=1)-1;
% Find the average length of a series of zeros
av_SB_2g_15 = mean(diff_indices);
%% Gemiddelde lengte in standby model 5g 1.5€/kg

% MOET BEGINNEN EN EINDIGEN MET EEN 1
% Example vector
vector = [1; objectives_SOEC_model5g_15.objective_values(1:158,5);objectives_SOEC_model5g_15.objective_values(160:365,5); 1]; % Your vector of length 8760

% Define the range
lower_limit = 0 - 10^-5;
upper_limit = 0 + 10^-5;

% Find the indices where the values are outside the range
nonzero_indices = find(vector < lower_limit | vector > upper_limit);
% Compute the differences between consecutive indices
diff_indices = diff(nonzero_indices);
diff_indices = diff_indices(diff_indices~=1)-1;
% Find the average length of a series of zeros
av_SB_5g_15 = mean(diff_indices);
%% Gemiddelde lengte in standby model 6g 1.5€/kg

% MOET BEGINNEN EN EINDIGEN MET EEN 1
% Example vector
vector = [objectives_SOEC_model6g_15.objective_values(1:158,5);objectives_SOEC_model6g_15.objective_values(160:365,5)]; % Your vector of length 8760

% Define the range
lower_limit = 0 - 10^-5;
upper_limit = 0 + 10^-5;

% Find the indices where the values are outside the range
nonzero_indices = find(vector < lower_limit | vector > upper_limit);
% Compute the differences between consecutive indices
diff_indices = diff(nonzero_indices);
diff_indices = diff_indices(diff_indices~=1)-1;
% Find the average length of a series of zeros
av_SB_6g_15 = mean(diff_indices);


%% Gemiddelde lengte in standby model 1g 2.5€/kg

% MOET BEGINNEN EN EINDIGEN MET EEN 1
% Example vector
vector = mH2_1g_2_5; % Your vector of length 8760

% Define the range
lower_limit = 0 - 10^-5;
upper_limit = 0 + 10^-5;

% Find the indices where the values are outside the range
nonzero_indices = find(vector < lower_limit | vector > upper_limit);
% Compute the differences between consecutive indices
diff_indices = diff(nonzero_indices);
diff_indices = diff_indices(diff_indices~=1)-1;
% Find the average length of a series of zeros
av_SB_1g_25 = mean(diff_indices);
%% Gemiddelde lengte in standby model 2g 2.5€/kg

% MOET BEGINNEN EN EINDIGEN MET EEN 1
% Example vector
vector = mH2_2g_2_5; % Your vector of length 8760

% Define the range
lower_limit = 0 - 10^-5;
upper_limit = 0 + 10^-5;

% Find the indices where the values are outside the range
nonzero_indices = find(vector < lower_limit | vector > upper_limit);
% Compute the differences between consecutive indices
diff_indices = diff(nonzero_indices);
diff_indices = diff_indices(diff_indices~=1)-1;
% Find the average length of a series of zeros
av_SB_2g_25 = mean(diff_indices);
%% Gemiddelde lengte in standby model 5g 2.5€/kg

% MOET BEGINNEN EN EINDIGEN MET EEN 1
% Example vector
vector = mH2_5g_2_5; % Your vector of length 8760

% Define the range
lower_limit = 0 - 10^-5;
upper_limit = 0 + 10^-5;

% Find the indices where the values are outside the range
nonzero_indices = find(vector < lower_limit | vector > upper_limit);
% Compute the differences between consecutive indices
diff_indices = diff(nonzero_indices);
diff_indices = diff_indices(diff_indices~=1)-1;
% Find the average length of a series of zeros
av_SB_5g_25 = mean(diff_indices);
%% Gemiddelde lengte in standby model 6g 2.5€/kg

% MOET BEGINNEN EN EINDIGEN MET EEN 1
% Example vector
vector = mH2_6g_2_5; % Your vector of length 8760

% Define the range
lower_limit = 0 - 10^-5;
upper_limit = 0 + 10^-5;

% Find the indices where the values are outside the range
nonzero_indices = find(vector < lower_limit | vector > upper_limit);
% Compute the differences between consecutive indices
diff_indices = diff(nonzero_indices);
diff_indices = diff_indices(diff_indices~=1)-1;
% Find the average length of a series of zeros
av_SB_6g_25 = mean(diff_indices);

%% Gemiddelde lengte in standby model 1g 3.5€/kg

% MOET BEGINNEN EN EINDIGEN MET EEN 1
% Example vector
vector = mH2_1g_3_5; % Your vector of length 8760

% Define the range
lower_limit = 0 - 10^-5;
upper_limit = 0 + 10^-5;

% Find the indices where the values are outside the range
nonzero_indices = find(vector < lower_limit | vector > upper_limit);
% Compute the differences between consecutive indices
diff_indices = diff(nonzero_indices);
diff_indices = diff_indices(diff_indices~=1)-1;
% Find the average length of a series of zeros
av_SB_1g_35 = mean(diff_indices);
%% Gemiddelde lengte in standby model 2g 3.5€/kg

% MOET BEGINNEN EN EINDIGEN MET EEN 1
% Example vector
vector = mH2_2g_3_5; % Your vector of length 8760

% Define the range
lower_limit = 0 - 10^-5;
upper_limit = 0 + 10^-5;

% Find the indices where the values are outside the range
nonzero_indices = find(vector < lower_limit | vector > upper_limit);
% Compute the differences between consecutive indices
diff_indices = diff(nonzero_indices);
diff_indices = diff_indices(diff_indices~=1)-1;
% Find the average length of a series of zeros
av_SB_2g_35 = mean(diff_indices);
%% Gemiddelde lengte in standby model 5g 3.5€/kg

% MOET BEGINNEN EN EINDIGEN MET EEN 1
% Example vector
vector = mH2_5g_3_5; % Your vector of length 8760

% Define the range
lower_limit = 0 - 10^-5;
upper_limit = 0 + 10^-5;

% Find the indices where the values are outside the range
nonzero_indices = find(vector < lower_limit | vector > upper_limit);
% Compute the differences between consecutive indices
diff_indices = diff(nonzero_indices);
diff_indices = diff_indices(diff_indices~=1)-1;
% Find the average length of a series of zeros
av_SB_5g_35 = mean(diff_indices);
%% Gemiddelde lengte in standby model 6g 3.5€/kg

% MOET BEGINNEN EN EINDIGEN MET EEN 1
% Example vector
vector = mH2_6g_3_5; % Your vector of length 8760

% Define the range
lower_limit = 0 - 10^-5;
upper_limit = 0 + 10^-5;

% Find the indices where the values are outside the range
nonzero_indices = find(vector < lower_limit | vector > upper_limit);
% Compute the differences between consecutive indices
diff_indices = diff(nonzero_indices);
diff_indices = diff_indices(diff_indices~=1)-1;
% Find the average length of a series of zeros
av_SB_6g_35 = mean(diff_indices);

%% Gemiddelde lengte in standby model 1g 4.5€/kg

% MOET BEGINNEN EN EINDIGEN MET EEN 1
% Example vector
vector = mH2_1g_4_5; % Your vector of length 8760

% Define the range
lower_limit = 0 - 10^-5;
upper_limit = 0 + 10^-5;

% Find the indices where the values are outside the range
nonzero_indices = find(vector < lower_limit | vector > upper_limit);
% Compute the differences between consecutive indices
diff_indices = diff(nonzero_indices);
diff_indices = diff_indices(diff_indices~=1)-1;
% Find the average length of a series of zeros
av_SB_1g_45 = mean(diff_indices);
%% Gemiddelde lengte in standby model 2g 4.5€/kg

% MOET BEGINNEN EN EINDIGEN MET EEN 1
% Example vector
vector = mH2_2g_4_5; % Your vector of length 8760

% Define the range
lower_limit = 0 - 10^-5;
upper_limit = 0 + 10^-5;

% Find the indices where the values are outside the range
nonzero_indices = find(vector < lower_limit | vector > upper_limit);
% Compute the differences between consecutive indices
diff_indices = diff(nonzero_indices);
diff_indices = diff_indices(diff_indices~=1)-1;
% Find the average length of a series of zeros
av_SB_2g_45 = mean(diff_indices);
%% Gemiddelde lengte in standby model 5g 4.5€/kg

% MOET BEGINNEN EN EINDIGEN MET EEN 1
% Example vector
vector = mH2_5g_4_5; % Your vector of length 8760

% Define the range
lower_limit = 0 - 10^-5;
upper_limit = 0 + 10^-5;

% Find the indices where the values are outside the range
nonzero_indices = find(vector < lower_limit | vector > upper_limit);
% Compute the differences between consecutive indices
diff_indices = diff(nonzero_indices);
diff_indices = diff_indices(diff_indices~=1)-1;
% Find the average length of a series of zeros
av_SB_5g_45 = mean(diff_indices);
%% Gemiddelde lengte in standby model 6g 4.5€/kg

% MOET BEGINNEN EN EINDIGEN MET EEN 1
% Example vector
vector = mH2_6g_4_5; % Your vector of length 8760

% Define the range
lower_limit = 0 - 10^-5;
upper_limit = 0 + 10^-5;

% Find the indices where the values are outside the range
nonzero_indices = find(vector < lower_limit | vector > upper_limit);
% Compute the differences between consecutive indices
diff_indices = diff(nonzero_indices);
diff_indices = diff_indices(diff_indices~=1)-1;
% Find the average length of a series of zeros
av_SB_6g_45 = mean(diff_indices);

%% Gemiddelde lengte in standby model 1g 5.5€/kg

% MOET BEGINNEN EN EINDIGEN MET EEN 1
% Example vector
vector = mH2_1g_5_5; % Your vector of length 8760

% Define the range
lower_limit = 0 - 10^-5;
upper_limit = 0 + 10^-5;

% Find the indices where the values are outside the range
nonzero_indices = find(vector < lower_limit | vector > upper_limit);
% Compute the differences between consecutive indices
diff_indices = diff(nonzero_indices);
diff_indices = diff_indices(diff_indices~=1)-1;
% Find the average length of a series of zeros
av_SB_1g_55 = mean(diff_indices);
%% Gemiddelde lengte in standby model 2g 5.5€/kg

% MOET BEGINNEN EN EINDIGEN MET EEN 1
% Example vector
vector = mH2_2g_5_5; % Your vector of length 8760

% Define the range
lower_limit = 0 - 10^-5;
upper_limit = 0 + 10^-5;

% Find the indices where the values are outside the range
nonzero_indices = find(vector < lower_limit | vector > upper_limit);
% Compute the differences between consecutive indices
diff_indices = diff(nonzero_indices);
diff_indices = diff_indices(diff_indices~=1)-1;
% Find the average length of a series of zeros
av_SB_2g_55 = mean(diff_indices);
%% Gemiddelde lengte in standby model 5g 5.5€/kg

% MOET BEGINNEN EN EINDIGEN MET EEN 1
% Example vector
vector = mH2_5g_5_5; % Your vector of length 8760

% Define the range
lower_limit = 0 - 10^-5;
upper_limit = 0 + 10^-5;

% Find the indices where the values are outside the range
nonzero_indices = find(vector < lower_limit | vector > upper_limit);
% Compute the differences between consecutive indices
diff_indices = diff(nonzero_indices);
diff_indices = diff_indices(diff_indices~=1)-1;
% Find the average length of a series of zeros
av_SB_5g_55 = mean(diff_indices);
%% Gemiddelde lengte in standby model 6g 5.5€/kg

% MOET BEGINNEN EN EINDIGEN MET EEN 1
% Example vector
vector = mH2_6g_5_5; % Your vector of length 8760

% Define the range
lower_limit = 0 - 10^-5;
upper_limit = 0 + 10^-5;

% Find the indices where the values are outside the range
nonzero_indices = find(vector < lower_limit | vector > upper_limit);
% Compute the differences between consecutive indices
diff_indices = diff(nonzero_indices);
diff_indices = diff_indices(diff_indices~=1)-1;
% Find the average length of a series of zeros
av_SB_6g_55 = mean(diff_indices);

%% figuur average efficiencies
hprice = [2.5 3.5 4.5 5.5];
figure
hold on
plot(hprice, [mean_eff_6g_2_5 mean_eff_6g_3_5 mean_eff_6g_4_5 mean_eff_6g_5_5]*100)
plot(hprice, [mean_eff_5g_2_5 mean_eff_5g_3_5 mean_eff_5g_4_5 mean_eff_5g_5_5]*100)
plot(hprice, [mean_eff_2g_2_5 mean_eff_2g_3_5 mean_eff_2g_4_5 mean_eff_2g_5_5]*100)
plot(hprice, [mean_eff_1g_2_5 mean_eff_1g_3_5 mean_eff_1g_4_5 mean_eff_1g_5_5]*100)



legend(fliplr({'1 power segment','4 power segments','4 power segments with heat from an LWR','4 power segments with heat from an HTGR'}), Location='east', FontSize=12)
ylabel('Average system efficiency [%]', FontSize=12)
xlabel('Hydrogen price [€/kg]', FontSize=12)
title(sprintf('Average system efficiency for an SOE for different\n models as a function of hydrogen price'), FontSize=14)
print -depsc averageefficiency_SOE.eps

