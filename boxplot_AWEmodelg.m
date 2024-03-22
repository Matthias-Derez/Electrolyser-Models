clearvars

%% Model 1g
objectives_AWE_model1g = load("objectives_AWE_model1g.mat");
objectives_PEM_model1g_15 = load("objectives_PEM_model1g_15.mat");

obj1g_1_5 = [objectives_PEM_model1g_15.objective_values(1:158,1);objectives_PEM_model1g_15.objective_values(160:365,1)];

obj1g_2_5 = objectives_AWE_model1g.objective_values(1:365,1);
obj1g_3_5 = objectives_AWE_model1g.objective_values(1:365,2);
obj1g_4_5 = objectives_AWE_model1g.objective_values(1:365,3);
obj1g_5_5 = objectives_AWE_model1g.objective_values(1:365,4);

figure
boxplot(obj1g_2_5);  % Create the boxplot
title({"Revenues per day for a AWE model1g in 2019"},{"with hydrogen price = 2.5€/kg"})
hold on
plot(mean(obj1g_2_5), 'xg')
hold off
print -depsc boxplot_AWE5e_2.5.eps

figure
boxplot(obj1g_3_5);  % Create the boxplot
title({"Revenues per day for a AWE model1g in 2019"},{"with hydrogen price = 3.5€/kg"})
hold on
plot(mean(obj1g_3_5), 'xg')
hold off
print -depsc boxplot_AWE1g_3.5.eps

figure
boxplot(obj1g_4_5);  % Create the boxplot
title({"Revenues per day for a AWE model1g in 2019"},{"with hydrogen price = 4.5€/kg"})
hold on
plot(mean(obj1g_4_5), 'xg')
hold off
print -depsc boxplot_AWE1g_4.5.eps

figure
boxplot(obj1g_5_5);  % Create the boxplot
title({"Revenues per day for a AWE model5e in 2019"},{"with hydrogen price = 5.5€/kg"})
hold on
plot(mean(obj1g_5_5), 'xg')
hold off
print -depsc boxplot_AWE1g_5.5.eps


mean_rev_1g_2_5 = mean(obj1g_2_5);
mean_rev_1g_3_5 = mean(obj1g_3_5);
mean_rev_1g_4_5 = mean(obj1g_4_5);
mean_rev_1g_5_5 = mean(obj1g_5_5);







HHV_H2 = 141.88e6; %J/kg Higher heating value hydrogen
LHV_H2 = 119.96e6; %J/kg Lower heating value hydrogen
mH2_1g_2_5 = objectives_AWE_model1g.objective_values(:,5);
mH2_1g_3_5 = objectives_AWE_model1g.objective_values(:,6);
mH2_1g_4_5 = objectives_AWE_model1g.objective_values(:,7);
mH2_1g_5_5 = objectives_AWE_model1g.objective_values(:,8);
p_E_1g_2_5 = objectives_AWE_model1g.objective_values(:,9);
p_E_1g_3_5 = objectives_AWE_model1g.objective_values(:,10);
p_E_1g_4_5 = objectives_AWE_model1g.objective_values(:,11);
p_E_1g_5_5 = objectives_AWE_model1g.objective_values(:,12);
eff_1g_2_5 = zeros(length(mH2_1g_5_5),1);
eff_1g_3_5 = zeros(length(mH2_1g_5_5),1);
eff_1g_4_5 = zeros(length(mH2_1g_5_5),1);
eff_1g_5_5 = zeros(length(mH2_1g_5_5),1);
for i = 1:8760
    eff_1g_2_5(i) =  mH2_1g_2_5(i)*HHV_H2/p_E_1g_2_5(i);
    eff_1g_3_5(i) =  mH2_1g_3_5(i)*HHV_H2/p_E_1g_3_5(i);
    eff_1g_4_5(i) =  mH2_1g_4_5(i)*HHV_H2/p_E_1g_4_5(i);
    eff_1g_5_5(i) =  mH2_1g_5_5(i)*HHV_H2/p_E_1g_5_5(i);
end
eff_1g_2_5 = eff_1g_2_5(~isnan(eff_1g_2_5));
eff_1g_3_5 = eff_1g_3_5(~isnan(eff_1g_3_5));
eff_1g_4_5 = eff_1g_4_5(~isnan(eff_1g_4_5));
eff_1g_5_5 = eff_1g_5_5(~isnan(eff_1g_5_5));
mean_eff_1g_2_5 =  mean(eff_1g_2_5(eff_1g_2_5 ~=0));
mean_eff_1g_3_5 =  mean(eff_1g_3_5(eff_1g_3_5 ~=0));
mean_eff_1g_4_5= mean(eff_1g_4_5(eff_1g_4_5 ~=0));
mean_eff_1g_5_5 =  mean(eff_1g_5_5(eff_1g_5_5 ~=0));
%% Model 2g
objectives_AWE_model2g = load("objectives_AWE_model2g.mat");
objectives_PEM_model2g_15 = load("objectives_PEM_model2g_15.mat");

obj2g_1_5 = [objectives_PEM_model2g_15.objective_values(1:158,1);objectives_PEM_model2g_15.objective_values(160:365,1)];

obj2g_2_5 = objectives_AWE_model2g.objective_values(1:365,1);
obj2g_3_5 = objectives_AWE_model2g.objective_values(1:365,2);
obj2g_4_5 = objectives_AWE_model2g.objective_values(1:365,3);
obj2g_5_5 = objectives_AWE_model2g.objective_values(1:365,4);

figure
boxplot(obj2g_2_5);  % Create the boxplot
title({"Revenues per day for a AWE model2g in 2019"},{"with hydrogen price = 2.5€/kg"})
hold on
plot(mean(obj2g_2_5), 'xg')
hold off
print -depsc boxplot_AWE5e_2.5.eps

figure
boxplot(obj2g_3_5);  % Create the boxplot
title({"Revenues per day for a AWE model2g in 2019"},{"with hydrogen price = 3.5€/kg"})
hold on
plot(mean(obj2g_3_5), 'xg')
hold off
print -depsc boxplot_AWE2g_3.5.eps

figure
boxplot(obj2g_4_5);  % Create the boxplot
title({"Revenues per day for a AWE model2g in 2019"},{"with hydrogen price = 4.5€/kg"})
hold on
plot(mean(obj2g_4_5), 'xg')
hold off
print -depsc boxplot_AWE2g_4.5.eps

figure
boxplot(obj2g_5_5);  % Create the boxplot
title({"Revenues per day for a AWE model5e in 2019"},{"with hydrogen price = 5.5€/kg"})
hold on
plot(mean(obj2g_5_5), 'xg')
hold off
print -depsc boxplot_AWE2g_5.5.eps


mean_rev_2g_2_5 = mean(obj2g_2_5);
mean_rev_2g_3_5 = mean(obj2g_3_5);
mean_rev_2g_4_5 = mean(obj2g_4_5);
mean_rev_2g_5_5 = mean(obj2g_5_5);







HHV_H2 = 141.88e6; %J/kg Higher heating value hydrogen
LHV_H2 = 119.96e6; %J/kg Lower heating value hydrogen
mH2_2g_2_5 = objectives_AWE_model2g.objective_values(:,5);
mH2_2g_3_5 = objectives_AWE_model2g.objective_values(:,6);
mH2_2g_4_5 = objectives_AWE_model2g.objective_values(:,7);
mH2_2g_5_5 = objectives_AWE_model2g.objective_values(:,8);
p_E_2g_2_5 = objectives_AWE_model2g.objective_values(:,9);
p_E_2g_3_5 = objectives_AWE_model2g.objective_values(:,10);
p_E_2g_4_5 = objectives_AWE_model2g.objective_values(:,11);
p_E_2g_5_5 = objectives_AWE_model2g.objective_values(:,12);
eff_2g_2_5 = zeros(length(mH2_2g_5_5),1);
eff_2g_3_5 = zeros(length(mH2_2g_5_5),1);
eff_2g_4_5 = zeros(length(mH2_2g_5_5),1);
eff_2g_5_5 = zeros(length(mH2_2g_5_5),1);
for i = 1:8760
    eff_2g_2_5(i) =  mH2_2g_2_5(i)*HHV_H2/p_E_2g_2_5(i);
    eff_2g_3_5(i) =  mH2_2g_3_5(i)*HHV_H2/p_E_2g_3_5(i);
    eff_2g_4_5(i) =  mH2_2g_4_5(i)*HHV_H2/p_E_2g_4_5(i);
    eff_2g_5_5(i) =  mH2_2g_5_5(i)*HHV_H2/p_E_2g_5_5(i);
end
eff_2g_2_5 = eff_2g_2_5(~isnan(eff_2g_2_5));
eff_2g_3_5 = eff_2g_3_5(~isnan(eff_2g_3_5));
eff_2g_4_5 = eff_2g_4_5(~isnan(eff_2g_4_5));
eff_2g_5_5 = eff_2g_5_5(~isnan(eff_2g_5_5));
mean_eff_2g_2_5 =  mean(eff_2g_2_5(eff_2g_2_5 ~=0));
mean_eff_2g_3_5 =  mean(eff_2g_3_5(eff_2g_3_5 ~=0));
mean_eff_2g_4_5= mean(eff_2g_4_5(eff_2g_4_5 ~=0));
mean_eff_2g_5_5 =  mean(eff_2g_5_5(eff_2g_5_5 ~=0));
revenue_week20_perday = sum(obj2g_2_5(141:147))
revenue_week20_perweek = 72713.63 % nog niet klaar met runnen incumbent al lang 72713.63, best bound 72741.3
rel_diff = (revenue_week20_perweek-revenue_week20_perday)/revenue_week20_perday

revenue_week20_perday = sum(obj2g_5_5(141:147))
revenue_week20_perweek = 302634.77029693074 % wel heel rap klaar met runnen 
rel_diff = (revenue_week20_perweek-revenue_week20_perday)/revenue_week20_perday


%% Model 5g
objectives_AWE_model5g = load("objectives_AWE_model5g.mat");
objectives_PEM_model5g_15 = load("objectives_PEM_model5g_15.mat");

obj5g_1_5 = [objectives_PEM_model5g_15.objective_values(1:158,1);objectives_PEM_model5g_15.objective_values(160:365,1)];

obj5g_2_5 = objectives_AWE_model5g.objective_values(1:365,1);
obj5g_3_5 = objectives_AWE_model5g.objective_values(1:365,2);
obj5g_4_5 = objectives_AWE_model5g.objective_values(1:365,3);
obj5g_5_5 = objectives_AWE_model5g.objective_values(1:365,4);

figure
boxplot(obj5g_2_5);  % Create the boxplot
title({"Revenues per day for a AWE model5g in 2019"},{"with hydrogen price = 2.5€/kg"})
hold on
plot(mean(obj5g_2_5), 'xg')
hold off
print -depsc boxplot_AWE5e_2.5.eps

figure
boxplot(obj5g_3_5);  % Create the boxplot
title({"Revenues per day for a AWE model5g in 2019"},{"with hydrogen price = 3.5€/kg"})
hold on
plot(mean(obj5g_3_5), 'xg')
hold off
print -depsc boxplot_AWE5g_3.5.eps

figure
boxplot(obj5g_4_5);  % Create the boxplot
title({"Revenues per day for a AWE model5g in 2019"},{"with hydrogen price = 4.5€/kg"})
hold on
plot(mean(obj5g_4_5), 'xg')
hold off
print -depsc boxplot_AWE5g_4.5.eps

figure
boxplot(obj5g_5_5);  % Create the boxplot
title({"Revenues per day for a AWE model5e in 2019"},{"with hydrogen price = 5.5€/kg"})
hold on
plot(mean(obj5g_5_5), 'xg')
hold off
print -depsc boxplot_AWE5g_5.5.eps


mean_rev_5g_2_5 = mean(obj5g_2_5);
mean_rev_5g_3_5 = mean(obj5g_3_5);
mean_rev_5g_4_5 = mean(obj5g_4_5);
mean_rev_5g_5_5 = mean(obj5g_5_5);







HHV_H2 = 141.88e6; %J/kg Higher heating value hydrogen
LHV_H2 = 119.96e6; %J/kg Lower heating value hydrogen
mH2_5g_2_5 = objectives_AWE_model5g.objective_values(:,5);
mH2_5g_3_5 = objectives_AWE_model5g.objective_values(:,6);
mH2_5g_4_5 = objectives_AWE_model5g.objective_values(:,7);
mH2_5g_5_5 = objectives_AWE_model5g.objective_values(:,8);
p_E_5g_2_5 = objectives_AWE_model5g.objective_values(:,9);
p_E_5g_3_5 = objectives_AWE_model5g.objective_values(:,10);
p_E_5g_4_5 = objectives_AWE_model5g.objective_values(:,11);
p_E_5g_5_5 = objectives_AWE_model5g.objective_values(:,12);
eff_5g_2_5 = zeros(length(mH2_5g_5_5),1);
eff_5g_3_5 = zeros(length(mH2_5g_5_5),1);
eff_5g_4_5 = zeros(length(mH2_5g_5_5),1);
eff_5g_5_5 = zeros(length(mH2_5g_5_5),1);
for i = 1:8760
    eff_5g_2_5(i) =  mH2_5g_2_5(i)*HHV_H2/p_E_5g_2_5(i);
    eff_5g_3_5(i) =  mH2_5g_3_5(i)*HHV_H2/p_E_5g_3_5(i);
    eff_5g_4_5(i) =  mH2_5g_4_5(i)*HHV_H2/p_E_5g_4_5(i);
    eff_5g_5_5(i) =  mH2_5g_5_5(i)*HHV_H2/p_E_5g_5_5(i);
end
eff_5g_2_5 = eff_5g_2_5(~isnan(eff_5g_2_5));
eff_5g_3_5 = eff_5g_3_5(~isnan(eff_5g_3_5));
eff_5g_4_5 = eff_5g_4_5(~isnan(eff_5g_4_5));
eff_5g_5_5 = eff_5g_5_5(~isnan(eff_5g_5_5));
mean_eff_5g_2_5 =  mean(eff_5g_2_5(eff_5g_2_5 ~=0));
mean_eff_5g_3_5 =  mean(eff_5g_3_5(eff_5g_3_5 ~=0));
mean_eff_5g_4_5= mean(eff_5g_4_5(eff_5g_4_5 ~=0));
mean_eff_5g_5_5 =  mean(eff_5g_5_5(eff_5g_5_5 ~=0));
revenue_week20_perday = sum(obj5g_2_5(141:147))
revenue_week20_perweek = 83506.65 % nog niet klaar met runnen incumbent al lang 72713.63, best bound 72741.3
rel_diff = (revenue_week20_perweek-revenue_week20_perday)/revenue_week20_perday

revenue_week20_perday = sum(obj5g_5_5(141:147))
revenue_week20_perweek = 313427.78 % wel heel rap klaar met runnen 
rel_diff = (revenue_week20_perweek-revenue_week20_perday)/revenue_week20_perday




%% Figuur efficienties
hprice = [2.5 3.5 4.5 5.5];
figure
hold on
plot(hprice, [mean_eff_1g_2_5 mean_eff_1g_3_5 mean_eff_1g_4_5 mean_eff_1g_5_5]*100)
plot(hprice, [mean_eff_5g_2_5 mean_eff_5g_3_5 mean_eff_5g_4_5 mean_eff_5g_5_5]*100)
plot(hprice, [mean_eff_6g_2_5 mean_eff_6g_3_5 mean_eff_6g_4_5 mean_eff_6g_5_5]*100)
legend('1 power segment','4 power segments','4 power segments with heat', Location='northeast')
ylabel('Average system efficiency [%]')
xlabel('Hydrogen price [€/kg]')
title(sprintf('Average system efficiency for an AWE for different models\n as a function of hydrogen price'))
print -depsc averageefficiency_AWE.eps
