% Heat Transfer - Spring 2023 Computational Assignment 1
% Michael Allen, Cullen Hirstius, Hayden Payne
% 3/5/23

clc; clear;  close all;
format short 

FIN_COUNT = [1, (3^2), (6^2), (9^2), (12^2)];
max_temp = zeros(size(FIN_COUNT));

for N = 1:length(FIN_COUNT)
    max_temp(N) = (T_base_solve(FIN_COUNT(N)));  
end

figure('Name', 'Size of fin array to Base Temp')  %initalize first figure
hold on

for i = 1:length(max_temp)
    bar(FIN_COUNT(i)^.5, max_temp(i))
end

xlabel('Size of Fin Array');
ylabel('Max Temperature at Base (°C)');
title('Size of fin array to Base Temp');

legendCell = strcat('N =',string(num2cell(FIN_COUNT .^ .5)))
legend(legendCell, "Location","eastoutside");

hold off


function T_base = T_base_solve(N)
    %System Properties
    P_chip = 30; %W
    h_avg = 125; %W/m^2K
    T_ambient = 20; %°C
    length_chip = 25 * 10^(-3); %m
    Area_chip = length_chip^2; %m^2
    
    
    %Fin Properties
    k_fin = 200; %W/m^2K
    L_fin = 15 * 10^(-3); %m
    diameter_fin = 2 * 10^(-3); %m
    P_fin = diameter_fin * pi; %m
    A_c = pi * (diameter_fin/2)^2; %m^2
    
    
    %Calcs for individual fin
    h_fin = h_avg;
    m = (h_fin*P_fin/(k_fin*A_c))^(1/2);
    mL = m*L_fin;
    M = @(theta_b) (h_fin * P_fin * k_fin * A_c)^(1/2) * theta_b; %f1, theta_b
    q_fin = @(theta_b) M(theta_b) * (sinh(mL) + (h_fin/(m*k_fin))*cosh(mL))/(cosh(mL) + (h_fin/m*k_fin)*sinh(mL)); %f2 <<< f1, theta_b
    
    
    %Calcs for N fins -- Fins are hottest at x = 0
    A_b = Area_chip - N*(A_c);
    Base_conv = @(theta_b) theta_b*h_avg*A_b; %f3 
    
    RHS = @(theta_b) q_fin(theta_b)*N + Base_conv(theta_b) - P_chip; %f4 <<< f3 & f2, theta b
    
    theta_b = fsolve(RHS, 0);
    T_base = theta_b + T_ambient; %Return base temp in °C where temp is hottest
    
end




