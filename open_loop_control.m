%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set up parameters needed in simulation            %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cap_per_m2 = 2500;
h_o = 0.3; % heat transfer coeff per m^2 for exterior walls
h_i = 0.9; % heat transfer coeff per m^2 for interior walls

areaA = 15; areaB = 8; areaC = 2;

capacA = cap_per_m2*areaA;
capacB = cap_per_m2*areaB;
capacC = cap_per_m2*areaC;

htAO = h_o * 11 * 3;  htBO = h_o * 9 * 3;  htCO = h_o * 0 * 3;
htAB = h_i * 3 * 3;   htAC = h_i * 2 * 3;  htBC = h_i * 4 * 3;

A = 100*[-1/capacA * (htAO + htAB + htAC), htAB / capacA,                    htAC / capacA;
         htAB / capacB,                    -1/capacB * (htBO + htAB + htBC), htBC / capacB;
         htAC / capacC,                    htBC / capacC,                    -1/capacC * (htCO + htAC + htBC)];
B = [1 / capacA, 0, 0; 0, 1/capacB, 0; 0, 0, 1/capacC];
C = [1, 0, 0; 0, 1, 0; 0, 0, 1];
D = [0, 0, 0; 0, 0, 0; 0, 0, 0];

T_conv = 300;
T_maint = 400;
t = 0:0.01:T_maint;
x_des = [20; 20; 20];
x0 = [0; 0; 0]; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run simulation with full-rank B matrix            %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
heat_sys = ss(A, B, C, D);
Wr_conv = gram(heat_sys, 'c', gramOptions('TimeIntervals', [0, T_conv]));
Wr_maint = gram(heat_sys, 'c', gramOptions('TimeIntervals', [0, T_maint]));
eta_conv = Wr_conv \ (x_des - expm(A * T_conv) * x0);
eta_maint = Wr_maint \ (x_des - expm(A * T_maint) * x_des);

u = zeros(3, sum(t <= T_maint));
for i=1:sum(t <= T_maint)
    u(:,i) = B'*expm(A'*(T_maint - t(i))) * eta_maint;
end
for i=1:sum(t <= T_conv)
    u(:,i) = B'*expm(A'*(T_conv - t(i))) * eta_conv;
end

figure(1);
animate_input_ux0(u, x0, t, heat_sys, x_des)

y = lsim(heat_sys, u', t, x0);

figure(2)
hold off; plot(t, y(:,1)); hold on; plot(t, y(:,2)); plot(t, y(:,3))
legend('Room 1', 'Room 2', 'Room 3')
title("Temperature in Rooms over Time with Full Heaters"); xlabel("Time (sec)"); ylabel("Temp in Room (°C)")
saveas(gcf, "figs/ol_fullB_temp.png")

figure(3)
hold off; plot(t, u(1,:)); hold on; plot(t, u(2,:)); plot(t, u(3,:))
legend('Room 1', 'Room 2', 'Room 3')
title("Heating Inputs over Time with Full Heaters"); xlabel("Time (sec)"); ylabel("Heating Input");
saveas(gcf, "figs/ol_fullB_pow.png")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run simulation with no heater in room 3           %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
B = [1 / capacA, 0, 0; 0, 1/capacB, 0; 0, 0, 0];
heat_sys = ss(A, B, C, D);
Wr_conv = gram(heat_sys, 'c', gramOptions('TimeIntervals', [0, T_conv]));
Wr_maint = gram(heat_sys, 'c', gramOptions('TimeIntervals', [0, T_maint]));
eta_conv = Wr_conv \ (x_des - expm(A * T_conv) * x0);
eta_maint = Wr_maint \ (x_des - expm(A * T_maint) * x_des);

u = zeros(3, sum(t <= T_maint));
for i=1:sum(t <= T_maint)
    u(:,i) = B'*expm(A'*(T_maint - t(i))) * eta_maint;
end
for i=1:sum(t <= T_conv)
    u(:,i) = B'*expm(A'*(T_conv - t(i))) * eta_conv;
end

figure(4)
animate_input_ux0(u, x0, t, heat_sys, x_des)
y = lsim(heat_sys, u', t, x0);

figure(5)
hold off; plot(t, y(:,1)); hold on; plot(t, y(:,2)); plot(t, y(:,3))
legend('Room 1', 'Room 2', 'Room 3')
title("Temperature in Rooms over Time"); xlabel("Time (sec)"); ylabel("Temp in Room (°C)")
saveas(gcf, "figs/ol_no3_temp.png")

figure(6)
hold off; plot(t, u(1,:)); hold on; plot(t, u(2,:)); plot(t, u(3,:))
legend('Room 1', 'Room 2', 'Room 3')
title("Heating Inputs over Time"); xlabel("Time (sec)"); ylabel("Heating Input")
saveas(gcf, "figs/ol_no3_pow.png")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run simulation with filtered input                %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T = 1000;
T_conv = 300;
T_maint = 100;
t = 0:0.01:T;

B = [1 / capacA, 0, 0; 0, 1/capacB, 0; 0, 0, 1/capacC];
heat_sys = ss(A, B, C, D);

Wr_conv = gram(heat_sys, 'c', gramOptions('TimeIntervals', [0, T_conv]));
Wr_maint = gram(heat_sys, 'c', gramOptions('TimeIntervals', [0, T_maint]));
eta_conv = Wr_conv \ (x_des - expm(A * T_conv) * x0);
eta_maint = Wr_maint \ (x_des - expm(A * T_maint) * x_des);

n_conv = sum(t < T_conv);
n_maint = sum(t < T_maint); 
u = zeros(3, n_conv);
u_conv = zeros(3, n_conv);

for i=1:n_conv
    u_conv(:,i) = B'*expm(A'*(T_conv - t(i))) * eta_conv;
end

n_maint = sum(t < T_maint);
u_maint = zeros(3, n_maint);
for i=1:n_maint
    u_maint(:,i) = B'*expm(A'*(T_maint - t(i))) * eta_maint;
end

u = zeros(3, n_pts);
for i=1:size(t, 2)
    if i <= n_conv
        u(:,i) = u_conv(:,i);
    else
        u(:,i) = u_maint(:, mod(i - n_conv, n_maint) + 1);
    end
end

filt = fir1(10000, 0.001);
u_lp = (filter(filt, 1, u'))';
u = u_lp;

figure(7);
disp(size(u)); disp(size(x0)); disp(size(t))
animate_input_ux0(u, x0, t, heat_sys, x_des)

y = lsim(heat_sys, u, t, x0);

figure(8)
hold off; plot(t, y(:,1)); hold on; plot(t, y(:,2)); plot(t, y(:,3))
lgd = legend('Room 1', 'Room 2', 'Room 3');
lgd.Location = 'northwest';
title("Temperature in Rooms over Time for Filtered u"); xlabel("Time (sec)"); ylabel("Temp in Room (°C)")
saveas(gcf, "figs/ol_filt_temp.png")

figure(9)
hold off; plot(t, u(1,:)); hold on; plot(t, u(2,:)); plot(t, u(3,:))
lgd = legend('Room 1', 'Room 2', 'Room 3');
lgd.Location = 'northwest';
title("Heating Inputs over Time for Filtered u"); xlabel("Time (sec)"); ylabel("Heating Input");
saveas(gcf, "figs/ol_filt_pow.png")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Helper functions                                  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = animate_input_ux0(u, x0, t, heat_sys, x_des) 
    n_pts = size(u, 2);
    y = lsim(heat_sys, u', t, x0);
    yt = y';

    max_t = max(max(y(:, 1:3)));
    colorbar('Ticks',[0, 1], 'TickLabels',{num2str(min(y, [], "all")), num2str(max(y, [], "all"))});
    energy_consumed = 0;
    for i = [1:max(floor(n_pts/100), 1):n_pts, n_pts]
        if i ~= 1
            energy_consumed = energy_consumed + sum(u(:,i).^2 * (t(i) - t(iprev)));
        end
        plot_rooms([temp2color(y(i, 1), max_t); temp2color(y(i, 2), max_t); temp2color(y(i, 3), max_t)])
        xlabel(strcat("Energy Consumed: ", formattedDisplayText(energy_consumed, 'NumericFormat', 'shortEng')));
        pause(0.01)
        title(strcat(num2str(t(i)), " seconds elapsed"));
        iprev = i;
    end

    tdiffs = t(2:end) - t(1:(end-1));
    disp(size(u)); 
    disp(size(t));

    disp(strcat("Total Energy:", num2str(sum(u(:,2:end).^2 * tdiffs'), '%0.4g')))
    disp(strcat("Max Overshoot:", num2str(max(y'-x_des*ones(1,n_pts), [],'all'), '%0.4g')))
    disp("End temperatures:"); disp(y(end, :))
end

function [] = plot_rooms(colors,outside_color)
    if nargin == 1
        outside_color = [0,0,0];
    end
    rectangle('Position', [-2,-1,14,7],'FaceColor', outside_color)
    rectangle('Position', [0,0,5,5], 'FaceColor', colors(1,:))
    text(2.5 , 2.5, '1')
    rectangle('Position', [5,0,4,5], 'FaceColor', colors(2,:))
    text(7, 3.75, '2')
    rectangle('Position', [5,0.5,2,2], 'FaceColor', colors(3,:))
    text(6, 1.5, '3')
%     xlim(-1,10)
%     ylim(-1,6)
    set(gca,'xtick',[])
    set(gca,'ytick',[])
end

function [c] = temp2color(t, t_max)
    idx = ceil(255 * max(0,t) / t_max) + 1;
    colors = parula(256);
    c = colors(min(idx, 256), :);
end