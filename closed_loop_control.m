%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set up coefficients needed in simulation          %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cap_per_m2 = 2500;
h_o = 0.3; % heat transfer coeff per m^2 for exterior walls
h_i = 0.9; % heat transfer coeff per m^2 for interior walls

areaA = 15; areaB = 8; areaC = 2;

capacA = cap_per_m2*areaA;
capacB = cap_per_m2*areaB;
capacC = cap_per_m2*areaC;


htAO = h_o * 11 * 3;  htBO = h_o * 9 * 3; htCO = h_o * 0 * 3;
htAB = h_i * 3 * 3;  htAC = h_i * 2 * 3; htBC = h_i * 4 * 3;

A = 100*[-1/capacA * (htAO + htAB + htAC), htAB / capacA,                    htAC / capacA;
     htAB / capacB,                    -1/capacB * (htBO + htAB + htBC), htBC / capacB;
     htAC / capacC,                    htBC / capacC,                    -1/capacC * (htCO + htAC + htBC)];
B = [1 / capacA, 0, 0; 0, 1/capacB, 0; 0, 0, 1/capacC];
C = [1, 0, 0; 0, 1, 0; 0, 0, 1];
D = [0, 0, 0; 0, 0, 0; 0, 0, 0];

Q = 100*eye(3);
R = eye(3);

[K, P, E] = lqr(A, B, Q, R);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run simulation with full-rank B matrix            %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x0 = [0;0;0];
x_des = [20;20;20];
u_des = -pinv(A\B)*x_des;

out = sim('closed_loop_control_sim.slx');

y = out.xout;
u = out.u';
t = out.t;

figure(1)
animate_input_y(y, t, u, x_des);

figure(2)
hold off; plot(t, y(:,1)); hold on; plot(t, y(:,2)); plot(t, y(:,3))
legend('Room 1', 'Room 2', 'Room 3')
title("Temperature in Rooms over Time"); xlabel("Time (sec)"); ylabel("Temp in Room (°C)")
saveas(gcf, "figs/cl_fullB_temp.png")

figure(3)
hold off; plot(t, u(1,:)); hold on; plot(t, u(2,:)); plot(t, u(3,:))
legend('Room 1', 'Room 2', 'Room 3')
title("Heating Inputs over Time"); xlabel("Time (sec)"); ylabel("Temp in Room (°C)")
saveas(gcf, "figs/cl_fullB_pow.png")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run simulation with no heater in room 3           %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
B = [1 / capacA, 0, 0; 0, 1/capacB, 0; 0, 0, 0];
[K, P, E] = lqr(A, B, Q, R);

x0 = [0;0;0];
x_des = [20;20;20];
u_des = -pinv(A\B)*x_des;

out = sim('closed_loop_control_sim.slx');

y = out.xout;
u = out.u';
t = out.t;

figure(4)
animate_input_y(y, t, u, x_des);

figure(5)
hold off; plot(t, y(:,1)); hold on; plot(t, y(:,2)); plot(t, y(:,3))
legend('Room 1', 'Room 2', 'Room 3')
title("Temperature in Rooms over Time"); xlabel("Time (sec)"); ylabel("Temp in Room (°C)")
saveas(gcf, "figs/cl_no3_temp.png")

figure(6)
hold off; plot(t, u(1,:)); hold on; plot(t, u(2,:)); plot(t, u(3,:))
legend('Room 1', 'Room 2', 'Room 3')
title("Heating Inputs over Time"); xlabel("Time (sec)"); ylabel("Temp in Room (°C)")
saveas(gcf, "figs/cl_no3_pow.png")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run simulation with no heater in room 2           %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
B = [1 / capacA, 0, 0; 0, 0, 0; 0, 0, 1/capacC];
[K, P, E] = lqr(A, B, Q, R);

x0 = [0;0;0];
x_des = [20;20;20];
u_des = -pinv(A\B)*x_des;

out = sim('closed_loop_control_sim.slx');

y = out.xout;
u = out.u';
t = out.t;

figure(4)
animate_input_y(y, t, u, x_des);

figure(5)
hold off; plot(t, y(:,1)); hold on; plot(t, y(:,2)); plot(t, y(:,3))
legend('Room 1', 'Room 2', 'Room 3')
title("Temperature in Rooms over Time"); xlabel("Time (sec)"); ylabel("Temp in Room (°C)")
saveas(gcf, "figs/cl_no2_temp.png")

figure(6)
hold off; plot(t, u(1,:)); hold on; plot(t, u(2,:)); plot(t, u(3,:))
legend('Room 1', 'Room 2', 'Room 3')
title("Heating Inputs over Time"); xlabel("Time (sec)"); ylabel("Temp in Room (°C)")
saveas(gcf, "figs/cl_no2_pow.png")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run simulation with only thermostat in room 1     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
B = [1 / capacA, 0, 0; 0, 1/capacB, 0; 0, 0, 1/capacC];
C = [1, 0, 0; 1, 0, 0; 1, 0, 0];
[K, P, E] = lqr(A, B, Q, R);
[kalmf, L, P] = kalman(ss(A, B, C, D), eye(3), eye(3));

x0 = [0;0;0];
x_des = [20;20;20];
u_des = -pinv(A\B)*x_des;

out = sim('closed_loop_control_observer_sim_2022b.slx');

y = out.yout;
xout = out.xout;
u = out.u';
t = out.t;

figure(7)
animate_input_y(xout, t, u, x_des);

figure(8)
hold off; plot(t, y(:,1)); hold on; plot(t, y(:,2)); plot(t, y(:,3))
legend('Room 1', 'Room 2', 'Room 3')
title("Temperature in Rooms over Time"); xlabel("Time (sec)"); ylabel("Temp in Room (°C)")
saveas(gcf, "figs/cl_obs1_temp.png")

figure(9)
hold off; plot(t, u(1,:)); hold on; plot(t, u(2,:)); plot(t, u(3,:))
legend('Room 1', 'Room 2', 'Room 3')
title("Heating Inputs over Time"); xlabel("Time (sec)"); ylabel("Temp in Room (°C)")
saveas(gcf, "figs/cl_obs1_pow.png")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Helper functions                                  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = animate_input_y(y, t, u, x_des) 
    n_pts = size(y,1);
    disp(n_pts)
    yt = y';
    colorbar('Ticks',[0, 1], 'TickLabels',{num2str(min(y, [], "all")), num2str(max(y, [], "all"))});

    max_t = max(max(y(:, 1:3)));
    energy_consumed = 0;
    for i = 1:ceil(n_pts/100):n_pts
        if i ~= 1
            energy_consumed = energy_consumed + sum(u(:,i).^2 * (t(i) - t(iprev)));
        end
        plot_rooms([temp2color(y(i, 1), max_t); temp2color(y(i, 2), max_t); temp2color(y(i, 3), max_t)])
        %annotation('textbox',[.2 .5 .3 .3],'String',strcat('t = ' + num2str(t(i))),'FitBoxToText','on');
%         xlabel(strcat(["Temp in 1 = ";"Temp in 2 = ";"Temp in 3 = "], string(yt(1:3, i))))
        xlabel(strcat("Energy Consumed: ", formattedDisplayText(energy_consumed, 'NumericFormat', 'shortEng')));
        title(strcat(num2str(t(i)), " seconds elapsed"));
        pause(0.001)
        iprev = i;
    end

    tdiffs = t(2:end) - t(1:(end-1));
    fprintf("Time diff: %1.9f\n", t(2)-t(1))
    disp(strcat("Total Input:", num2str(sum(u(:,2:end).^2 * tdiffs), '%0.4g')))
    disp(strcat("Max Overshoot:", num2str(max(y'-x_des*ones(1,n_pts), [],'all'), '%0.4g')))
    %disp(strcat("Max Undershoot:", num2str(min(y(n_conv+5000:n_pts,:)-ones(n_pts-(n_conv+5000) + 1,1)*x_des', [],'all'), '%0.4g')))
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
    idx = round(255 * max(0,t) / t_max) + 1;
    colors = parula(256);
    c = colors(idx, :);
end