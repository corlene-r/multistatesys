%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set up parameters needed in simulation            %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

effs = ones(1, 8) * 0.9;

% Capacity numbers are how many tiles on floor plan room is
% Times 2500 to convert to capacity
capacs = [2; 22; 4; 6; 2; 13; 33; 9] * 2500;

% * 3 for height of house; * 0.3 W/m^2 * K for exterior wall heat conduc.
out_hts = [3, 8, 0, 5, 0, 5, 16, 3]' * 3 * 0.3; 

hts = [0, 3, 0, 0, 0, 0, 0, 0;
       0, 0, 6, 1, 0, 2, 4, 0;
       0, 0, 0, 1, 2, 1, 0, 0;
       0, 0, 0, 0, 1, 2, 0, 0;
       0, 0, 0, 0, 0, 3, 0, 0;
       0, 0, 0, 0, 0, 0, 0, 3;
       0, 0, 0, 0, 0, 0, 0, 6;
       0, 0, 0, 0, 0, 0, 0, 0] * 3 * 0.9;
hts = hts + hts'; 

A = hts - diag(sum(hts, 2) + out_hts);
A = A .* (1 ./ capacs);
B = diag([0 1 0 1 0 1 1 1]) .* (1 ./ capacs); 
C = eye(8);
D = zeros(8, 8);

Q = 100*eye(8);
R = eye(8);
[K, P, E] = lqr(A, B, Q, R);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run simulation with 8 rooms, not full heaters     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x0 = [0;0;0;0;0;0;0;0];
x_des = [20;20;20;20;20;20;20;20];
u_des = -pinv(A\B)*x_des;
disp("Room temps we're converging to:")
disp(-(A^(-1)) * B * u_des)

out = sim('closed_loop_8room_sim.slx');

y = out.xout;
u = out.u';
t = out.t;
figure(1)
animate_input_y(y, t, u, x_des);

t = t';
figure(2)
hold off; plot(t, y(:,1)); hold on; plot(t, y(:,2)); plot(t, y(:,3)); 
plot(t, y(:,4)); plot(t, y(:,5)); plot(t, y(:,6)); plot(t, y(:,7)); plot(t, y(:,8))
lgd = legend('Room 1', 'Room 2', 'Room 3', 'Room 4', 'Room 5', 'Room 6', 'Room 7', 'Room 8');
lgd.Location = "southeast";
title("Temperature in Rooms over Time with Full Heaters"); xlabel("Time (sec)"); ylabel("Temp in Room (°C)")
saveas(gcf, "figs/cl_8room_temp.png")

% t = t';

figure(3)
hold off; plot(t, u(1,:)); hold on; plot(t, u(2,:)); plot(t, u(3,:)); plot(t, u(3, :)); 
plot(t, u(4,:)); plot(t, u(5,:)); plot(t, u(6,:)); plot(t, u(7,:)); plot(t, u(8,:))
legend('Room 1', 'Room 2', 'Room 3', 'Room 4', 'Room 5', 'Room 6', 'Room 7', 'Room 8')
title("Heating Inputs over Time with Full Heaters"); xlabel("Time (sec)"); ylabel("Heating Input");
saveas(gcf, "figs/cl_8room_pow.png")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run simulation with 8 rooms with Less Observers   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C = diag([0 1 0 1 1 0 0 1]); 
x0 = [0;0;0;0;0;0;0;0];
x_des = [20;20;20;20;20;20;20;20];
u_des = -pinv(A\B)*x_des;
[kalmf, L, P] = kalman(ss(A, B, C, D), eye(8), eye(8));

disp("Room temps we're converging to:")
disp(-(A^(-1)) * B * u_des)

out = sim('closed_loop_8rooms_observer.slx');

y = out.yout;
xout = out.xout;
u = out.u';
t = out.t;
figure(1)
animate_input_y(xout, t, u, x_des);

t = t';
figure(2)
hold off; plot(t, xout(:,1)); hold on; plot(t, xout(:,2)); plot(t, xout(:,3)); 
plot(t, xout(:,4)); plot(t, xout(:,5)); plot(t, xout(:,6)); plot(t, xout(:,7)); plot(t, xout(:,8))
lgd = legend('Room 1', 'Room 2', 'Room 3', 'Room 4', 'Room 5', 'Room 6', 'Room 7', 'Room 8');
lgd.Location = "southeast";
title("Temperature in Rooms over Time with Missing Heaters & Thermostats"); xlabel("Time (sec)"); ylabel("Temp in Room (°C)")
saveas(gcf, "figs/cl_8room_temp.png")

% t = t';

figure(3)
hold off; plot(t, u(1,:)); hold on; plot(t, u(2,:)); plot(t, u(3,:)); plot(t, u(3, :)); 
plot(t, u(4,:)); plot(t, u(5,:)); plot(t, u(6,:)); plot(t, u(7,:)); plot(t, u(8,:))
legend('Room 1', 'Room 2', 'Room 3', 'Room 4', 'Room 5', 'Room 6', 'Room 7', 'Room 8')
title("Heating Inputs over Time with Missing Heaters & Thermostats"); xlabel("Time (sec)"); ylabel("Heating Input");
saveas(gcf, "figs/cl_8room_pow.png")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Helper Functions                                  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = animate_input_y(y, t, u, x_des) 
    n_pts = size(u, 2);
    yt = y';

    max_t = max(y, [], "all"); min_t = min(y, [], "all");
    colorbar('Ticks',[0, 1], 'TickLabels',{num2str(min(y, [], "all")), num2str(max(y, [], "all"))});
    energy_consumed = 0;
    for i = [1:max(1, floor(n_pts/100)):n_pts, n_pts]
        if i ~= 1
            energy_consumed = energy_consumed + sum(u(:,i).^2 * (t(i) - t(iprev)));
        end
        title(strcat(num2str(t(i)), " seconds elapsed"));
        colors = arrayfun(@(y) temp2color(y, max_t, min_t), y(i, :), 'uniformoutput',false);
        %disp(i)   
        plot_rooms(colors)
        %annotation('textbox',[.2 .5 .3 .3],'String',strcat('t = ' + num2str(t(i))),'FitBoxToText','on');
        % xlabel(strcat(["Temp in 1 = ";"Temp in 2 = ";"Temp in 3 = ";"Temp in 4 = ";"Temp in 5 = ";"Temp in 6 = ";"Temp in 7 = ";"Temp in 8 = "], string(yt(1:8, i))))
        xlabel(strcat("Energy Consumed: ", formattedDisplayText(energy_consumed, 'NumericFormat', 'shortEng')));
        pause(0.01)
        iprev = i;
    end

    tdiffs = t(2:end) - t(1:(end-1));
    fprintf("Time diff: %1.9f\n", t(2)-t(1))
    disp(strcat("Total Input:", num2str(sum(u(:,2:end).^2 * tdiffs), '%0.4g')))
    %disp(strcat("Max Overshoot:", num2str(max(y'-x_des*ones(1,n_pts), [],'all'), '%0.4g')))
end

function [] = plot_rooms(colors,outside_color)
    if nargin == 1
        outside_color = [0,0,0];
    end
    rectangle('Position', [-2,-1,15+2,9+1], 'FaceColor', outside_color)
    rectangle('Position', [0,3,7,4], 'FaceColor', colors{2})
    text(3.5-0.1, 5, '2')
    rectangle('Position', [0,5,1,2], 'FaceColor', colors{1})
    text(0.5-0.1, 6, '1')
    rectangle('Position', [1,3,4,1], 'FaceColor', colors{3})
    text(3-0.1,3.5,'3')
    rectangle('Position', [0,0,2,3],'FaceColor', colors{4})
    text(1-0.1,1.5,'4')
    rectangle('Position', [2,0,5,3], 'FaceColor', colors{6})
    text(4.5-0.1,1.5,'6')
    rectangle('Position', [2,2,2,1], 'FaceColor', colors{5})
    text(3-0.1,2.5,'5')
    rectangle('Position', [7,0,6,7], 'FaceColor', colors{7})
    text(10-0.1,4.5,'7')
    rectangle('Position', [7,0,3,3], 'FaceColor', colors{8})
    text(8.5-0.1,1.5,'8')
    xlim([-2,15]); ylim([-1,8])
    set(gca,'xtick',[], 'ytick', [])
end

function [c] = temp2color(t, t_max, t_min)
    idx = ceil(255 * (max(0,t) - t_min) / (t_max-t_min)) + 1;
    colors = parula(256);
    c = colors(min(idx, 256), :);
end