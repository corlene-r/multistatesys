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

htAO = h_o * 11 * 3;  htBO = h_o * 9 * 3; htCO = h_o * 0 * 3;
htAB = h_i * 3 * 3;  htAC = h_i * 2 * 3; htBC = h_i * 4 * 3;

A = 100*[-1/capacA * (htAO + htAB + htAC), htAB / capacA,                    htAC / capacA;
     htAB / capacB,                    -1/capacB * (htBO + htAB + htBC), htBC / capacB;
     htAC / capacC,                    htBC / capacC,                    -1/capacC * (htCO + htAC + htBC)];
B = [1 / capacA, 0, 0; 0, 1/capacB, 0; 0, 0, 1/capacC];
C = [1, 0, 0; 0, 1, 0; 0, 0, 1];
D = [0, 0, 0; 0, 0, 0; 0, 0, 0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define u based on linear model of system          %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
heat_sys = ss(A, B, C, D);

t = 0:0.01:2;
npts = size(t, 2);
T = 2;
Wr = gram(heat_sys, 'c', gramOptions('TimeIntervals', [0, T]));
x_des = [20; 20; 20];
x0 = [0,0,0]; 
eta = Wr \ (x_des - expm(A * T) * x0');

u = zeros(3, 201);
for i=1:npts
    u(:,i) = B'*expm(A'*(T - t(i))) * eta;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Model an A(t) using a discrete system             %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disc_sys = c2d(heat_sys, t(2), c2dOptions('Method','tustin','FractDelayApproxOrder',3)); % NOTE: Assumes t starts at 0
Ad = disc_sys.A; Bd = disc_sys.B; Cd = disc_sys.C; Dd = disc_sys.D;

nonlinear_weighting = 1; % in [0, 1]; higher if want the nonlinear components to have larger influence

ys = zeros(3, npts);
xs = zeros(3, npts);
xs(:, 1) = x0';
ys(:, 1) = Cd * xs(:, 1) + Dd * u(:, 1);
for i=2:npts
    % Reset the h parameters based on temperature differences between rooms and outside
    htAO = 2.2 * (xs(1, i-1))^0.22 * 11 * h_o;  
    htBO = 2.2 * xs(2, i-1)^0.22 * 9 * h_o; 
    htCO = 2.2 * xs(3, i-1)^0.22 * 0 * h_o;
    htAB = 2.2 * abs(xs(1, i-1) - xs(2, i-1))^0.22 * 3 * h_i;  
    htAC = 2.2 * abs(xs(1, i-1) - xs(3, i-1))^0.22 * 2 * h_i; 
    htBC = 2.2 * abs(xs(2, i-1) - xs(3, i-1))^0.22 * 4 * h_i;
    
    % Set up A in reference to these new parameters
    Adis = 100*[-1/capacA * (htAO + htAB + htAC), htAB / capacA,                    htAC / capacA;
                htAB / capacB,                    -1/capacB * (htBO + htAB + htBC), htBC / capacB;
                htAC / capacC,                    htBC / capacC,                    -1/capacC * (htCO + htAC + htBC)];
    Adis = expm(Adis * t(2));

    % Use Adis-Ad in case if want 
    xs(:, i) = Ad * xs(:, i-1) + Bd * u(:, i) + nonlinear_weighting * (Adis - Ad) * xs(:, i-1);
    ys(:, i) = Cd * xs(:, i) + Dd * u(:, i);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get plots                                         %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1);
animate_input_disc(x_des, ys, u, t)

ys = ys';
figure(2)
hold off; plot(t, ys(:,1)); hold on; plot(t, ys(:,2)); plot(t, ys(:,3))
legend('Room 1', 'Room 2', 'Room 3')
title("Temperature in Rooms over Time with Full Heaters"); xlabel("Time (sec)"); ylabel("Temp in Room (°C)")
saveas(gcf, "figs/ol_AofT_temp.png")

figure(3)
hold off; plot(t, u(1,:)); hold on; plot(t, u(2,:)); plot(t, u(3,:))
legend('Room 1', 'Room 2', 'Room 3')
title("Heating Inputs over Time with Full Heaters"); xlabel("Time (sec)"); ylabel("Heating Input");
saveas(gcf, "figs/ol_AofT_pow.png")


function [] = animate_input_disc(x_des, ys, u, t) 
    y = ys'; 
    max_t = max(max(y(:, 1:3)));
    colorbar('Ticks',[0, 1], 'TickLabels',{num2str(min(y, [], "all")), num2str(max(y, [], "all"))});
    n_pts = size(t, 2);
    energy_consumed = 0;
    for i = [1:floor(n_pts/100):n_pts, n_pts]
        if i ~= 1
            energy_consumed = energy_consumed + sum(u(:,i).^2 * (t(i) - t(iprev)));
        end
        plot_rooms([temp2color(y(i, 1), max_t); temp2color(y(i, 2), max_t); temp2color(y(i, 3), max_t)])
        xlabel(strcat("Energy Consumed: ", formattedDisplayText(energy_consumed, 'NumericFormat', 'shortEng')));
        title(strcat(num2str(t(i)), " seconds elapsed"));
        pause(0.01)
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
    set(gca,'xtick',[])
    set(gca,'ytick',[])
end

function [c] = temp2color(t, t_max)
    idx = ceil(255 * max(0,t) / t_max) + 1;
    colors = parula(256);
    c = colors(min(256, idx), :);
end