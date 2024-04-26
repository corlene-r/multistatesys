figure(1);
plot_rooms([1,1,1;1,1,1;1,1,1])

%% State system; I put some random values below, but I feel like something is wrong
%% Temperature will blow up over time. 
Inter_Room_Heat_Bleeds = [0.7, 0.2, 0.1;
                          0.2, 0.6, 0.2;
                          0.1, 0.2, 0.7];
Outside_Heat_Bleeds = [0.2; 0.2; 0];

t = linspace(0,5,200);
u = ones(size(t)) * 70; 
x0 = [0,0,0];
heat_sys = ss(Inter_Room_Heat_Bleeds, Outside_Heat_Bleeds, eye(3), 0);
figure(2)
lsim(heat_sys, u, t);
y = lsim(heat_sys, u, t, x0);

%% 
max_t = max(max(max(y)), max(max(u)));
figure(3)
for i = 1:max(size(t))
    %disp(i)
    if i == 167
        disp([y(i, 1), y(i, 2), y(i, 3), u(i)])
    end
    plot_rooms([temp2color(y(i, 1), max_t); temp2color(y(i, 2), max_t); temp2color(y(i, 3), max_t)], temp2color(u(i), max_t))
    pause(0.05)
end
% dotTrooms = Inter_Room_Heat_Bleeds * Trooms + Outside_Heat_Bleeds * Toutside;
% Output is heat in each room, which should be just given by identity

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
    r = [1, 0, 0];
    b = [0, 0, 1];
    c = (max(0,t) * r + (t_max - max(0,t)) * b) / t_max; 
end