clc
clear 
close all
%% Parameters
dt = 0.01;
x1 =1;
x2 =10;
y1=1;
y2 =10;
tf = 6.37;
t = 0:dt:tf;
Y = zeros(4, length(t));
D = zeros(1, length(t));

Y(:,1) = [1.0 ;10.0; 1.0; 10.0];
%% Euler Test
for i = 1: length(t)-1
    
dx2dt = (x1 - x2)/(sqrt ((x2 -x1)^2 + (y2-y1)^2));

dy2dt = (y1 - y2)/(sqrt ((x2 -x1)^2 + (y2-y1)^2));
	
dx1dt = (x2 - x1)/(sqrt ((x2 -x1)^2 + (y2-y1)^2));
	
dy1dt = (y2 - y1)/(sqrt ((x2 -x1)^2 + (y2-y1)^2));
    
    dY = [dx1dt; dx2dt; dy1dt; dy2dt];
    Y(:, i+1) = Y(:, i)+ dt *dY;
    

end


%uncomment this for animation plot
for i = 1:2800
 figure(1)
plot (Y(1, i), Y(3, i), 'bo');
hold on
xlim([0 10])
ylim ([0 10])
plot (Y(2,i), Y(4,i), 'ro');

hold off;

end

% figure(1)
% plot (Y(1,:), Y(3,:), 'b');
% hold on
% plot (Y(2,:), Y(4,:), 'r');
