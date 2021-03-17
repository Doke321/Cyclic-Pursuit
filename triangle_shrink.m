clc
clear 
close all
%% Parameters
dt_shrink = 0.00009;
epsilom = 0.0009;
mu =1.0;

% h is a horizontal shift factor for the plot. My initial values plot some
% points in the negative number region. The h value helps shift the plot in
% the positive region. 


%P1 P1 (1, 1.7321)  P2 ( 0, 0 )  P3 (2, 0 )

% Temporal length 2, 4,4
x1_shrink = 5.7857;
x2_shrink = 0;
x3_shrink = 7;


y1_shrink= 1.5892;
y2_shrink = 0.0;
y3_shrink =  0.0;

tf_shrink = 20;
t_shrink = 0:dt_shrink:tf_shrink;
Y = zeros(6, length(t_shrink));

%% Let check to find out initial r lengths 

initialR1 = sqrt((x3_shrink -x2_shrink)^2  + (y3_shrink-y2_shrink)^2);
initialR2 = sqrt((x1_shrink -x3_shrink)^2 + (y1_shrink-y3_shrink)^2);
initialR3 = sqrt((x2_shrink -x1_shrink)^2 + (y2_shrink-y1_shrink)^2);


disp ("initial lengths are");
disp (initialR1);
disp (initialR2);
disp (initialR3);


% Y(:,1) = [0.0 ; 4.8038; -  5.1962 ;17.8551;0.0 ; 0.0];
Y(:,1) = [x1_shrink ; x2_shrink; x3_shrink ; y1_shrink; y2_shrink ; y3_shrink];
%% Euler Test
for i = 1: length(t_shrink)-1

%find the side lengths    
R1 = sqrt((x3_shrink -x2_shrink)^2  + (y3_shrink-y2_shrink)^2);
R2 = sqrt((x1_shrink -x3_shrink)^2 + (y1_shrink-y3_shrink)^2);
R3 = sqrt((x2_shrink -x1_shrink)^2 + (y2_shrink-y1_shrink)^2);

%stop process when at least one side is less than epsilom
 if (R1 && R2 && R3 > epsilom)   
     
     
%compute you instantaneous movements 

x1_shrinkdt_shrink = (x2_shrink - x1_shrink)/(sqrt ((x2_shrink -x1_shrink)^2 + (y2_shrink-y1_shrink)^2));
	
y1_shrinkdt_shrink = (y2_shrink - y1_shrink)/(sqrt ((x2_shrink -x1_shrink)^2 + (y2_shrink-y1_shrink)^2));

x2_shrinkdt_shrink = (x3_shrink - x2_shrink)/(sqrt ((x3_shrink -x2_shrink)^2 + (y3_shrink-y2_shrink)^2));

y2_shrinkdt_shrink = (y3_shrink - y2_shrink)/(sqrt ((x3_shrink -x2_shrink)^2 + (y3_shrink-y2_shrink)^2));

x3_shrinkdt_shrink = (x1_shrink - x3_shrink)/(sqrt ((x1_shrink -x3_shrink)^2 + (y1_shrink-y3_shrink)^2));

y3_shrinkdt_shrink = (y1_shrink - y3_shrink)/(sqrt ((x1_shrink -x3_shrink)^2 + (y1_shrink-y3_shrink)^2));






 % populate respective incremental changes to the array below
 
 dY = [x1_shrinkdt_shrink; x2_shrinkdt_shrink; x3_shrinkdt_shrink; y1_shrinkdt_shrink; y2_shrinkdt_shrink; y3_shrinkdt_shrink];

 % multiply array above by dt_shrink and add to previous value and resulting sum
 % to next array respectively
 Y(:, i+1) = Y(:, i)+ dt_shrink *dY;
    

% operation to add resulting sum above to next array respectively

x1_shrink =Y(1,i+1);
x2_shrink = Y(2,i+1);
x3_shrink = Y(3,i+1);
y1_shrink = Y(4,i+1);
y2_shrink = Y(5,i+1); 
y3_shrink = Y(6,i+1);

 else
     break;
 end


end

figure (1)

n = 6980;
%plot triangle for every nth loop but plot first 20 points
for i = 1: length(t_shrink)-1

    if ((mod(i,n)==0)|| i<20)
%Lets make our triangle!    
%Let extract the x axis for P1 to P3    
D = Y(1:3,i );
%Lets add P1 a second time to complete the loop for a triangle
Dcom = [D ;Y(1,i )];

%repeating the same steps above
U = Y(4:6,i);
Ucom = [U ;Y(4,i )];
daspect([1 1 1]);
%plot triangle just one loop
plot(Dcom,Ucom);
% legend ('r1', 'r2', 'r3');
% view([90 90])
hold on
% disp (i);
    end
end
% xlim([0 20])
% ylim([0 20])
hold off
