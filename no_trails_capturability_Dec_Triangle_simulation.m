clc
clear
close all
%% Parameters
dt = 0.0001;
l =0.5000;





dt_shrink = 0.00002;
epsilom = 0.00009;
x1_shrink = 0;
x2_shrink = 1;
x3_shrink = -1;




y1_shrink= sqrt(3);
y2_shrink = 0.0;
y3_shrink =  0.0;

tf_shrink = 20;
t_shrink = 0:dt_shrink:tf_shrink;
Y = zeros(6, length(t_shrink)); %200001

%Escape path parameters
x1 =(x1_shrink + x2_shrink + x3_shrink)/3;
y1= 1/ sqrt(3);
mu  =1/2;
mu_pursuer = 1;
tf = 10;
t = 0:dt:tf;
W = zeros(2, length(t_shrink));
D = zeros(2, length(t_shrink));
W(:,1) = [0 ; 1/sqrt(3)];
r_o = 2;
minimum_mu = (sqrt(3)/2)* (l/(r_o - (2*l)));
t_e = (2/sqrt(3)) *(r_o /((4 * mu) + sqrt(3)));
t_er = round(t_e, 4) *10000;
% t_e = 2/3; %4/9;
t_c = (2/3) * (r_o - (2 *l));
t_cr = round(t_c, 4) *10000;
t_et = 0:dt_shrink:t_e;

k_time=0;
j = t_e/dt_shrink;

G = zeros(1, length(t_shrink));
capture_range = l;
W(:,1) = [x1, y1];
A = [(W(1, 1)-x1_shrink)^2 + (W(2, 1)- y1_shrink)^2 (W(1, 1)-x2_shrink)^2 + (W(2, 1)- y2_shrink)^2 (W(1, 1)-x3_shrink)^2 + (W(2, 1) - y3_shrink)^2];
G(1,1) = sqrt(min (A))/capture_range;
disp(G(1,1))
i=1;
final_tf = i;
 final_tf_te = i; 

th_t = 0:pi/50:2*pi;
%% Let check to find out initial r lengths

initialR1 = sqrt((x3_shrink -x2_shrink)^2  + (y3_shrink-y2_shrink)^2);
initialR2 = sqrt((x1_shrink -x3_shrink)^2 + (y1_shrink-y3_shrink)^2);
initialR3 = sqrt((x2_shrink -x1_shrink)^2 + (y2_shrink-y1_shrink)^2);




Y(:,1) = [x1_shrink ; x2_shrink; x3_shrink ; y1_shrink; y2_shrink ; y3_shrink];

%% Euler Test 2



x1dt = 1/2 * mu *(sin((1/sqrt(3))* log(1 +(sqrt(3)/(4 * mu))))    +    sqrt(3)* cos((1/sqrt(3))* log(1 + (sqrt(3)/(4 * mu)))));
	
y1dt = 1/2 * mu *(cos((1/sqrt(3))* log(1 + (sqrt(3)/(4 * mu))))    -    sqrt(3)* sin((1/sqrt(3))* log(1 + (sqrt(3)/(4 * mu)))));


%% Euler Test 1- We will simulate the cyclic pursuit motion here
% for i = 1: length(t_shrink)-1 %200001
%     while  G(1,i )>  1  && i < 90000
%     while (dt_shrink * i)< t_e + 0.0004
       while i < j+15000 && G(1,i )>  1 
        if i == 1
            A = [(W(1, 1)-x1_shrink)^2 + (W(2, 1)- y1_shrink)^2 (W(1, 1)-x2_shrink)^2 + (W(2, 1)- y2_shrink)^2 (W(1, 1)-x3_shrink)^2 + (W(2, 1) - y3_shrink)^2];
            G(1,1) = sqrt(min (A))/capture_range;
        disp('test')
        end 

    R1 = sqrt((x3_shrink -W(1, i))^2  + (y3_shrink-W(2, i))^2);
    R2 = sqrt((x1_shrink -W(1, i))^2 + (y1_shrink-W(2, i))^2);
    R3 = sqrt((x2_shrink -W(1, i))^2 + (y2_shrink-W(2, i))^2);
%     disp (R1);
%     % % plot values for average length of epsilom and above
%     if (((R1 + R2 + R3)/3) >= epsilom)
        
        %compute you instantaneous movements
        
x1_shrinkdt_shrink = ((x2_shrink - x1_shrink)/(sqrt ((x2_shrink -x1_shrink)^2 + (y2_shrink-y1_shrink)^2)))*mu_pursuer;
	
y1_shrinkdt_shrink = ((y2_shrink - y1_shrink)/(sqrt ((x2_shrink -x1_shrink)^2 + (y2_shrink-y1_shrink)^2)))*mu_pursuer;

x2_shrinkdt_shrink = ((x3_shrink - x2_shrink)/(sqrt ((x3_shrink -x2_shrink)^2 + (y3_shrink-y2_shrink)^2)))*mu_pursuer;

y2_shrinkdt_shrink = ((y3_shrink - y2_shrink)/(sqrt ((x3_shrink -x2_shrink)^2 + (y3_shrink-y2_shrink)^2)))*mu_pursuer;

x3_shrinkdt_shrink = ((x1_shrink - x3_shrink)/(sqrt ((x1_shrink -x3_shrink)^2 + (y1_shrink-y3_shrink)^2)))*mu_pursuer;

y3_shrinkdt_shrink = ((y1_shrink - y3_shrink)/(sqrt ((x1_shrink -x3_shrink)^2 + (y1_shrink-y3_shrink)^2)))*mu_pursuer;
        
        
        
        
        
        
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
        
%         k_time = i;
        %compute change in Evader location at time dt and place in an array
        dW = [x1dt; y1dt];
        
        %update Evader's location information into array 'W'
        W(:, i+1) =W(:, i) + (dt_shrink *dW);
        A = [(W(1, i)-x1_shrink)^2 + (W(2, i)- y1_shrink)^2 (W(1, i)-x2_shrink)^2 + (W(2, i)- y2_shrink)^2 (W(1, i)-x3_shrink)^2 + (W(2, i) - y3_shrink)^2];
       G(1,i+1) = sqrt(min (A))/capture_range;
        
%     end
    final_tf = i;
    if (dt_shrink * i)< t_e
       final_tf_te = i; 
    end 
i= i+1;
    
end
Y(:, final_tf  +1:end)= [];
W(:, final_tf  +1 :end)= [];
G(:, final_tf  +1 :end)= [];
% k_time = 10000;

%vary location of info on evader plot points
b = -1;
a = 1;

v =2;
n =50;
j = 4020;
k_time =  final_tf ;

 
for i = 1:final_tf
    
%     disp(i);
    % Simulating Evader's movement
    if ((mod (i,n)==0)||(i==1))
        figure(1);
        plot (W(1,i), W(2, i), '.r','linewidth',1);
        % daspect([1 1 1]);
        % end
        hold on
%         daspect([1 1 1]);
        xlim([-2 2])
        ylim ([-1 3])
        
     str = string(t_er/10000);
 msg = ' t=' + str;
%  text (Y(1,t_er),Y(4,t_er),msg);    
        
        
        
        
        
        
        
        
        
   
        %Lets make our triangle for t_er
        %Let extract the x axis for P1 to P3
        D = Y(1:3, final_tf_te  );
        %Lets add P1 a second time to complete the loop for a triangle
        Dcom = [D ;Y(1,final_tf_te )];
        
        %repeating the same steps above
        U = Y(4:6,final_tf_te);
        Ucom = [U ;Y(4,final_tf_te )]; 
        disp(j);
%         if (i == t_er)
% %         disp(t_er);
%         end
         plot(Dcom,Ucom, 'b', 'LineWidth',2);
        
        
        
        
        
        
        
        %Lets make our triangle!
        %Let extract the x axis for P1 to P3
        D = Y(1:3,i );
        %Lets add P1 a second time to complete the loop for a triangle
        Dcom = [D ;Y(1,i )];
        
        %repeating the same steps above
        U = Y(4:6,i);
        Ucom = [U ;Y(4,i )];
        
        %Lets make sure the both x and y axis are the same
        % daspect([1 1 1]);
        
        %add time steps to each plot
        if (i==1)
            str = string(0);
            % elseif(i ==   5774)
            %     str = "1/ \surd{3}";
        elseif(i == t_cr)
            str = "2/3";
        elseif(i == t_er)
            str = "4/(2\surd{3}+3)";
        else
            str = string(i *dt);
        end
        msg = ' t=' + str ;
        
        
        if (i==t_er)
            plot(Dcom,Ucom, 'LineWidth',2);
            % hold on;
            %display current time slot
            % msg1 = 't_e = 4/9';
%             text (Y(1,i),Y(4,i),msg);
        else
            plot(Dcom,Ucom, 'k');
            
            
            
        end
        
        % hold off;
        
        
        
        
        %     else
        %Add some texts to the figures
        txt1 = 'P_1';
        txt2 = 'P_3';
        txt3 = 'P_2';
        txt4 = 'O';
        
        
        
        %plot the capture range arround the pursuer at each time
        circleouta = circle(Y(1, i), Y(4, i), l, 'b') ;
        circleoutb = circle(Y(2, i), Y(5, i), l, 'b') ;
        circleoutc = circle(Y(3, i), Y(6, i), l, 'b') ;
        
        %Put time stamps on the circumference of the circles
        %compute x and y location
        x_t1 = l * cos(th_t(63)) + Y(1, i);
        y_t1 = l * sin(th_t(63)) + Y(4, i);
        x_t2 = l * cos(th_t(63)) + Y(2, i);
        y_t2 = l * sin(th_t(63)) + Y(5, i);
        x_t3 = l * cos(th_t(63)) + Y(3, i);
        y_t3 = l * sin(th_t(63)) + Y(6, i);
        
        %add time steps to each plot
        if (i==1)
            str = string(0);
            
        elseif(i == t_cr)
            str = "2/3";
        elseif(i == t_er)
            str = "4/(2\surd{3}+3)";
        else
            str = string(i *dt);
        end
        
        
        msg = ' t=' + str ;
        
        hold off
    end
    
end

% legend(  'Escape Trajectory, r_o = 2.0, l=1/2, \mu = \surd{3}/4, t_c =2/3, t^{*}_e= 4/(2\surd{3}+3)');
% legend(  'Escape Trajectory, r_o = 2.0, l= 3/4, \mu = (4/5) * \surd{3}, t_c =2/3, t^{*}_e= 4/9');
% legend(  'Escape Trajectory, r_o = 2.0, l=3/4, \mu = \surd{3}, t_c =1/3, t^{*}_e= 4/15');
%  legend(  'Escape Trajectory, r_o = 2.0, l=1/2, \mu = \surd{3}/2, t_c =2/3, t^{*}_e= 4/9');
% legend(  'Escape Trajectory, r_o = 2.0, l=1/4, \mu = \surd{3}, t_c = 1, t^{*}_e= 4/15)');
legend(  'Escape Trajectory, r_o = 2.0, l= 1/2, \mu =1/2, t_c =2/3, t^{*}_e= 4/(2 \surd{3} + 3)');
%function that draw circle
function circles = circle(x,y,r,c)
%hold on
th = 0:pi/50:2*pi;
x_circle = r * cos(th) + x;
y_circle = r * sin(th) + y;
circles = plot(x_circle, y_circle, 'b');
%fill(x_circle, y_circle, c)
%hold off
axis equal
end
