clc
clear 
close all
%% Parameters
% make sure dt is smaller that epsilom and rho
dt = 0.03;
%epsilom was chosen at 0.04 because at mu =0, this was the numerical
%adjustment needed to ensure capture, rho as an additional variable whose
%single job was to tell when evader was nearing the midpoint in order to
%switch speed direction. Epsilom purpose is to create a big enough room for
%rho. rho MUST be less than epsilom
epsilom =  0.09;
rho = 0.04;
 capture_value = 5*sqrt(3) - 5/sqrt(3);
 % R records capture time, capure range, and mu for each loop
 R = zeros(3, length(25));
R_counter = 1;
%capture range
L = 0.4;

%Evader's initial radius
encirclement = zeros(2, 50);
evader_r = 0.3;

%value for generating lenght of array
 tf = 40;
 t = 0:dt:tf;
    % B is a list of ones used for the plot as a threashold for capture, the
% first array contains 1 while the second contains dt. 
    B = ones(2, length(t));
%% Euler Test
%reminder to change mu in both areas- also see line 70
% mu set to 0.0 gives an NAN error so we need mu as close to 0 as possible
for mu = 0.8;
% for mu = 0.0:0.2:2.0
%    for i = 1: length(t)-1

 sel_cone =0;
 beta_1 =0;
 beta_2 = 0;
 beta_3 = 0;
% create a circle at given dimmension and radius and return a point on the circle

% [evader_x, evader_y, encirclement] = drawCircle(evader_r,0,5 * sqrt(3)/3);
[evader_x, evader_y, encirclement] = drawCircle(evader_r,0,5/sqrt(3));
% circleoute = circle(encirclement(1,:), encirclement(2,:), evader_r) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for capture_range = capture_value * L
%     for capture_range = 1.9799    3.9598    5.9397    7.9196
%     for capture_range = capture_value * 0.2 : + capture_value * 0.2 : capture_value * 0.8
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% disp(capture_range);


Y = zeros(8, length(t));
    D = zeros(3, length(t));
    Z = zeros(1, length(t));
    xM = zeros(1, length(t));
    yM = zeros(1, length(t));
    xm=0;
    ym=0;

     
% New point generated on the circle with raduius r becomes x5 and y5

%     x1 =8.216;
%     x2 = 16.432;
%     x3 = 0.0;
        
    x1 =0;
    x2 = 5;
    x3 = -5;
%     x4 = -7.0;
    x4= evader_x;
%     x5
% 
%     y1=18.235;
%     y2 = 0.0;
%     y3 = 0.0;
    
    y1=5 *sqrt(3) ;
    y2 = 0.0;
    y3 = 0.0;
%     y4 =-7.0;
    y4 = evader_y;

Y(:,1) = [x1 ; x2 ; x3 ; x4; y1; y2 ; y3; y4];



   count = 0;
       % capture variable
  captured = false;
%           disp("capture turned off");
 reverse = false;

 count_capture = 0;
% %  mu = 0.8;
% disp(mu);
    A = [(x4-x1)^2 + (y4 - y1)^2 (x4-x2)^2 + (y4 - y2)^2 (x4-x3)^2 + (y4 - y3)^2];
    D(1,1) = sqrt(min (A))/capture_range;
i=1;
final_tf = i;
%    disp(D(1,i ));
% while evader lives 
while  D(1,i )>  1  && i < 3000
%     for i = 1: length(t)-1
        D(2,i) = i * dt; 
     %if for initial variable for one time only then the iteration takes
     %over
%      disp(D(1,i ));
     
%           disp(i);
    if i == 1
    A = [(x4-x1)^2 + (y4 - y1)^2 (x4-x2)^2 + (y4 - y2)^2 (x4-x3)^2 + (y4 - y3)^2];
    D(1,1) = sqrt(min (A))/capture_range;
    D(2,1) = 1 * dt;  
%     D(3,1) = sqrt(min (A));

    %compute to see if evader has passed the midpoint (when z < 0)
%     Z(1,i) = (1/(x1 -x3)) *((y1-y3)* (x4) + (y3 *x1) - (y1 *x3))- y4;
%      if i == 1

    
    %select cone with largest angle
%     disp( beta_1);
%     disp( 'beta_1');
%     disp( beta_2);
%     disp( 'beta_2');
%     disp( beta_3);
%     disp( 'beta_3'); 


%         disp ('selected code')
%         disp(sel_cone);
%          disp ('beta_1')
%         disp(beta_1);
%         disp ('beta_2')
%         disp(beta_2);
%         disp ('beta_3')
%         disp(beta_3);
%     switch (sel_cone)
%        case beta_1
%        %find the midpoint
%             xM(1,i) = (x1 + x2 )/2;
%             yM(1,i) = (y1 + y2)/2;
% %             disp ('choosing x1 and x2');
%             xm = xM(1,i);
%             ym = yM(1,i); 
%              Z(1,i) = (1/(x1 -x2)) *((y1-y2)* (x4) + (y2 *x1) - (y1 *x2))- y4;
%        case beta_2
%             xM(1,i) = (x2 + x3 )/2;
%             yM(1,i) = (y2 + y3)/2;
% %             disp ('choosing x2 and x3');
%             xm = xM(1,i);
%             ym = yM(1,i); 
%              Z(1,i) = (1/(x2 -x3)) *((y2-y3)* (x4) + (y3 *x2) - (y2 *x3))- y4;
%        case beta_3
%             xM(1,i) = (x1 + x3 )/2;
%             yM(1,i) = (y1 + y3)/2;
% %             disp ('choosing x1 and x3');
%             xm = xM(1,i);
%             ym = yM(1,i);  
% %             disp((x1 + x3 )/2);
% %             disp((y1 + y3)/2);
%             Z(1,i) = (1/(x1 -x3)) *((y1-y3)* (x4) + (y3 *x1) - (y1 *x3))- y4;
%        otherwise
%            disp('error with selecting midpoint');
%     end
   
       
    end
    
        %compute Escape cone
    dist_E_P1 = sqrt ((x4 -x1)^2 + (y4-y1)^2);
    dist_E_P2 = sqrt ((x4 -x2)^2 + (y4-y2)^2);
    dist_E_P3 = sqrt ((x4 -x3)^2 + (y4-y3)^2);
    dist_P1_P2 = sqrt((x2 -x1)^2 + (y2-y1)^2);
    dist_P1_P3 = sqrt ((x3 -x1)^2 + (y3-y1)^2);
    dist_P3_P2 = sqrt ((x2 -x3)^2 + (y2-y3)^2);
    
    %using the law of cosine, compute escape cone angles
    beta_1 = acos((dist_E_P1^2 + dist_E_P2 ^2 - dist_P1_P2^2)/(2 * dist_E_P1 * dist_E_P2));
    beta_2 = acos((dist_E_P2^2 + dist_E_P3 ^2 - dist_P3_P2^2)/(2 * dist_E_P2 * dist_E_P3));
    beta_3 = acos((dist_E_P1^2 + dist_E_P3 ^2 - dist_P1_P3^2)/(2 * dist_E_P1 * dist_E_P3));
    %at first iteration compute the escape cone and pick which midpoint to
 %pursue
   beta = [ beta_1 beta_2 beta_3];
   sel_cone = max(beta); 
%       disp(sel_cone) ;  
    switch (sel_cone)
       case beta_1
       %find the midpoint
            xM(1,i) = (x1 + x2 )/2;
            yM(1,i) = (y1 + y2)/2;
%             disp ('continuing x1 and x2');
            xm = xM(1,i);
            ym = yM(1,i);
%             disp(xm);
%             disp(ym);
%             Z(1,i) = (1/(x1 -x2)) *((y1-y2)* (x4) + (y2 *x1) - (y1 *x2))- y4;
       case beta_2
            xM(1,i) = (x2 + x3 )/2;
            yM(1,i) = (y2 + y3)/2;
%             disp ('continuing x2 and x3');
            xm = xM(1,i);
            ym = yM(1,i);
%             disp(xm);
%             disp(ym);
%             Z(1,i) = (1/(x2 -x3)) *((y2-y3)* (x4) + (y3 *x2) - (y2 *x3))- y4;
       case beta_3
            xM(1,i) = (x1 + x3 )/2;
            yM(1,i) = (y1 + y3)/2;
%             disp ('continuing x1 and x3');
            xm = xM(1,i);
            ym = yM(1,i);
%             disp(xm);
%             disp(ym);
%             disp((x1 + x3 )/2);
%             disp((y1 + y3)/2);
%             Z(1,i) = (1/(x1 -x3)) *((y1-y3)* (x4) + (y3 *x1) - (y1 *x3))- y4;
       otherwise
           disp('error with selecting midpoint');
   end
%     disp(D(1,i ));


%   Only run while evader has not been captured
        if D(1,i )>  1 - epsilom
            
            
%         disp ('boolean before capture loop' );
%         disp (captured);
    %if capture has not occured
   % always check for capture range that takes the longest and plot the
 % threashold line with its dt. Check the R array to see the respective durations 
%      if capture_range == 1.15470053837925
  if capture_range == capture_value * 0.6 
    B(2,i) = i * dt;
%     disp(i);
  end
    
    
   end 
% This stores the value of dt for the shortest capture range because it
% takes the longest time and is need to plot the dashed line that extends
% all the way through the plots

    
    
    if i < 10000

%if capture occurs set boolean to true        
     if D(1,i )<=  1
%         disp(D(1,i));
%         disp (captured);
        %ensure capture only get called once per loop
        count_capture = count_capture + 1;
        if count_capture <= 1
           captured = true;
%            disp("capture turned on");
%            disp(i);
%            disp ('L value inside capture loop' );
%            disp(D(2,i));
            R(1, R_counter) = i * dt;
            R(2, R_counter) = capture_range;
            R(3, R_counter) = mu;
%                disp (R_counter);
            R_counter = R_counter + 1;
                %                       disp(mu);
%            disp ('capture range');
%            disp (capture_range);
        end

    end
    %compute you instantaneous movements 

x2dt = (x3 - x2)/(sqrt ((x3 -x2)^2 + (y3-y2)^2));

y2dt = (y3 - y2)/(sqrt ((x3 -x2)^2 + (y3-y2)^2));
	
x1dt = (x2 - x1)/(sqrt ((x2 -x1)^2 + (y2-y1)^2));
	
y1dt = (y2 - y1)/(sqrt ((x2 -x1)^2 + (y2-y1)^2));

x3dt = (x1 - x3)/(sqrt ((x3 -x1)^2 + (y3-y1)^2));

y3dt = (y1 - y3)/(sqrt ((x3 -x1)^2 + (y3-y1)^2));
    




   




    if i > 1
    %    disp (Z(1,i-1));
        if Z(1,i) <= 0 && count == 0
            reverse = true;
            count = count + 1;

    %         disp (i);
        end

    end

%reverse midpoint if not already reversed
    if (reverse == true) && (count == 1)
        mu = -1.0 * mu;
        count = count + 1;


  disp (i);
    end

    
    %change midpoint as you appproach, make sure to deduct epsilom from 0
    %to get a more smoother change in heading. 

    if  Z(1,i)<= 0 - rho


      xM(1,i)= (x1 + x2 + x3)/3;
      yM(1,i)= (y1 + y2+ y3)/3;

        xm = xM(1,i);
        ym = yM(1,i);
%          disp(i);
    %      if i == 114
    %          disp(xm);
    %          disp(ym);
    %      end
    end

%     disp(xm);
%     disp(ym);
    x4dt =( mu *(xm - x4))/(sqrt ((xm -x4)^2 + (ym-y4)^2));

    y4dt = (mu * (ym - y4))/(sqrt ((xm -x4)^2 + (ym-y4)^2));






     % populate respective incremental changes to the array below

     dY = [x1dt; x2dt; x3dt; x4dt; y1dt; y2dt; y3dt; y4dt];

     % multiply array above by dt and add to previous value and resulting sum
     % to next array respectively
     Y(:, i+1) = Y(:, i)+ dt *dY;


    % operation to add resulting sum above to next array respectively

    x1 =Y(1,i+1);
    x2 = Y(2,i+1);
    x3 = Y(3,i+1);
    x4 = Y(4, i+1);
%     x5 = Y(5,i+1);
    y1 = Y(5,i+1);
    y2 = Y(6,i+1); 
    y3 = Y(7,i+1);
    y4 = Y(8,i+1); 
%     y5 = Y(10,i+1);

    A = [(x4-x1)^2 + (y4 - y1)^2 (x4-x2)^2 + (y4 - y2)^2 (x4-x3)^2 + (y4 - y3)^2];
    D(1,i+1) = sqrt(min (A))/capture_range;
    D(2,i+1) = (i+1) * dt;  
    D(3,i+1) = sqrt(min (A));
%     Z(1,i+1) = (1/(x1 -x3)) *((y1-y3)* (x4) + (y3 *x1) - (y1 *x3))- y4;

    switch (sel_cone)
       case beta_1

            Z(1,i+1) = (1/(x1 -x2)) *((y1-y2)* (x4) + (y2 *x1) - (y1 *x2))- y4;
       case beta_2

            Z(1,i+1) = (1/(x2 -x3)) *((y2-y3)* (x4) + (y3 *x2) - (y2 *x3))- y4;
       case beta_3

            Z(1,i+1) = (1/(x1 -x3)) *((y1-y3)* (x4) + (y3 *x1) - (y1 *x3))- y4;
       otherwise
           disp('error with selecting midpoint');
   end

    %     disp('x')
    %     disp (x5);
    %     disp ('y')
    %     disp (y5);



    % D(2,i) = dt * i;
    % 
    % %plot function later on
    % D(3,i) = sqrt((dt)* i)/(D(1,i));

    % disp((Z(i)));
    end
    

final_tf = i;
i= i+1;

end;

% disp(D(1,i ));
% disp(i);


%reset mu to +ve after while loop is done
if mu < 0
   mu = -mu; 
%    disp(i)
end
%Trim array by elimination unused arrays positions
D(:,final_tf:end)= [];
%  D = D(1: end-final_tf,:);
%  D = D(2: end-final_tf,:);
% disp(final_tf);
    %  if capture has occured, plot graph
%                      if captured == true
%                         figure(2) 
% 
%                         plot(D(2,:),D(1,: ));
%     %                     plot(D(2,:),a(:,:,count-1));
%                         xlim([0 30])
%                         ylim ([0 10])
%                         title('mu = 0.7, r=0.1')
%                         xlabel('t') 
%                         ylabel('sqrt(d(t))/l') 
%                         % figure(1)
%                         % plot (Y(1,:), Y(3,:), 'b');
%                         hold on
% %                         disp ('plotting inside')
%                         % plot (Y(2,:), Y(4,:), 'r');
%                          captured = false;
% %                          disp("capture turned off");
%                          %decrement count back to zero
%                           
% %                     plot(D(2,:),B(1,:), 'b--');
%                      end
                        
    end
end
   
    plot(B(2,:),B(1,:), 'b--');
    
%     %3D Plot
%     figure (3)
%     scatter3 (abs (R(3,:)),R(2,:), R(1,:) );
%     zlabel ('Capture time, tf');
%     ylabel(' Pursuer capture ranges');
%     xlabel ('Pursuer speed');

% end
%  legend( 'l=60% initial capture range','l=80%  initial capture range' , 'sqrt(d(t))/l = 1');
%'l=20% initial capture range','l=40% initial capture range',
k = 1;
% %uncomment this section for animation plot
% %uncomment line below to have the code print immediately once capture range
% %is attainded
while ((D(1,k)>1) && k < 200)
%while k < 150 
 figure(1)

%plot(Y(1, k), Y(5, k), 'o', 'MarkerSize', 10);
circleouta = circle(Y(1, k ), Y(5, k),  capture_value * L , 'b') ;
hold on
plot (Y(1, k), Y(5, k), '.b');
% daspect([1 1 1])
xlim([-10 10])
ylim ([-5 15])
%plot(Y(2,k), Y(6,k), 'o', 'MarkerSize', 10);
plot (Y(2,k), Y(6,k), '.b');
circleoutb = circle(Y(2,k), Y(6,k),  capture_value * L , 'b') ;

plot (Y(3,k), Y(7,k), '.b');
circleoutc = circle(Y(3,k), Y(7,k),  capture_value * L, 'b') ;
%plot(Y(3,k), Y(7,k), 'o', 'MarkerSize', 10);

% plot (Y(4,k), Y(9,k), 'go');
% circleoutd = circle(Y(4,k), Y(9,k), capture_value * 0.1 , 'y') ;
plot(encirclement(1,:), encirclement(2,:), 'r');
% xlim([-20 20])
% ylim ([-20 20])

plot (Y(4,k), Y(8,k), '.r');
k = k + 1;
% disp(k);
hold off;

legend(  'Equilateral triangle formation of side length 10, r = 0.3 L = 0.4 max capture range, speed \mu = 0.8 ');
end

% % function that draw circle and fills the circle
function circles = circle(x,y,r,c)
%hold on
th = 0:pi/50:2*pi;
x_circle = r * cos(th) + x;
y_circle = r * sin(th) + y;
circles = plot(x_circle, y_circle);
% fill(x_circle, y_circle, c)
%hold off
axis equal
end
 
% %function that draws circle and return randox x,y point on circumference
function [u, w,m] = drawCircle(r,x,y)
% count = 0;
th = 0:pi/60:2*pi;

encirclement(1,:) = r * cos(th) + x;
encirclement(2,:) = r * sin(th) + y;
% encirclement(2,:)
% circles = plot(q, r);
r=70;
%uncomment r below to randomize point on circle
% r = randi([0 50],1,1);
u =encirclement(1,r);
w =encirclement(2,r);
m =  encirclement(:,:);
% fill(q, z, r);
% disp(r);
end

% figure(1)
% plot (Y(1,:), Y(3,:), 'b');
% hold on
% plot (Y(2,:), Y(4,:), 'r');
