clc
clear 
close all
%% Parameters
dt = 0.00020;
epsilom2 = 0.01;
epsilom = 10;
mu =1.0;




%P1 (5.7857, 1.5892)  P2 ( 0, 0 )  P3 (7, 0 )

% Temporal length 2, 4,4
x1 = 7;
x2 = 7;
x3 = 0;
x4 = 0;


y1= 7.0;
y2 = 0.0;
y3 = 0.0;
y4 = 7.0;

tf = 11;
t = 0:dt:tf;
Y = zeros(8, length(t));
Z = zeros(4, length(t));

k_time = 0;

%% Let check to find out initial r lengths 

P3_P2 = sqrt((x3 -x2)^2  + (y3-y2)^2);
P4_P3 = sqrt((x4 -x3)^2 + (y4-y3)^2);
P2_P1 = sqrt((x2 -x1)^2 + (y2-y1)^2);
P1_P4 = sqrt((x1 -x4)^2 + (y1-y4)^2);


% disp ("initial lengths are");
% disp (P3_P2);
% disp (P4_P3);
% disp (P2_P1);
% disp (P1_P4);


% Y(:,1) = [0.0 ; 4.8038; -  5.1962 ;17.8551;0.0 ; 0.0];
Y(:,1) = [x1 ; x2; x3 ;x4; y1; y2 ; y3; y4];
%% Euler Test
for i = 1: length(t)-1

%find the side lengths    
R1 = sqrt((x3 -x2)^2  + (y3-y2)^2);
R2 = sqrt((x4 -x3)^2 + (y4-y3)^2);
R3 = sqrt((x2 -x1)^2 + (y2-y1)^2);
R4 = sqrt((x1 -x4)^2 + (y1-y4)^2);
%stop process when at least one side is less than epsilom
%  if (R1 || R2 || R3 || R4 > epsilom)   
     
     
%compute you instantaneous movements 

x1dt = (x2 - x1)/(sqrt ((x2 -x1)^2 + (y2-y1)^2));
	
y1dt = (y2 - y1)/(sqrt ((x2 -x1)^2 + (y2-y1)^2));

x2dt = (x3 - x2)/(sqrt ((x3 -x2)^2 + (y3-y2)^2));

y2dt = (y3 - y2)/(sqrt ((x3 -x2)^2 + (y3-y2)^2));

x3dt = (x4 - x3)/sqrt((x4 -x3)^2 + (y4-y3)^2);

y3dt = (y4 - y3)/sqrt((x4 -x3)^2 + (y4-y3)^2);

x4dt = (x1 - x4)/(sqrt((x1 -x4)^2 + (y1-y4)^2));

y4dt = (y1 - y4)/(sqrt((x1 -x4)^2 + (y1-y4)^2));

%          disp('test');
%                 x = [R1, R2,R3, R4]; 




 % populate respective incremental changes to the array below
 
 dY = [x1dt; x2dt; x3dt; x4dt; y1dt; y2dt; y3dt; y4dt];

 % multiply array above by dt and add to previous value and resulting sum
 % to next array respectively
 Y(:, i+1) = Y(:, i)+ dt *dY;
    

% operation to add resulting sum above to next array respectively

x1 =Y(1,i+1);
x2 = Y(2,i+1);
x3 = Y(3,i+1);
x4 = Y(4,i+1);
y1 = Y(5,i+1);
y2 = Y(6,i+1); 
y3 = Y(7,i+1);
y4 = Y(8,i+1);

%  else
%      break;
k_time = i;
%  else
%      break;
%  end


end
% k_time = 90000;
%trim the array that contain the shrinking motion by removing empty row and
%colums
Y(:,k_time+1:end)= [];
Z(:,k_time+1 :end)= [];

 figure (1);

n =1326;
%plot triangle for every nth loop but plot first 20 points
for i = 1: k_time
   
R1 = sqrt((Y(3,i) -Y(2,i))^2  + (Y(7,i)-Y(6,i))^2);
R2 = sqrt((Y(4,i)-Y(3,i))^2 + (Y(8,i)-Y(7,i))^2);
R3 = sqrt((Y(2,i) -Y(1,i))^2 + (Y(6,i)-Y(5,i))^2);
R4 = sqrt((Y(1,i) -Y(4,i))^2 + (Y(5,i)-Y(8,i))^2);

Z(1,i) = R1;
Z(2,i) = R2;
Z(3,i) = R3;
Z(4,i) = R4;


% if i is a multiple of n
 if ( mod(i,n)==0|| (i==1))
       
%          disp('test');
%                 x = [R1, R2,R3, R4];   
                
                
% also (nested 'if') if at least one of the sides is less than epsilom      
    if((R1< epsilom)&&(R2< epsilom)&&(R3< epsilom)&&(R4< epsilom))
       x = [R1, R2,R3, R4];   
 disp(x)
 

   %Lets make our triangle!    
%Let extract the x axis for P1 to P3    
D = Y(1:4,i );
%Lets add P1 a second time to complete the loop for a triangle
Dcom = [D ;Y(1,i )];

%repeating the same steps above
U = Y(5:8,i);
Ucom = [U ;Y(5,i )];

%plot triangle just one loop
% plot(Dcom,Ucom, 'k', 'linewidth',2);
plot(Dcom,Ucom);
hold on;  
daspect([1 1 1]);
% legend ('r1', 'r2', 'r3');
% view([90 90])
% xlim([6 8])
% ylim ([3 5])

     
     
%  also (double nested 'if') if all sides are less the epsipom draw a bold
%  point at the location (which should be the midpoint)
 
 if ((R1< epsilom2)&&(R2< epsilom2)&&(R3< epsilom2)&&(R4< epsilom2))
 
%Lets make our triangle!    
%Let extract the x axis for P1 to P3    
D = Y(1:4,i );
%Lets add P1 a second time to complete the loop for a triangle
Dcom = [D ;Y(1,i )];

%repeating the same steps above
U = Y(5:8,i);
Ucom = [U ;Y(5,i )];
% daspect([1 1 1]);

%plot triangle just one loop
plot(Dcom,Ucom, '.k', 'MarkerSize', 4.9);
daspect([1 1 1]);
% hold on 
% xlim([0 14])
% ylim ([0 10])
%  
% % legend ('r1', 'r2', 'r3');
% view([90 90])
% hold on
% disp (i);


   end
 
 end 
 end
% hold off
end
% xlim([0 20])
% ylim([0 20])
% hold off
