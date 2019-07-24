%%%%%%%%%%%%%%%%%%%%%%      SRM Design       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%   R. Krishnan Procedure(Rev #0) %%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%% Developed power P_d equation 

 %Sigma_s  = L_as/L_au;              % L_as : Inductance defined at the deep saturated aligned operating point
 %Sigma_u = L_au/L_u;                % L_au : Inductance defined at the unsaturated aligned operating point
                                    % L_u  : Inductance defined at the unsaturated  unaligned operating point  
                                  
 P_d = 10000;                       % Developed Power in Watts
 I_max =  20;                       % Maximum amperes allowed from the supply or the source
 W_m = 523.59;                      % Base speed in rad/s 
 N_r = (60/pi)*(W_m)                % Speed in RPM (Speed is controlled by frequency of switching )
 Beta_s = 30;                       % Stator pole arc angle degrees
 time =  Beta_s /(W_m);             % Time taken for rotor to move from unaligned to alligend 
 
 m = 1;                             % Number of phases conducting simultaniously
 A_s = 50000 ;                      % Specific electric loading A_s = (2*T_ph*i*m)/(pi*D)
                                    % 25000 < A_s < 90000
                                    
 K = 0.5;                              % Non Servo applications 0.25 < K < 0.7
                                    % Servo applications 1 < K < 3
 L = k*D ;                          % Axial length of the machine or stack length taken as multiple or sub multiple of 'D'
 Theta_i =   ;                      % Current conduction angle in degrees 
 q = 3  ;                           % Number of stator phases(P_s /2)
 P_r = 6 ;                          % Number of rotor poles
 k_e = 0.95 ;                       % Efficinecy of the machine
 Theta_i =  30;                     % Current pulse angle
 K_d = (Theta_i*q*P_r)/(360) ;      % Maximum value it can take is '1'
 K_1 =  (pi^2)/120;
 K_2 = 0.5;                                     % Variable dependant on the operating point(need to fix at design stage)
                                                % 0.65 < K_2 < 0.75  (1-(1/( Sigma_s * Sigma_u )))
 K_3 = pi/4 ;
 
 D = (P_d/(K_e*K_d*K_1*K_2*K*B*A_s*N_r))^(1/3);                   
                                                % Bore diameter (rotor Pole diameter)
                                                % Calculated from P_d = K_e*K_d*K_1*K_2*B*A_s*D^2*L* N_r; % Developed power equation
 T_ph = (A_s*pi*D)/(2*I_max*m);                 % Number of Turns per phase
                                                                                                                             
                                              
 
 T_d = K_e*K_d*K_3*K_2*B*A_s*D^2*L;             % Developed torque 
 
 %%%%%%%%%%%%%%%%% Selection of Dimensions 
  
 
 
 
 
 
 
 

