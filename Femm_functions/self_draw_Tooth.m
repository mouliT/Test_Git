function [  ] = self_draw_Tooth(init_geo,femm_handle)
%% Input geometry to the femm from init_geo structure

% Number of stator poles
Ns = init_geo.Ns ;
% Number of rotor poles
Nr = init_geo.Nr ;
% Number of turns
N  = init_geo.N ;
% Stator pole angle
Thetea_sp = init_geo.Th_sp ;
% Rotor pole angle
Thetea_rp = init_geo.Th_rp ;
% Stator slot angle
%Thetea_ss = init_geo.Th_ss ;
% Rotor slot angle
%Thetea_rs = init_geo.Th_rs ;
% Variable rotar angle
Theta = init_geo.Th;
% Air gap of the machine 
g = init_geo.g ;
% Stator pole thickness
%t_s = init_geo.ts ;
% Rotor pole thickness
%t_r = init_geo.tr ;
% Rotor outer radius
R1 = init_geo.R1 ;
% stator inner radius
R2 = init_geo.R2 ;
% Rotor pole height
h_r = init_geo.hr ;
% Stack Length
L = init_geo.L_stk;
% Stator yoke width
y_s = init_geo.ys;
% Rotor yoke width
y_r = init_geo.yr;
%%%% Teeth changes added here
% Tooth depth in stator
h_teeth = init_geo.h_teeth;

Theta_teeth_p = init_geo.Theta_teeth_p ;

x_out_teeth_p = init_geo.x_out_teeth_p;
y_out_teeth_p = init_geo.y_out_teeth_p;

x_out_teeth_n = init_geo.x_out_teeth_n;
y_out_teeth_n = init_geo.y_out_teeth_n;


x_in_teeth_p  = init_geo.x_in_teeth_p;
y_in_teeth_p  = init_geo.y_in_teeth_p;

x_in_teeth_n  = init_geo.x_in_teeth_n;
y_in_teeth_n  = init_geo.y_in_teeth_n;


%% Drawing in FEMM
global HandleToFEMM

% Use correct handle to draw geomertry (Dont know why !!!)
HandleToFEMM = femm_handle;
%%
% Satrt a magnetostatic problelm
%openfemm;
newdocument(0);

%define problem
mi_probdef(0, 'meters', 'planar', 1.e-8, L, 30);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Building block of Stator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%PA- Pole arc , SC- stator Coil

% Co-Ordinates for start of stator arc(x,y)(Outer tooth circle)(x-axis centre)
Co_Stat_PA_start = [(R1+g+h_teeth), 0];

% Co-Ordinates for end of stator arc(x,y)
Co_Stat_PA_end = [(R1+g+h_teeth)*cos(Thetea_sp/2) (R1+g+h_teeth)*sin(Thetea_sp/2)];


% Co-Ordinates for start of stator coil(x,y)
Co_Stat_Co_start = [sqrt(R2^2 - (R1+g+h_teeth)^2*( sin(Thetea_sp/2) )^2 ) ...
    (R1 +g+h_teeth)*sin(Thetea_sp/2)];

% Co-Ordinates for end of stator coil(x,y)(mid poin)
Co_Stat_Co_end =  [R2*cos((pi/Ns)) R2*sin((pi/Ns))];

% Angle of stator coil 
Ang_Stat_Co_Arc = acos( dot(Co_Stat_Co_start,...
    Co_Stat_Co_end)/R2^2 );

% Coil Gap mid point between two stator poles
coil_gap = [(R1+g+h_teeth)*cos(pi/Ns) (R1+g+h_teeth)*sin(pi/Ns)];

% Coil arc angle
coil_angle = (pi/Ns - Thetea_sp/2);

% % draw half of stator pole arc centered in x-axis
% mi_drawarc(Co_Stat_PA_start(1), Co_Stat_PA_start(2), ...
%     Co_Stat_PA_end(1), Co_Stat_PA_end(2), ...
%     (Thetea_sp/2)*180/pi, 2.5);
% Add node for the pont on outer stator teeth radius 
mi_addnode(x_out_teeth_p(3), y_out_teeth_p(3))

% Draw arc for first tooth of  stator inner portion
mi_drawarc(x_in_teeth_p(2), y_in_teeth_p(2), ...
    x_in_teeth_p(3), y_in_teeth_p(3), ...
    (Thetea_rp)*180/pi, 2.5);


% Draw arc for first slot of stator outer portion
mi_drawarc(x_out_teeth_p(1), y_out_teeth_p(1), ...
    x_out_teeth_p(2), y_out_teeth_p(2), ...
    ((2*pi/Nr)-Thetea_rp)*(180/pi), 2.5);

% Draw arc half of second tooth inner portion centered with x-axis
mi_drawarc(R1+g ,0, ...
    x_in_teeth_p(1), y_in_teeth_p(1), ...
    (Thetea_rp/2)*180/pi, 2.5);

% draw stator coil arc
mi_drawarc(Co_Stat_Co_start(1), Co_Stat_Co_start(2), ...
    Co_Stat_Co_end(1), Co_Stat_Co_end(2), ...
    Ang_Stat_Co_Arc*180/pi, 2.5);


% join the first tooth end edge to pole 
mi_addsegment(x_in_teeth_p(2), y_in_teeth_p(2), ...
    x_out_teeth_p(2), y_out_teeth_p(2));

% join the outer edge of pole to the inner teeth edge
mi_addsegment(x_in_teeth_p(3), y_in_teeth_p(3), ...
    x_out_teeth_p(3), y_out_teeth_p(3));

% join the first slot end point to second tooth start point 
mi_addsegment(x_in_teeth_p(1), y_in_teeth_p(1), ...
    x_out_teeth_p(1), y_out_teeth_p(1));


% draw line from end of stator pole arc to start of stator coil arc
mi_addsegment(Co_Stat_PA_end(1), Co_Stat_PA_end(2), ...
    Co_Stat_Co_start(1), Co_Stat_Co_start(2));

% add node at stator coil point midway between two stator poles
mi_addnode(coil_gap(1), coil_gap(2));

% join the nodes to form stator coil unit
mi_addsegment(coil_gap(1), coil_gap(2), ...
    Co_Stat_Co_end(1), Co_Stat_Co_end(2));
mi_addarc(Co_Stat_PA_end(1), Co_Stat_PA_end(2), ...
    coil_gap(1), coil_gap(2), ...
    coil_angle*180/pi, 2.5);





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% copy the building block to form the entire stator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% select the basic unit arcs
mi_selectarcsegment(x_out_teeth_p(2), y_out_teeth_p(2));
mi_selectarcsegment(x_in_teeth_p(2), y_in_teeth_p(2));
mi_selectarcsegment(coil_gap(1), coil_gap(2));
mi_selectarcsegment(Co_Stat_Co_start(1), Co_Stat_Co_start(2));
mi_selectarcsegment(R1+g ,0);

% select the basic unit segments
mi_selectsegment(Co_Stat_PA_end(1), Co_Stat_PA_end(2));
mi_selectsegment(x_in_teeth_p(2), y_in_teeth_p(2));
mi_selectsegment(x_in_teeth_p(1), y_in_teeth_p(1));

mi_selectsegment(Co_Stat_Co_start(1), Co_Stat_Co_start(2));
mi_selectsegment(coil_gap(1), coil_gap(2));

% set the group no of basic unit to 1
mi_setgroup(1);

% mirror the objects about x axis
mi_selectgroup(1);
mi_mirror(0, 0, R1, 0);

% copy around axis to create stator
mi_selectgroup(1);
mi_copyrotate(0, 0, (2*pi/Ns)*180/pi, Ns);
mi_clearselected;

% draw outer circle to complete the stator
mi_addnode( R2 + y_s, 0);
mi_addnode( -(R2 + y_s), 0);
mi_addarc(R2 + y_s, 0, -(R2 + y_s), 0, 180, 2.5);
mi_addarc( -(R2 + y_s), 0, (R2 + y_s), 0, 180, 2.5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% draw rotor building block
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% point to start rotor pole arc
Co_Rot_PA_start = [R1 0];

% end of rotor pole arc
Co_Rot_PA_end = [R1*cos(Thetea_rp/2) R1*sin(Thetea_rp/2)];

% start of rotor yoke arc
Co_Rot_YA_start = [(R1*cos(Thetea_rp/2) - h_r) (R1)*sin(Thetea_rp/2)];

% rotor yoke arc end
rotor_yoke_rad = sqrt( R1^2 - 2*R1*h_r*cos(Thetea_rp/2) + h_r^2 );
Co_Rot_YA_end = [rotor_yoke_rad*cos(pi/Nr) rotor_yoke_rad*sin(pi/Nr)];

% angle of rotor yoke arc
Ang_Rot_YA = acos( dot(Co_Rot_YA_start,...
    Co_Rot_YA_end)/rotor_yoke_rad^2 );

%%%%%%%%%%%% Draw 
%draw half of rotor pole arc centered in x-axis
mi_drawarc(Co_Rot_PA_start(1), Co_Rot_PA_start(2), ...
   Co_Rot_PA_end(1), Co_Rot_PA_end(2), ...
    (Thetea_rp/2)*180/pi, 2.5);

%draw rotot yoke arc
mi_drawarc(Co_Rot_YA_start(1), Co_Rot_YA_start(2), ...
    Co_Rot_YA_end(1), Co_Rot_YA_end(2), ...
    Ang_Rot_YA*180/pi, 2.5);

%draw line from end of rotor pole arc to start of rotor yoke arc
mi_addsegment(Co_Rot_PA_end(1), Co_Rot_PA_end(2), ...
    Co_Rot_YA_start(1), Co_Rot_YA_start(2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% copy rotor basic unit to make the whole rotor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% select the basic unit arcs
mi_selectarcsegment(Co_Rot_PA_start(1), Co_Rot_PA_start(2));
mi_selectarcsegment(Co_Rot_YA_start(1), Co_Rot_YA_start(2));

% select the basic unit segments
mi_selectsegment(Co_Rot_PA_end(1), Co_Rot_PA_end(2));

% set the group no of basic unit to 1
mi_setgroup(2);

% mirror the objects about x axis
mi_selectgroup(2);
mi_mirror(0, 0, R1, 0);

% copy around axis to create rotor
mi_selectgroup(2);
mi_copyrotate(0, 0, (2*pi/Nr)*180/pi, Nr);
mi_clearselected;

% draw inner circle to complete the rotor
mi_addnode( R1 - y_r - h_r, 0);
mi_addnode( -(R1 - y_r - h_r), 0);
mi_addarc(R1 - y_r - h_r, 0, -(R1 - y_r - h_r), 0, 180, 2.5);
mi_addarc( -(R1 - y_r - h_r), 0, (R1 - y_r - h_r), 0, 180, 2.5);

% rotate the rotor to theta
mi_selectgroup(2);
mi_moverotate(0, 0, Theta*180/pi);

end