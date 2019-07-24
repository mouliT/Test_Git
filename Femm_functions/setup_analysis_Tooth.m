function [ ] = setup_analysis_Tooth( init_geo, mmf, femm_handle)

%SETUP_ANALYSIS setup the analysis of srm
% setup materials, boundary and excitation

% declare global FEMM handle to track FEMM process
global HandleToFEMM

% use correct handle to femm process
HandleToFEMM = femm_handle;

%% get few relevant dimensions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R1 = init_geo.R1;
h_r = init_geo.hr;
y_r = init_geo.yr;
theta = init_geo.Th;
Nr = init_geo.Nr;
R2 = init_geo.R2;
y_s = init_geo.ys;
Ns = init_geo.Ns;
Thetea_sp = init_geo.Th_sp;
N = init_geo.N;
g = init_geo.g;
% %%%% barrier parameters
% h_u = init_geo.hu;
% R11 = init_geo.R11;
% t_m = init_geo.tm;
% t_s = init_geo.ts ;
% t_b = init_geo.tb ;
% beta_ri = init_geo.beta_ri;
% beta_ro = init_geo.beta_ro ;
% h_r1 = init_geo.h_r1;
% h_r2 = init_geo.h_r2 ;
% % R_By_o = init_geo.R_By_o;
% % R_By_i = init_geo.R_By_i;

%% set-up and assign materials to stator, rotor, air-gap and shaft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get material iron m-19
mi_getmaterial('M-19 Steel');
%mi_addmaterial('M-19 Steel', 4416, 4416, 0, 0, 0, 0, 0, 1, 0, 0, 0);
% air material
mi_getmaterial('Air');

% copper material for coil
mi_getmaterial('Copper');

% set shaft to air
mi_addblocklabel(0, 0);
mi_selectlabel(0,0);
mi_setblockprop('Air', 1, 0, '<None>', 0, 0, 0);
mi_clearselected;

% set rotor to M-19 steel
mi_addblocklabel( ( (R1 - h_r)+(R1 - h_r - y_r) )/2, 0 );
mi_selectlabel( ( (R1 - h_r)+(R1 - h_r - y_r) )/2, 0 );
% mi_addblocklabel( ( (R_By_i)+(R1 - h_r - y_r) )/2, 0 );
% mi_selectlabel( ( (R_By_i)+(R1 - h_r - y_r) )/2, 0 );
mi_setblockprop('M-19 Steel', 1, 0, '<None>', 0, 0, 0);
mi_clearselected;

% set air-gap to air
% first select air-gap independent of rotor position
% air_x = (R1 - h_r/2)*cos(theta + pi/Nr);
% air_y = (R1 - h_r/2)*sin(theta + pi/Nr);

air_x = (R1 + g/2)*cos(theta + pi/Nr);
air_y = (R1 + g/2)*sin(theta + pi/Nr);
% set block properties
mi_addblocklabel(air_x, air_y);
mi_selectlabel(air_x, air_y);
mi_setblockprop('Air', 1, 0, '<None>', 0, 0, 0);
mi_clearselected;
% set stator to M-19 Steel
mi_addblocklabel(R2 + y_s/2, 0);
mi_selectlabel(R2 + y_s/2, 0);
mi_setblockprop('M-19 Steel', 1, 0, '<None>', 0, 0, 0);
mi_clearselected;

%% set-up windings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% first set-up circuit properties
for idx = 1:1:init_geo.N
    
    str = sprintf('phase_%d', idx);
    mi_addcircprop(str, mmf(idx), 1);
    
end

% first x and y cooridinate of middle of stator coil basic unit.
% start from here
x = (R2 + R1)*cos( (pi/Ns + Thetea_sp/2)/2 )/2;
y = (R2 + R1)*sin( (pi/Ns + Thetea_sp/2)/2 )/2;

% angle between end of coil and block label
alpha = pi/(2*Ns) - Thetea_sp/4;

% start assigning the coils following logic in copy phd-02
% count circuit number and coil number as well
% going clockwise
cir_num = 0;
turns = 80;
for coil_num = 1:1:2*Ns
    
    % if coil number is odd change circuit number
    if ( mod(coil_num,2) == 1)
        cir_num = cir_num + 1;
        % if circuit number exceeds phase number then reset
        % and if phases is even then reset turns as well
        if(cir_num == N + 1)
            cir_num = 1;
            if( mod(N,2) == 0)
                turns = - turns;
            end
        end
    end
    
    % if coil number is even change turns number
    if ( mod(coil_num,2) == 0)
        turns = -turns;
    end
    
    % add material to this coil
    mi_addblocklabel(x, y);
    mi_selectlabel(x, y);
    cir_name = sprintf('phase_%d', cir_num);
    mi_setblockprop('Copper', 1, 0, cir_name, 0, 0, turns);
    mi_clearselected;
    
    % go to next coil
    
    % if present coil is odd then there is a pole section next
    if (mod(coil_num,2) == 1)
         
        % rotation angle
        gamma = - (2*alpha + Thetea_sp);
        
        % new coordinates
        x1 = x*cos(gamma) - y*sin(gamma);
        y1 = x*sin(gamma) + y*cos(gamma);
        x = x1;
        y = y1;
       
    end
    
    % if present coil is even then next coil is directly ahead
    if (mod(coil_num,2) == 0)
        
        % rotation angle
        gamma = - (2*alpha);
        
        % new coordinates
        x1 = x*cos(gamma) - y*sin(gamma);
        y1 = x*sin(gamma) + y*cos(gamma);
        x = x1;
        y = y1;
       
    end
end
%% set-up boundary
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% new boundary property with A=0
mi_addboundprop('bnd', 0, 0, 0, 0, 0, 0, 0, 0, 0)

% set upper half boundary
mi_selectarcsegment(0, R2+y_s);
mi_setarcsegmentprop(2.5, 'bnd', 0, 0);
mi_clearselected;

% set lower half boundary
mi_selectarcsegment(0, -R2-y_s);
mi_setarcsegmentprop(2.5, 'bnd', 0, 0);
mi_clearselected;
%mi_makeABC();

