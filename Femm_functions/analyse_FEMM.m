function [ res ] = analyse_FEMM( init_geo, femm_handle, save_file )
%ANALYSE_SRM analyses the problem and provides output
% input is srm details, femm program handle and save dir
% output is torque, phase current, phase flux linkage
% directly uses save file name

% declare global FEMM handle to track FEMM process
global HandleToFEMM

% use correct handle to femm process
HandleToFEMM = femm_handle;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% analyze
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% save, analyse and load solution
mi_saveas(save_file);
mi_analyze(1);
mi_loadsolution;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% get output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get srm details
R1 = init_geo.R1;
h_r = init_geo.hr;
y_r = init_geo.yr;

g = init_geo.g;
r_middle = R1+(g/2);
theta_point = init_geo.theta_point;
x = init_geo.x;
y = init_geo.y;

%Draw a contour at the midpoint of the airgap covering all 2 poles
% mo_clearcontour() ;        % clear rany previous contours
% mo_seteditmode('contour'); 
% mo_addcontour(0,r_middle);
% mo_addcontour(0,-r_middle);
% mo_bendcontour(180,1);
% mo_addcontour(0,r_middle);
% mo_bendcontour(180,1);
% mo_makeplot(2,360);
% res.Bn = mo_lineintegral(0);
% mo_clearblock;
% Get flux density at different points in the air gap

for idx =1 :1:length(theta_point)
Bval_out(:,idx) = mo_getpointvalues(x(idx),y(idx)).'
end
res.p_x(:,1).Bx = Bval_out(2,:).' % Store as a column vector
res.p_y(:,1).By = Bval_out(3,:).' % Store as a column vector

% get torque

% select rotor
mo_selectblock( ( (R1 - h_r)+(R1 - h_r - y_r) )/2, 0);
%mo_selectblock( ( (R_By_i)+(R1 - h_r - y_r) )/2, 0);
% do block integral to find torque

res.torque = mo_blockintegral(22);
mo_clearblock;

% next phase details
% loop for all phases
for idx = 1:1:init_geo.N
    
    % find circuit name
    cir_name = sprintf('phase_%d', idx);
    
    % obtain properties
    cir_out = mo_getcircuitproperties(cir_name);
    
    % put it in results
    res.phase(idx).current = cir_out(1);
    res.phase(idx).flux_linkage = cir_out(3);
    
end

end

