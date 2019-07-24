function [ result ] = Pass_Par_femm(init_geo, mmf, save_file)
%SRM_FEMM_PAR Test to make femm work in parallel

% global FEMM handle to track FEMM process
global HandleToFEMM
% open femm and store the process handle
openfemm;

femm_h = HandleToFEMM;

% draw the geometry. select correct draw function
self_draw_Tooth(init_geo, femm_h);

% set-up materials, boundary and excitation
setup_analysis_Tooth(init_geo, mmf, femm_h);

% analyze the problem
result = analyse_FEMM(init_geo, femm_h, save_file);

% finally close femm
closefemm;

end

