%% obtain complete flux linkage and torque characteristics using FEMM
set(groot,'defaultfigureposition',[100 100 1000 550])
% add path to femm functions
addpath( [pwd '\Femm_functions']);
%addpath( [pwd '\Femm_store']);
%addpath( [pwd '\mfiles']);
%% obtain torque characteristics and Si-i char
% m_idx = 0;
% j =0;
% %for m_idx = m_idx:1:20
%     j = j+1;
% load nominal geometry of the machine
load('init_geo.mat');
% set correct angle
init_geo.Th = 0*pi/180;
% the mmf vector
mmf_femm_vec = linspace(0,70,10);
%mmf_femm_vec = i_m2_t(j);
% set-up the multiple geometries -> vary angle. note the negative starting
% point
theta_femm = (-9)*(pi/180):1*(pi/180):0*(pi/180);%linspace((-180/geom.Nr)*pi/180,0,25);
% Initialize the same theta_femm in analyse_FEMM 
%theta_femm = Th_m2_t(j)*(pi/180);
for idx = 1:1:length(theta_femm)
    
    geom_vec(idx) = init_geo();
    geom_vec(idx).Th = theta_femm(idx);
    
end

% analyze all the variants

% each row is for one theta
res_flux = zeros(length(theta_femm),length(mmf_femm_vec)); % flux linkage
res_torque = zeros(length(theta_femm),length(mmf_femm_vec)); % torque

% outer loop is for angle
% inner parallel loop for mmf

for idx_theta = 1:length(theta_femm)

    % the geometry for femm
    geom_femm = geom_vec(idx_theta);
    
    % temp array to store the torque and flux linkage
    temp_flux = zeros(length(mmf_femm_vec),1);
    temp_torque = zeros(length(mmf_femm_vec),1);
    
    parfor idx_mmf = 1:length(mmf_femm_vec)
    
        % generate save file name
    save_file = sprintf('Femm_functions/femm_temp/femm_%d_%d.fem', idx_theta,...
        idx_mmf);

    % call femm
    result = Pass_Par_femm(geom_femm, [mmf_femm_vec(idx_mmf) 0 0],...
        save_file);


    % get torque
    temp_torque(idx_mmf) = result.torque;
    % get flux linkage
    temp_flux(idx_mmf) = result.phase(1).flux_linkage;
    % get current 
    temp_curr(idx_mmf) = result.phase(1).current;
    %calc inductance
    temp_L(idx_mmf) =  temp_flux(idx_mmf)/temp_curr(idx_mmf);
    % get normal component of torque
    temp_Bx(:,idx_mmf) = result.p_x.Bx; % Store as a  column vector for fixed mmf
    %disp(temp_Bx(:,idx_mmf));
    %disp('idx_mmf');
    temp_By(:,idx_mmf) = result.p_y.By; % Store as a column vector for fixed mmf
    end
    
    % store torque in matrix
    res_torque(idx_theta,:) = temp_torque;
    %res_torque_mat(j,:) =  temp_torque;
    % store flux linkage in matrix
    res_flux(idx_theta,:) = temp_flux;
    %res_flux_mat(j,:) = temp_flux;
    % store current in a matrix
    res_curr(idx_theta,:) = temp_curr;
    % Compute inductance
    res_L(idx_theta,:) = temp_L;
    %res_L_mat(j,:) =  temp_L;
    res_Bx(:,:,idx_theta) = temp_Bx;
    res_By(:,:,idx_theta) = temp_By;
    
end
%end

%% Plotting

% Si _ i _ theta

for p =1:1:length(res_flux)
  figure(1)
  %if p == 1||p ==11
  plot(res_curr(1,:),res_flux(p,:)*1e3,'-','markersize',2,'linewidth',3)
  %end
  hold on
  grid on
  xlabel('Phase Current(A)','Fontsize',18,'Fontname','Times');
  ylabel('Flux linkages(mWb)','Fontsize',18,'Fontname','Times');
  title('Flux linkage vs. Phase current','Fontsize',24,'Fontname','Times');
  set(gcf,'color','white');
  set(gca,'Fontsize',24);
  set(gca,'linewidth',3,'Fontsize',24);
  

%   %legend({'y = 0cm','y = 0.1cm','y = 0.2cm','y = 0.3cm','y = 0.4cm',...
%                 'y = 0.5cm','y = 0.6cm','y = 0.7cm','y = 0.75cm','y = 79cm','y = 0.799cm'},'Location','northwest')
end
%% Plotting

%T _ i _ theta

for p =1:1:length(res_torque)
    %colormap(parula)
    %cmap=colormap;
    %Plot_color=cmap((d),:);
  figure(2)
  %if  p == 1 || p==11
  plot((theta_femm)*(180/pi),res_torque(:,p),'-','markersize',2,'linewidth',3)
  hold on
  grid on
  %end
  xlabel('Theta,(Degrees)','Fontsize',18,'Fontname','Times');
  ylabel('Torque(N-m)','Fontsize',18,'Fontname','Times');
  title('Torque vs. Theta','Fontsize',24,'Fontname','Times');
  set(gcf,'color','white');
  set(gca,'Fontsize',24);
  set(gca,'linewidth',3,'Fontsize',24);

  %legend({'1A','70A'},'Location','northwest')
                 
end
%% Matlab inductance
figure(2)
 subplot(2,2,1);
 plot(Th_m1Copy,res_L_mat_M1Copy*1e3,':','markersize',2,'linewidth',3)
 axis([-45 0 0 6])
% xlabel('Theta,(Degrees)','Fontsize',18,'Fontname','Times');
 ylabel('Inductance(mH)','Fontsize',18,'Fontname','Times');
 title('Inductance vs. Theta','Fontsize',24,'Fontname','Times');
 set(gcf,'color','white');
 set(gca,'Fontsize',24);
 set(gca,'linewidth',3,'Fontsize',24);
 hold on 
 grid on
 subplot(2,2,3);
 plot(Th_m1Copy,i_m1Copy,':','markersize',2,'linewidth',3)
 axis([-45 0 0 15])
 hold on 
 grid on
 xlabel('Theta,(Degrees)','Fontsize',18,'Fontname','Times');
 ylabel('Current(A)','Fontsize',18,'Fontname','Times');
 title('Current vs. Theta','Fontsize',24,'Fontname','Times');
 set(gcf,'color','white');
 set(gca,'Fontsize',24);
 set(gca,'linewidth',3,'Fontsize',24);

 subplot(2,2,2);
 
  plot(Th_m1Copy,res_torque_mat_M1Copy,':','markersize',2,'linewidth',3)
 axis([-45 0 0 0.75])
 hold on 
 grid on
 
 ylabel('Torque(N-m)','Fontsize',18,'Fontname','Times');
 title('Torque vs. Theta','Fontsize',24,'Fontname','Times');
 set(gcf,'color','white');
 set(gca,'Fontsize',24);
 set(gca,'linewidth',3,'Fontsize',24);

  subplot(2,2,4);
 
 plot(Th_m2,V_pulse,':','markersize',2,'linewidth',3)
 axis([-45 0 -100 100])
 hold on 
 grid on
 xlabel('Theta,(Degrees)','Fontsize',18,'Fontname','Times');
 ylabel('Voltage(V)','Fontsize',18,'Fontname','Times');
 title('Voltage vs. Theta','Fontsize',24,'Fontname','Times');
 set(gcf,'color','white');
 set(gca,'Fontsize',24);
 set(gca,'linewidth',3,'Fontsize',24);
 
%% Inductance

for p =1:1:length(res_L)
    
  figure(3)  
  if p == 11 
  plot(theta_femm*(180/pi),res_L(:,p)*1e3,'-','markersize',2,'linewidth',3)
  hold on
  grid on
  end
  xlabel('Angle(degrees)','Fontsize',18,'Fontname','Times');
  ylabel('Inductnce(mH)','Fontsize',18,'Fontname','Times');
  title('Inductance  vs.  rotor angle','Fontsize',24,'Fontname','Times');
  set(gcf,'color','white');
  set(gca,'Fontsize',24);
  set(gca,'linewidth',3,'Fontsize',24);
  %legend({'10A','20A','70A'},'Location','northwest')
end
%%
 plot(Th_m1Copy,res_L_mat_M1Copy*1e3,'-+','markersize',6,'linewidth',3)
%% slope of inductance
r=1;
dL_dTh = zeros(length(res_L),length(res_L));
for d = 1:1:length(res_L)
    colormap(parula)
    cmap=colormap;
    Plot_color=cmap((d),:);
    r = r+1;
    L_f = res_L(:,r); % Loaded sI-I plot at a fixed position 
    N = length(L_f);
    dL_f = diff(L_f);
    dTh = diff(theta_femm);
    dL_f_dTh = (dL_f)./(dTh.'); % Keeping position fixed
    dL_f_dTh(N) = dL_f_dTh(N-1);
    dL_dTh(r,:) = dL_f_dTh;  % Again storing the result as a row vector ...
                                % for a consatnt position 
    figure(1)
    subplot(1,2,1);
    plot((theta_femm*(180/pi)),L_f,'--','markersize',5,'linewidth',1.5,'Color', Plot_color);
    xlabel('Theta,radians','Fontsize',18,'Fontname','Times');
    ylabel('Inductance,Henry','Fontsize',18,'Fontname','Times');
    title('Si vs. current','Fontsize',24,'Fontname','Times');
    set(gcf,'color','white');
    set(gca,'Fontsize',24);
    set(gca,'linewidth',3,'Fontsize',24);
    %legend({'x = 2cm','x = 1.75cm','x = 1.5cm','x = 1.25','x = 1cm',...
                %'x = 0.75cm','x = 0.5cm','x = 0.25cm','x = 0cm'},'Location','northwest')
    hold on
    grid on
    subplot(1,2,2);
    plot((theta_femm*(180/pi)),dL_dTh(r,:),'--','markersize',5,'linewidth',1.5,'Color', Plot_color);
    xlabel('Theta(radians)','Fontsize',18,'Fontname','Times');
    ylabel('Slope,Henry/rad','Fontsize',18,'Fontname','Times');
    title('slope of inductance vs. Theta','Fontsize',24,'Fontname','Times');
    set(gcf,'color','white');
    set(gca,'Fontsize',24);
    set(gca,'linewidth',3,'Fontsize',24);
    %legend({'x = 2cm','x = 1.75cm','x = 1.5cm','x = 1.25','x = 1cm',...
                %'x = 0.75cm','x = 0.5cm','x = 0.25cm','x = 0cm'},'Location','northeast')
    hold on
    grid on
end

%% Differnetiation of Si-I graph
% In res_flux(Th,i)......Row fix means position fix  // column fix means current fix
r=0;
dsi_di = zeros(length(res_flux),length(res_flux));
for d = 1:1:length(res_flux)
    colormap(parula)
    cmap=colormap;
    Plot_color=cmap((d),:);
    r = r+1;
    si_f = res_flux(r,:); % Loaded sI-I plot at a fixed position 
    N = length(si_f);
    dsi_f = diff(si_f);
    dcurr = diff(mmf_femm_vec);
    dsi_f_dcurr = (dsi_f)./(dcurr); % Keeping position fixed
    dsi_f_dcurr(N) = dsi_f_dcurr(N-1);
    dsi_di(r,:) = dsi_f_dcurr;  % Again storing the result as a row vector ...
                                % for a consatnt position 
    figure(1)
    subplot(1,2,1);
    plot(mmf_femm_vec,si_f,'-','markersize',5,'linewidth',1.5,'Color', Plot_color);
    xlabel('Current,Amperes','Fontsize',18,'Fontname','Times');
    ylabel('Flux linkages,weber turns','Fontsize',18,'Fontname','Times');
    title('Si vs. current','Fontsize',24,'Fontname','Times');
    set(gcf,'color','white');
    set(gca,'Fontsize',24);
    set(gca,'linewidth',3,'Fontsize',24);
    %legend({'x = 2cm','x = 1.75cm','x = 1.5cm','x = 1.25','x = 1cm',...
                %'x = 0.75cm','x = 0.5cm','x = 0.25cm','x = 0cm'},'Location','northwest')
    hold on
    grid on
    subplot(1,2,2);
    plot(mmf_femm_vec,dsi_di(r,:),'-','markersize',5,'linewidth',1.5,'Color', Plot_color);
    xlabel('Current,Amperes','Fontsize',18,'Fontname','Times');
    ylabel('Slope,Henry','Fontsize',18,'Fontname','Times');
    title('Incremental Inductance vs. current','Fontsize',24,'Fontname','Times');
    set(gcf,'color','white');
    set(gca,'Fontsize',24);
    set(gca,'linewidth',3,'Fontsize',24);
    %legend({'x = 2cm','x = 1.75cm','x = 1.5cm','x = 1.25','x = 1cm',...
                %'x = 0.75cm','x = 0.5cm','x = 0.25cm','x = 0cm'},'Location','northeast')
    hold on
    grid on
end

%% Differnetiation of Si-Theta graph
r=0;
dsi_dTh = zeros(length(res_flux),length(res_flux));
for d = 1:1:length(res_flux)
    colormap(parula)
    cmap=colormap;
    Plot_color=cmap((d),:);
    r = r+1;
    si_Th = res_flux(:,r); % Fixing the coloumns implies fixing the curren...
                           % and varying the poisitions 
    N = length(si_Th);
    dsi_f = diff(si_Th);
    dTh = diff(theta_femm);
    dsi_f_dTh = (dsi_f)./(dTh.');
    dsi_f_dTh(N) = dsi_f_dTh(N-1);
    dsi_dTh(:,r) = dsi_f_dTh; % Again storing the result in a column....
                              % rows being position is preserved
    figure(8)
    subplot(1,2,1);
    plot(theta_femm,si_Th,'-','markersize',5,'linewidth',1.5,'Color', Plot_color);
    xlabel('Theta(radians)','Fontsize',18,'Fontname','Times');
    ylabel('Flux linkages,weber turns','Fontsize',18,'Fontname','Times');
    title('Si vs. Position','Fontsize',24,'Fontname','Times');
    set(gcf,'color','white');
    set(gca,'Fontsize',24);
    set(gca,'linewidth',3,'Fontsize',24);
    %legend({'x = 2cm','x = 1.75cm','x = 1.5cm','x = 1.25','x = 1cm',...
     %           'x = 0.75cm','x = 0.5cm','x = 0.25cm','x = 0cm'},'Location','northwest')
    hold on
    grid on
    subplot(1,2,2);
    plot(theta_femm,dsi_dTh(:,r),'-','markersize',5,'linewidth',1.5,'Color', Plot_color);
    xlabel('Theta(Radians)','Fontsize',18,'Fontname','Times');
    ylabel('dsi by dTheta','Fontsize',18,'Fontname','Times');
    title('dSi by dTheta vs. Theta','Fontsize',24,'Fontname','Times');
    set(gcf,'color','white');
    set(gca,'Fontsize',24);
    set(gca,'linewidth',3,'Fontsize',24);
    %legend({'x = 2cm','x = 1.75cm','x = 1.5cm','x = 1.25','x = 1cm',...
     %           'x = 0.75cm','x = 0.5cm','x = 0.25cm','x = 0cm'},'Location','northeast')
    hold on
    grid on
end

%% Lookup data 10 by 10
% dsi/di
dsi_di_LT = zeros(length(res_flux),length(res_flux));
dsi_di_LT(:,:) = [dsi_di(1,:);dsi_di(2,:);dsi_di(3,:);dsi_di(4,:);dsi_di(5,:);...
                        dsi_di(6,:);dsi_di(7,:);dsi_di(8,:);dsi_di(9,:);dsi_di(10,:)];
dsi_di_LT_T = dsi_di_LT(:,:).';
% dsi/dTh
dsi_dTh_LT = zeros(length(res_flux),length(res_flux));
dsi_dTh_LT(:,:) = [dsi_dTh(:,1) dsi_dTh(:,2) dsi_dTh(:,3) dsi_dTh(:,4) dsi_dTh(:,5) ...
                        dsi_dTh(:,6) dsi_dTh(:,7) dsi_dTh(:,8) dsi_dTh(:,9) dsi_dTh(:,10)];
Si_i_Th = zeros(length(res_flux),length(res_flux));
Si_i_Th_LT(:,:) =  [res_flux(1,:);res_flux(2,:);res_flux(3,:);res_flux(4,:);res_flux(5,:);...
                        res_flux(6,:);res_flux(7,:);res_flux(8,:);res_flux(9,:);res_flux(10,:)];
                    
Si_i_Th_LT_T = Si_i_Th_LT(:,:).';
                    
T_i_Th = [res_torque(1,:);res_torque(2,:);res_torque(3,:);res_torque(4,:);res_torque(5,:);...
                        res_torque(6,:);res_torque(7,:);res_torque(8,:);res_torque(9,:);res_torque(10,:)];
                    
bp2d_i =mmf_femm_vec;
bp2d_Th =theta_femm;
%% Look up 20 by 20

dsi_di_LT = zeros(length(res_flux),length(res_flux));
dsi_di_LT(:,:) = [dsi_di(1,:);dsi_di(2,:);dsi_di(3,:);dsi_di(4,:);dsi_di(5,:);...
                        dsi_di(6,:);dsi_di(7,:);dsi_di(8,:);dsi_di(9,:);dsi_di(10,:);...
                        dsi_di(11,:);dsi_di(12,:);dsi_di(13,:);dsi_di(14,:);dsi_di(15,:);...
                        dsi_di(16,:);dsi_di(17,:);dsi_di(18,:);dsi_di(19,:);dsi_di(20,:)];
                  
dsi_di_LT_T = dsi_di_LT(:,:).';
% dsi/dTh

dsi_dTh_LT = zeros(length(res_flux),length(res_flux));
dsi_dTh_LT(:,:) = [dsi_dTh(:,1) dsi_dTh(:,2) dsi_dTh(:,3) dsi_dTh(:,4) dsi_dTh(:,5) ...
                        dsi_dTh(:,6) dsi_dTh(:,7) dsi_dTh(:,8) dsi_dTh(:,9) dsi_dTh(:,10)...
                        dsi_dTh(:,11) dsi_dTh(:,12) dsi_dTh(:,13) dsi_dTh(:,14) dsi_dTh(:,15)...
                        dsi_dTh(:,16) dsi_dTh(:,17) dsi_dTh(:,18) dsi_dTh(:,19) dsi_dTh(:,20)];
                  
Si_i_Th = zeros(length(res_flux),length(res_flux));
Si_i_Th_LT(:,:) =  [res_flux(1,:);res_flux(2,:);res_flux(3,:);res_flux(4,:);res_flux(5,:);...
                        res_flux(6,:);res_flux(7,:);res_flux(8,:);res_flux(9,:);res_flux(10,:);...
                        res_flux(11,:);res_flux(12,:);res_flux(13,:);res_flux(14,:);res_flux(15,:);...
                        res_flux(16,:);res_flux(17,:);res_flux(18,:);res_flux(19,:);res_flux(20,:)];
                    
Si_i_Th_LT_T = Si_i_Th_LT(:,:).';

                    
T_i_Th = [res_torque(1,:);res_torque(2,:);res_torque(3,:);res_torque(4,:);res_torque(5,:);...
                        res_torque(6,:);res_torque(7,:);res_torque(8,:);res_torque(9,:);res_torque(10,:);...
                        res_torque(11,:);res_torque(12,:);res_torque(13,:);res_torque(14,:);res_torque(15,:);...
                        res_torque(16,:);res_torque(17,:);res_torque(18,:);res_torque(19,:);res_torque(20,:)];
                    
bp2d_i =mmf_femm_vec;
bp2d_Th =theta_femm;