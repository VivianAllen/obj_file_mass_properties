function [Smass, SCoM, Stensor_origin, Stensor_com] = ...
    getMassPropSystem(varargin)
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   AUTHOR Vivian Allen (mrvivianallen@gmail.com)                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% VERSION: 2.3 (06 March 2017)
% fixed negative CoM bug
 

% VERSION: 2.2 (14 Sept 2016)
% version changes (2.2) PLANNED - add option to make density a vector of
% equal length to positive body list. Allow variable body density.
% version changes (2.1) - switched output order of origin and com tensors 
% to match older code 
% version changes (2.0) - general cleanup, integration with new version of
% get_body_mass_props.mmasses=

% loops through a list of body objects (.obj files) defining a mass system 
% and estimates the compound mass properties of the system by calling 
% get_body_mass_props.m and summing output. Mass properties for  cavities 
% (negative mass objects), if listed in a seperate cell array of strings 
% (see below) are calculated seperately and then subtracted from the 
% system.

%%%%%%%%%%%%%%%%%%%% INPUT EXPLANATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 	pos_wobjs = cell array of strings listing target names (.obj) that    %
%   are positive masses                                                   %
% 	neg_wobjs = (OPTIONAL) cell array of strings listing target names     %
%   (.obj) that are negative masses (e.g. cavities)                       %
% 	density = density of body in kgs per m^3 (1060 is standard)           % 
%   fidelity = control for size of ray array (array will be fidility^2    %
% 	rays). Higher is more accurate, exponentially slower. 50 seems good.  %
% 	units = 1 for metres, 0.01 for centimetres, etc. MESH WILL BE SCALED  %
%   TO METERS BASED ON UNITS                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%% OUTPUT EXPLANATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 	Smass = estimated mass of target system (kg)                          %
% 	SCoM = estimated CoM of target syste, (m)                             % 
% 	Stensor_com = inertial tensor about CoM of target system              %
% 	Stensor_origin = inertial tensor about origin of target system        %
% (NB use ~ to suppress output of any variable)                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% REQUIRED FUNCTIONS (must be in same folder):

% get_body_mass_props (AUTHOR: Vivian Allen)

%% ------------------------ Detect Input mode ----------------------------%
% CASE 1: positive AND negative bodies in system (segments and cavities).
% INPUT EXAMPLE: get_system_mass_props(pos_wobjs, neg_wobjs, fidelity, units, zeroName)

if nargin == 5
    runmode = 1;
    pos_wobjs = varargin{1};
    neg_wobjs = varargin{2};
    density = varargin{3};
    fidelity = varargin{4};
    units = varargin{5};
    
% CASE 2: positive ONLY bodies in system (segments).
% INPUT EXAMPLE: get_system_mass_props(pos_wobjs, fidelity, units, zeroName)  
elseif nargin == 4
    runmode = 2;
    pos_wobjs = varargin{1};
    density = varargin{2};
    fidelity = varargin{3};
    units = varargin{4};

else
    error('WRONG NUMBER OF INPUTS!')
end

%% ---------- get and collate mass props for pos and neg wobjlists -------%

% 1: positive segments

% intialise segment mass storage
segment_masses = 0;

% intialise segment CoM storage
segment_CoMs = [0 0 0];

% intialise segment MOI storage
segment_MOIs = [0 0 0];

% intialise segment POI storage
segment_POIs = [0 0 0];

for loop1 = 1:length(pos_wobjs);
    
    % detect if variable density or not
    if size(density,1)>1
        curr_density = density(loop1);
    else
        curr_density = density;
    end
    
    % get data for current body
    [seg_mass, seg_CoM, seg_tensorO, ~] = ...
        getMassProps( cell2mat(pos_wobjs(loop1)), curr_density, fidelity, units);
    
    % append mass
    segment_masses = [segment_masses; seg_mass];
    
    % append CoM
    segment_CoMs = [segment_CoMs; seg_CoM];
    
    % deconstruct tensor into MOI[x y z] & POI[XY XZ YZ]
    seg_MOI = [seg_tensorO(1,1) seg_tensorO(2,2) seg_tensorO(3,3)];
    seg_POI = [seg_tensorO(2,1) seg_tensorO(3,1) seg_tensorO(3,2)];
    
    % append MOI
    segment_MOIs = [segment_MOIs; seg_MOI];
    
    % append POI
    segment_POIs = [segment_POIs; seg_POI];
    
end

% trim!
segment_masses(1) = [];
segment_CoMs(1,:) = [];
segment_MOIs(1,:) = [];
segment_POIs(1,:) = [];

% 2: negative segments (if present

if runmode==1

    % intialise cavity mass storage
    cavity_masses = 0;

    % intialise cavity CoM storage
    cavity_CoMs = [0 0 0];

    % intialise cavity MOI storage
    cavity_MOIs = [0 0 0];

    % intialise cavity POI storage
    cavity_POIs = [0 0 0];
    
    
    for loop1 = 1:length(neg_wobjs);
        
        % get data for current cavity
        [cav_mass, cav_CoM, cav_tensorO, ~] = ...
            getMassProps( cell2mat(neg_wobjs(loop1)), density, fidelity, units);
    
        % append mass
        cavity_masses = [cavity_masses; cav_mass];
    
        % append CoM
        cavity_CoMs = [cavity_CoMs; cav_CoM];
    
        % deconstruct tensor into MOI[x y z] & POI[XY XZ YZ]
        cav_MOI = [cav_tensorO(1,1) cav_tensorO(2,2) cav_tensorO(3,3)];
        cav_POI = [cav_tensorO(2,1) cav_tensorO(3,1) cav_tensorO(3,2)];
    
        % append MOI
        cavity_MOIs = [cavity_MOIs; cav_MOI];
    
        % append POI
        cavity_POIs = [cavity_POIs; cav_POI];
    
    end
    
% trim!
cavity_masses(1) = [];
cavity_CoMs(1,:) = [];
cavity_MOIs(1,:) = [];
cavity_POIs(1,:) = [];

end

%% ------------------- collate pos and neg mass props --------------------%

% sum positive mass props

% sum mass
sys_pos_mass = sum(segment_masses);

% caculate pos system CoM
% EXCEPTION CATCH - If segment_CoMs has a single line,
% 'sum' will sum it to a single value rather than a single row.
if size(segment_CoMs,1)>1
sys_pos_CoM = (1/sys_pos_mass) * sum(segment_CoMs .* ...
    [segment_masses segment_masses segment_masses]);
else
sys_pos_CoM = (1/sys_pos_mass) * (segment_CoMs .* ...
    [segment_masses segment_masses segment_masses]);
end

% sum pos system MOIs (straight sum, same axis) about origin
sys_pos_X_MOI_origin = sum(segment_MOIs(:,1));
sys_pos_Y_MOI_origin = sum(segment_MOIs(:,2));
sys_pos_Z_MOI_origin = sum(segment_MOIs(:,3));

% sum pos system POIs (straight sum, same axis) about origin
sys_pos_XY_POI_origin = sum(segment_POIs(:,1));
sys_pos_XZ_POI_origin = sum(segment_POIs(:,2));
sys_pos_YZ_POI_origin = sum(segment_POIs(:,3));

% build pos system origin inertial tensor
sys_pos_origin_tensor = ...
    [sys_pos_X_MOI_origin sys_pos_XY_POI_origin sys_pos_XZ_POI_origin; ...
    sys_pos_XY_POI_origin sys_pos_Y_MOI_origin sys_pos_YZ_POI_origin; ...
    sys_pos_XZ_POI_origin sys_pos_YZ_POI_origin sys_pos_Z_MOI_origin];


if runmode==1
    
    % if cavities present:
    
    % sum negative mass props

    % sum mass
    sys_neg_mass = sum(cavity_masses);

    % caculate pos system CoM
    % EXCEPTION CATCH - If segment_CoMs has a single line,
    % 'sum' will sum it to a single value rather than a single row.
    if size(cavity_CoMs,1)>1
        sys_neg_CoM = (1/sys_neg_mass) * sum(cavity_CoMs .* ...
        [cavity_masses cavity_masses cavity_masses]);
    else
        sys_neg_CoM = (1/sys_neg_mass) * (cavity_CoMs .* ...
        [cavity_masses cavity_masses cavity_masses]);
    end

    % sum pos system MOIs (straight sum, same axis) about origin
    sys_neg_X_MOI_origin = sum(cavity_MOIs(:,1));
    sys_neg_Y_MOI_origin = sum(cavity_MOIs(:,2));
    sys_neg_Z_MOI_origin = sum(cavity_MOIs(:,3));

    % sum pos system POIs (straight sum, same axis) about origin
    sys_neg_XY_POI_origin = sum(cavity_POIs(:,1));
    sys_neg_XZ_POI_origin = sum(cavity_POIs(:,2));
    sys_neg_YZ_POI_origin = sum(cavity_POIs(:,3));

    % build pos system origin inertial tensor
    sys_neg_origin_tensor = ...
    [sys_neg_X_MOI_origin sys_neg_XY_POI_origin sys_neg_XZ_POI_origin; ...
    sys_neg_XY_POI_origin sys_neg_Y_MOI_origin sys_neg_YZ_POI_origin; ...
    sys_neg_XZ_POI_origin sys_neg_YZ_POI_origin sys_neg_Z_MOI_origin];

    % sum pos and neg to make system values
    
    % sum mass
    Smass = sys_pos_mass - sys_neg_mass;
    
    % calculate CoM
    SCoM = (1/Smass) * sum([sys_pos_CoM; sys_neg_CoM] .* ...
        [sys_pos_mass sys_pos_mass sys_pos_mass; ...
        sys_neg_mass*-1 sys_neg_mass*-1 sys_neg_mass*-1]);
    
    % sum origin tensor
    Stensor_origin = sys_pos_origin_tensor - sys_neg_origin_tensor;
    
    % calculate MOIs about CoM
    sys_X_MOI_CoM = ...
        (sys_pos_X_MOI_origin - (sum(SCoM([2 3]).^2,2).* sys_pos_mass)) - ...
        (sys_neg_X_MOI_origin - (sum(SCoM([2 3]).^2,2).* sys_neg_mass));
    sys_Y_MOI_CoM = ...
        (sys_pos_Y_MOI_origin - (sum(SCoM([1 3]).^2,2).* sys_pos_mass)) - ...
        (sys_neg_Y_MOI_origin - (sum(SCoM([1 3]).^2,2).* sys_neg_mass));
    sys_Z_MOI_CoM = ...
        (sys_pos_Z_MOI_origin - (sum(SCoM([1 2]).^2,2).* sys_pos_mass)) - ...
        (sys_neg_Z_MOI_origin - (sum(SCoM([1 2]).^2,2).* sys_neg_mass));   
    
    % calculate POIs about CoM
    sys_XY_POI_CoM = ...
        (sys_pos_XY_POI_origin + (sys_pos_mass * SCoM(1) * SCoM(2))) - ...
        (sys_neg_XY_POI_origin + (sys_neg_mass * SCoM(1) * SCoM(2)));
    sys_XZ_POI_CoM = ...
        (sys_pos_XZ_POI_origin + (sys_pos_mass * SCoM(1) * SCoM(3))) - ...
        (sys_neg_XZ_POI_origin + (sys_neg_mass * SCoM(1) * SCoM(3)));
    sys_YZ_POI_CoM = ...
        (sys_pos_YZ_POI_origin + (sys_pos_mass * SCoM(2) * SCoM(3))) - ...
        (sys_neg_YZ_POI_origin + (sys_neg_mass * SCoM(2) * SCoM(3)));
    
    % build CoM inertial tensor
    Stensor_com = ...
        [sys_X_MOI_CoM sys_XY_POI_CoM sys_XZ_POI_CoM; ...
        sys_XY_POI_CoM sys_Y_MOI_CoM sys_YZ_POI_CoM; ...
        sys_XZ_POI_CoM sys_YZ_POI_CoM sys_Z_MOI_CoM];
    
    
elseif runmode==2;
    
    % if no cavities present, just use pos values
    Smass = sys_pos_mass;
    SCoM = sys_pos_CoM;
    Stensor_origin = sys_pos_origin_tensor;

    % calculate MOIs about CoM
    sys_X_MOI_CoM = ...
        (sys_pos_X_MOI_origin - (sum(SCoM([2 3]).^2,2).* sys_pos_mass));
    sys_Y_MOI_CoM = ...
        (sys_pos_Y_MOI_origin - (sum(SCoM([1 3]).^2,2).* sys_pos_mass));
    sys_Z_MOI_CoM = ...
        (sys_pos_Z_MOI_origin - (sum(SCoM([1 2]).^2,2).* sys_pos_mass));   
    
    % calculate POIs about CoM
    sys_XY_POI_CoM = ...
        (sys_pos_XY_POI_origin + (sys_pos_mass * SCoM(1) * SCoM(2)));
    sys_XZ_POI_CoM = ...
        (sys_pos_XZ_POI_origin + (sys_pos_mass * SCoM(1) * SCoM(3)));
    sys_YZ_POI_CoM = ...
        (sys_pos_XZ_POI_origin + (sys_pos_mass * SCoM(2) * SCoM(3)));
    
    % build CoM inertial tensor
    Stensor_com = ...
        [sys_X_MOI_CoM sys_XY_POI_CoM sys_XZ_POI_CoM; ...
        sys_XY_POI_CoM sys_Y_MOI_CoM sys_YZ_POI_CoM; ...
        sys_XZ_POI_CoM sys_YZ_POI_CoM sys_Z_MOI_CoM];
    
end

end
