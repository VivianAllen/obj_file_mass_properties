function [mass, CoM, tensor_origin, tensor_com] ...
= getMassProps(Wobjname, density, fidelity, units)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   AUTHOR Vivian Allen (mrvivianallen@gmail.com)                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% VERSION: 2.4(ish) (09 February 2015)

% Version changes (2.5) - added option to specify a zero coordinate, which
% the mass properties will be expressed relative to. This is a hard-coded
% switch, rather than an option (laziness)

% version changes (2.4) - fixed bug when dealing with small meshes where
% ray array would not be correct dimensions due to rounding error
% version changes (2.3) - switched output order of origing and com tensors 
% to match older code
% version changes (2.2) - removed density scaling - scales mesh to metres 
% instead (makes more sense to have output always be in SI)
% version changes (2.1) - switched from line-plane intersection method of 
% ray/poly intersection detection to logic and edge-based homebrew version. 
% Quicker, and doesn't need any 3rd party code.

% populates a geometric mesh body (vertlist v, facelist f) with cuboid
% elements using a ray tracing method - i.e. plots an array of lines 
% through the longest dimension of the mesh, works out where each line 
% enters and exits the mesh, and then uses the 'internal' part of that line
% as the basis for a cuboid element.

% SPECIAL VERSION - USESE PRE-EXISTING V,F INSTEAD OF LOADING OBJ FILE,
% ASSUMES THAT V IS ALREADY SCALED TO METRES

%%%%%%%%%%%%%%%%%%%% INPUT EXPLANATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 	v = vertlist                                                          %
%   f = facelist                                                          %
% 	density = density of body in kgs per m^3 (1060 is kinda standard)     % 
% 	fidelity = control for size of ray array (array will be fidility^2    %
% 	rays). Higher is more accurate, exponentially slower. 50 seems good.  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%% OUTPUT EXPLANATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 	mass = estimated mass of target (kg)                                  %
% 	CoM = estimated CoM of target (m)                                     % 
% 	tensor_com = inertial tensor about CoM of target                      %
% 	tensor_origin = inertial tensor about origin of target                %
% (NB use ~ to suppress output of any variable)                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% REQUIRED FUNCTIONS (must be in same folder):

% NONE! 

%% DEBUG MODE CELL
% set debug mode to '1', set input options below, and run cell-by-cell 
% like a script) to debug. Debug mode plots the loaded mesh and the
% estimated intersecting polygons and intersect points for each ray
% plotted. Debug Figure will plot successive ray starts/stops as '*' with
% ray number below them, and ray/polygon intersects as 'o' with polygon 
% number. While each ray and intersect group is plotted, the vertices of
% the relevant polygons will flash as green '*'. The speed of this is
% controlled by the variable 'debug pause' below.


% set full loop if not in debug
loop_start = 1; % ray to start debug loop
loop_stop = fidelity^2; % ray to stop debug loop

    
% EXTRA SETTINGS
identity_threshold = 1/100000; % this sets how close intersect points can be
% before they get merged to avoid duplications arising from rare exact
% ray/edge interactions. Due to rounding errors (probably), 'unique' does 
% not catch all of these duplications

doZero = 0; % specify a zero (rather than the world one) for mass 
% properties expression
zeroCoord = [-0.037	-2.049	0.463];
%% LOAD WOBJ FILE
%ADAPTED FROM: 
%read_vertices_and_faces_from_obj_file(filename)
%BY: Alex Jacobson (http://www.alecjacobson.com/weblog/?p=917)
% Reads a .obj mesh file and outputs the vertex and face list
% assumes a 3D triangle mesh and ignores everything but:
% v x y z and f i j k lines: NB VIV EDITS TO ALSO DEAL WITH QUADS (NOT
% TOO HARD!)
%
v = zeros(0,3);
f = zeros(0,4); %ADDED AN EXTRA VALUE HERE (WORKS WITH QUADS)
vertex_index = 1;
face_index = 1;
fid = fopen(Wobjname,'rt');
line = fgets(fid);
while ischar(line)
    vertex = sscanf(line,'v %f %f %f');
    face = sscanf(line,'f %d %d %d %d'); % ADDED AN EXTRA VALUE HERE (WORKS WITH QUADS)
    face_long = sscanf(line,'f %d/%d/%d %d/%d/%d %d/%d/%d');
    % see if line is vertex command if so add to vertices
    if(size(vertex)>0)
        v(vertex_index,:) = vertex;
        vertex_index = vertex_index+1;
        % see if line is simple face command if so add to faces
    elseif(size(face)>0)
        f(face_index,1:size(face)) = face; %ADDED DYNAMIC INDEX size(face) TO AVOID DIMENSION MISMATCH
        face_index = face_index+1;
        % see if line is a long face command if so add to faces
    elseif(size(face_long)>0)
        % remove normal and texture indices
        face_long = face_long(1:3:end);
        f(face_index,:) = face_long;
        face_index = face_index+1;
    end
    line = fgets(fid);
end
fclose(fid);
  
% IF TRIS NOT QUADS, STRIPS LAST COLUMN (WILL BE ALL ZEROS) FROM FACE
% LIST, f
if sum(f(:,4))==0
    f(:,4)=[];
end

% IF NOT TRIS OR QUADS GIVE ERROR
if size(f,2)>4
    error('FACES LIST MORE THAN FOUR VERTS - NOT A TRI OR QUAD MESH, CANNOT BE PROCESSED')
end

%% ZERO IF WANTED
if doZero==1;
    zeroMat = repmat(zeroCoord,size(v,1),1);
   v = v -  zeroMat;
end

%% CONVERT UNITS

% convert coordinates to metres
v = v.*units;

%% GET RAY PARAMETERS

% determine body dimensions 
body_box = [max(v(:,1)) max(v(:,2)) max(v(:,3)); ...
     min(v(:,1)) min(v(:,2)) min(v(:,3))];
body_dims = [(body_box(1,1)-body_box(2,1)) (body_box(1,2)-body_box(2,2)) (body_box(1,3)-body_box(2,3))];

% get longest dimension
[~,ray_dim] = max(body_dims);

% check for equally large dimensions, use first max equal dimension if 
% found
if sum(body_dims==body_dims(ray_dim))~=1
    ray_dim = find((body_dims==ray_dim)==1, 1,'first');
end

%% BUILD FIDELITY * FIDELITY (2D) ARRAY OF RAYS POINTING ALONG MAX DIMENSION

array_dims = [1 2 3];
array_dims(ray_dim) = [];

array_1_unit = body_dims(array_dims(1))/(fidelity);
array_1_start = body_box(2,array_dims(1)) + (array_1_unit/2);
array_1_stop = body_box(1,array_dims(1)) - (array_1_unit/2);
array_1 = array_1_start:array_1_unit:array_1_stop;
% EXCEPTION CATCH - rounding error can sometimes cause array setup to stop
% before reaching proper stop point
if length(array_1)<fidelity
    array_1 = [array_1 array_1_unit];
end

array_2_unit = body_dims(array_dims(2))/(fidelity);
array_2_start = body_box(2,array_dims(2)) + (array_2_unit/2);
array_2_stop = body_box(1,array_dims(2)) - (array_2_unit/2);
array_2 = array_2_start:array_2_unit:array_2_stop;
% EXCEPTION CATCH - rounding error can sometimes cause array setup to stop
% before reaching proper stop point
if length(array_2)<fidelity
    array_2 = [array_2 array_2_unit];
end

rays = zeros(fidelity^2,2);

loop1_count = 1;
for loop1 = fidelity:fidelity:fidelity^2
rays(loop1-fidelity+1:loop1,1) = array_1;
rays(loop1-fidelity+1:loop1,2) = array_2(loop1_count);
loop1_count = loop1_count+1;
end

%% MAIN LOOP: GO THROUGH RAYS, FIND POLY INTERSECTIONS, POPULATE VOLUME

% intialise element mass storage
element_masses = 0;

% intialise element CoM storage
element_CoMs = [0 0 0];

% intialise element MOI storage
element_MOIs = [0 0 0];

% get polygon bounding boxes (makes searching faster
f_ids = 1:size(f,1);
searchable_f = [f_ids' f];
if sum(sum(f==0))>0;
    % if compound mesh, process tris and quads seperately
    % get tris
    trindex = searchable_f(:,end)==0;
    tris = searchable_f(trindex,:);
    tris(:,end) = [];
    % get tri verts
    tri_verts = zeros(size(tris,1),9);
    tri_verts(:,1:3) = v(tris(:,2),:);
    tri_verts(:,4:6) = v(tris(:,3),:);
    tri_verts(:,7:9) = v(tris(:,4),:);
    % get tri vert max mins
    tri_maxX = max(tri_verts(:,1:3:end),[],2);
    tri_maxY = max(tri_verts(:,2:3:end),[],2);
    tri_maxZ = max(tri_verts(:,3:3:end),[],2);
    tri_minX = min(tri_verts(:,1:3:end),[],2);
    tri_minY = min(tri_verts(:,2:3:end),[],2);
    tri_minZ = min(tri_verts(:,3:3:end),[],2);
    % assemble as bounding box
    tri_boxes = [tris(:,1) tri_maxX tri_maxY tri_maxZ tri_minX tri_minY...
        tri_minZ];
   
    % get quads
    quadex = searchable_f(:,end)>0;
    quads = searchable_f(quadex,:);
    % get quad verts
    quad_verts = zeros(size(quads,1),12);
    quad_verts(:,1:3) = v(quads(:,2),:);
    quad_verts(:,4:6) = v(quads(:,3),:);
    quad_verts(:,7:9) = v(quads(:,4),:);
    quad_verts(:,10:12) = v(quads(:,5),:);
    % get quad vert max mins
    quad_maxX = max(quad_verts(:,1:3:end),[],2);
    quad_maxY = max(quad_verts(:,2:3:end),[],2);
    quad_maxZ = max(quad_verts(:,3:3:end),[],2);
    quad_minX = min(quad_verts(:,1:3:end),[],2);
    quad_minY = min(quad_verts(:,2:3:end),[],2);
    quad_minZ = min(quad_verts(:,3:3:end),[],2);
    % assemble as bounding box
    quad_boxes = [quads(:,1) quad_maxX quad_maxY quad_maxZ quad_minX...
        quad_minY quad_minZ];
    
    % assemble and sort as final polybox
    poly_boxes = [tri_boxes; quad_boxes];
    poly_boxes = sortrows(poly_boxes,1);
    
else
    
    % if homegenous mesh, process all
    % get tri verts
    all_poly_verts = zeros(size(f,1),(size(f,2)*3));
    all_poly_verts(:,1:3) = v(f(:,1),:);
    all_poly_verts(:,4:6) = v(f(:,2),:);
    all_poly_verts(:,7:9) = v(f(:,3),:);
    if size(f,2)==4
        all_poly_verts(:,10:12) = v(f(:,4),:);
    end
    % get all_poly_vert max mins
    all_poly_maxX = max(all_poly_verts(:,1:3:end),[],2);
    all_poly_maxY = max(all_poly_verts(:,2:3:end),[],2);
    all_poly_maxZ = max(all_poly_verts(:,3:3:end),[],2);
    all_poly_minX = min(all_poly_verts(:,1:3:end),[],2);
    all_poly_minY = min(all_poly_verts(:,2:3:end),[],2);
    all_poly_minZ = min(all_poly_verts(:,3:3:end),[],2);
    % assemble as bounding box
    poly_boxes = [searchable_f(:,1) all_poly_maxX all_poly_maxY...
        all_poly_maxZ all_poly_minX all_poly_minY all_poly_minZ];
end

% we will only be searching the bounding boxes in the ray array plane, so
% remove the box elements in the ray dimensions itself
delindex = zeros(1,7);
delindex([ray_dim+1 ray_dim+1+3]) = 1;
poly_boxes(:,logical(delindex)) = [];

for loop1 = loop_start:loop_stop;
    
    % get ray
    ray = rays(loop1, :);
    
    % identify intersecting polygons - first stage (bounding box search)
    inter_polys1 = poly_boxes((poly_boxes(:,2)>ray(1) & poly_boxes(:,4)<ray(1)) & ...
        (poly_boxes(:,3)>ray(2) & poly_boxes(:,5)<ray(2)),1);
        
    % if polys found, proceed, check all properly for intersections, and
    % get intersection points if any
    if isempty(inter_polys1)~=1
    
        % prep/reset storage
        inter_points = [0 0 0];
        write2 = 1;
        
        for loop2 = 1:length(inter_polys1)
    
            % get poly
            poly = f(inter_polys1(loop2), :);
            % trim zeros (tris in a list also containing quads)
            poly(poly==0) = [];
            % get poly verts
            poly_verts = v(poly,:);
    
            % secondary detection - does ray intersect poly itself?
            % slice-like methodology: Rationale
            %   - treat ray (point) as slice horizon @ ray's first array dim
            %   - identify polys edges with verts both above and below slice 
            %   horizon (only two in quads/tris)
            %   - use simple constant to find point where slice hits edges
            %   - use simple constant to find point where line between slice/edge
            %   interactions == ray.
            %   - if this point is not between slice/edge interactions, ray is
            %   not intersecting with polygon
    
            poly_verts_horproj = poly_verts(:, array_dims(1));
            % get vert numbers above and below, first as index in poly list
            [verts_abovedex,~] = find((poly_verts_horproj>=ray(1))==1);
            [verts_belowdex,~] = find((poly_verts_horproj<=ray(1))==1);
            % then as vert numbers
            verts_above = poly(verts_abovedex);
            verts_below = poly(verts_belowdex);
    
            % use logic on poly vert order to turn into edge pairs
            % if either above or below is a singleton, procedure is simple. The
            % pairs will be the verts on either side in the poly list (poly
            % list is circular, so first and last entries are special cases
            if length(verts_above)==1 || length(verts_below)==1
                if length(verts_above)==1
                    vert1dex = verts_abovedex;
                    vert3dex = verts_abovedex;
                else
                    vert1dex = verts_belowdex;
                    vert3dex = verts_belowdex;
                end
                if vert1dex==1
                    vert2dex = length(poly);
                    vert4dex = 2;
                elseif vert1dex==length(poly)
                    vert2dex = length(poly)-1;
                    vert4dex = 1;
                else
                    vert2dex = vert1dex+1;
                    vert4dex = vert1dex-1;
                end
            else 
                % if above and below are binary two things are true -
                % 1, the poly is a quad, 2, it has been sliced in two. No other
                % outcomes are possible because this code will have stopped if
                % anything other than a tri or quad is present. Doesnt matter
                % if we start above or below, so start above (christian-like)
                vert1dex = verts_abovedex(1);
                vert3dex = verts_abovedex(2);
                % the corresponding below vert must be either one up or down on
                % poly list: test (NB circular, first and last special)
                if vert1dex==1 
                    % if first is one, second must be either second or last
                    vert2dex = verts_belowdex(verts_belowdex==vert1dex+1 |...
                        verts_belowdex==length(poly));
                elseif vert1dex==length(poly); 
                    % if first is last, second must be either second from last
                    % or first
                    vert2dex = verts_belowdex(verts_belowdex==vert1dex-1 |...
                        verts_belowdex==1);
                else
                    % else it must be either up or down one
                    vert2dex = verts_belowdex(verts_belowdex==vert1dex+1 |...
                        verts_belowdex==vert1dex-1);
                end
                if vert3dex==1 
                    % if first is one, second must be either second or last
                    vert4dex = verts_belowdex(verts_belowdex==vert3dex+1 |...
                        verts_belowdex==length(poly));
                elseif vert3dex==length(poly); 
                    % if first is last, second must be either second from last
                    % or first
                    vert4dex = verts_belowdex(verts_belowdex==vert3dex-1 |...
                        verts_belowdex==1);
                else
                    % else it must be either up or down one
                    vert4dex = verts_belowdex(verts_belowdex==vert3dex+1 |...
                        verts_belowdex==vert3dex-1);
                end                   
            end
            % now we have two pairs of edges that cross ray horizon
            edge1 = [poly_verts(vert1dex,:); poly_verts(vert2dex,:)];
            edge2 = [poly_verts(vert3dex,:); poly_verts(vert4dex,:)];
    
            % find the points on each edge equal to the ray horizon
    
            % get ratio of size-to-horizon to total edge size in horizon axis
            edge1_factor = (edge1(1, array_dims(1)) - ray(1)) / ...
                (edge1(1, array_dims(1)) - edge1(2, array_dims(1))); 
            edge2_factor = (edge2(1, array_dims(1)) - ray(1)) / ...
                (edge2(1, array_dims(1)) - edge2(2, array_dims(1)));      
    
            % get vector from edge start to intersect point
            edge1_intervec = (edge1(2,:) - edge1(1,:)) * edge1_factor;
            edge2_intervec = (edge2(2,:) - edge2(1,:)) * edge2_factor;
    
            % get intesect points
            edge1_point = edge1(1,:) + edge1_intervec;
            edge2_point = edge2(1,:) + edge2_intervec;
    
            % INTERSECTION TEST TWO: DO THE TWO POINTS BOUND THE RAY POINT?
            testmat =  sort([edge1_point(array_dims(2)) edge2_point(array_dims(2))]); 
            if ray(2)>=testmat(1) && ray(2)<=testmat(2)
    
                % IF SO, FIND INTERSECTION
                ray_factor = (edge1_point(array_dims(2)) - ray(2)) / ...
                (edge1_point(array_dims(2)) - edge2_point(array_dims(2)));
                ray_intervec = (edge2_point - edge1_point) * ray_factor;
                ray_point = edge1_point + ray_intervec;
    
                % write point
                inter_points(write2,:) = ray_point;
                write2 = write2+1;                               
                            
            end
                       
        end
        
        % set proceed (overwrite if condition met)
        proceed = 0;
        
        % ERROR CHECKING - if any intersection points found, check for
        % dupes
        if write2~=1
            
            % sort intersecting points along ray axis
            inter_points = sortrows(inter_points, ray_dim);
            
            % EXCEPTION CATCH: if the ray coincides exactly with a polygon
            % edge, you will get duplication in the interpoints list, which is
            % bad. Simple to remove though.
            % ah.. maybe not THAT simple to remove. unique does not catch
            % everything. Use threshold value instead...
            dupedex = diff(inter_points(:,ray_dim))<identity_threshold;
            inter_points(dupedex,:) = [];
            
            % EXCEPTION CATCH: if the ray only intersects with one edge,
            % i.e just 'grazes' mesh (could happen at very high fidelity)
            % then duplicate check could leave singleton intersect point.
            % check for this, and if it happens, short loop. This can be
            % differentiated from simply hitting part of a non-manifold
            % mesh by the presence of duplicates (i.e. non-zero dupedex
            % entries), as hopefully (hopefully) these will not occur on
            % the same ray
            % reset proceed (overwrite if condition met)
            proceed = 1;
            if sum(dupedex(dupedex==1))~=0
                if mod(size(inter_points,1),2)==1
                    proceed = 0; % do not proceed
                end
            end
        end
            
        % if any intersection points where found, proceed to build cuboids
        if proceed==1
            
            % use intersecting points to build cuboid(s)
            
            % error state - if there are an odd number of intersection points, the
            % ray must entered the mesh and did not leave, i.e. the MESH IS NON
            % MANIFOLD
            if mod(size(inter_points,1),2)==1
                ray_address = zeros(1,3);
                ray_address(array_dims) = ray;
                ray_string = ['(', num2str(ray_address(1)),',', ...
                    num2str(ray_address(2)),',',num2str(ray_address(3)),')'];
                
                % draw mesh
                if sum(sum(f==0))>0;
                    compoundMesh = 1;
                    [rows,~] = find(f==0);
                    tris = f(rows,:);
                    tris(:,4) = [];
                    quads = f;
                    quads(rows,:) = [];
                else
                    compoundMesh = 0;
                end
                
                figure(1)
                clf
                hold on
                if compoundMesh ==0;
                    trimesh(f,v(:,1), v(:,2), v(:,3),'edgecolor','m', 'facecolor','none')
                else
                    trimesh(quads,v(:,1), v(:,2), v(:,3),'edgecolor','m', 'facecolor','none')
                    trimesh(tris,v(:,1), v(:,2), v(:,3),'edgecolor','m', 'facecolor','none')
                end
                view ([0,0])
                axis('equal')
                
                % draw ray
                ray_plot = zeros(2,3);
                ray_plot(1,array_dims) = ray;
                ray_plot(2,array_dims) = ray;
                ray_plot(1,ray_dim) = body_box(2,ray_dim)-(ray_dim*0.25);
                ray_plot(2,ray_dim) = body_box(1,ray_dim)+(ray_dim*0.25);
                plot3(ray_plot(:,1), ray_plot(:,2), ray_plot(:,3), '-','color','k')
                text(ray_plot(1,1), ray_plot(1,2), ray_plot(1, 3), num2str(loop1), 'VerticalAlignment', 'top')
                
                title('ERROR IN MESH PROCESSING')
                
                error(['ERROR PROCESSING ',Wobjname,' RAY ',num2str(loop1),' AT ', ray_string,...
                    ' ENTERED MESH AND DID NOT LEAVE, I.E MESH IS NON-MANIFOLD (NOT WATERTIGHT)'])                
            end
            
            % break into pairs (inside/outside) that define length of cuboid
            % elements. Other dimensions are determined by the array setup.
    
            start_index = 1:2:size(inter_points,1);
            stop_index = 2:2:size(inter_points,1);
            % get element lengths
            element_lengths = sqrt(sum(abs(inter_points(start_index,:) -  inter_points(stop_index,:)).^2,2));
    
            % matrix-asize-ationate
            length_matrix = zeros(length(element_lengths),3);
            length_matrix(:,ray_dim) = element_lengths/2;
    
            % apppend element CoMs
            TEMP_CoMs = ((inter_points(start_index,:) + length_matrix) + ...
                (inter_points(stop_index,:) - length_matrix))./2;
    
            element_CoMs = [element_CoMs; TEMP_CoMs];
    
    
            % append element masses
            TEMP_masses = element_lengths * array_1_unit * array_2_unit * density;
            element_masses = [element_masses; TEMP_masses];
    
            % append principle element moments of inertia (cuboids, so no products)
            % NB - these are about element CoM
            % taken from: http://en.wikipedia.org/wiki/List_of_moments_of_inertia
            TEMP_MOIs = zeros(length(element_lengths),3);
            TEMP_MOIs(:,ray_dim) = (1/12).*TEMP_masses.*(array_1_unit^2 + array_2_unit^2);
            TEMP_MOIs(:,array_dims(1)) = (1/12).*TEMP_masses.*(element_lengths.^2 + array_2_unit^2);
            TEMP_MOIs(:,array_dims(2)) = (1/12).*TEMP_masses.*(element_lengths.^2 + array_1_unit^2);
            element_MOIs = [element_MOIs; TEMP_MOIs];
        end
    
    end
    
end

% trim!
element_masses(1) = [];
element_CoMs(1,:) = [];
element_MOIs(1,:) = [];

%% SUM BODY MASS PROPERTIES

% sum mass
body_mass = sum(element_masses);

% calculate CoM
body_CoM = (1/body_mass) * sum(element_CoMs .* ...
    [element_masses element_masses element_masses]);

% calculate MOIs about origin
body_X_MOI_origin = sum(element_MOIs(:,1) + (sum(element_CoMs(:,[2 3]).^2,2) .* element_masses));
body_Y_MOI_origin = sum(element_MOIs(:,2) + (sum(element_CoMs(:,[1 3]).^2,2) .* element_masses));
body_Z_MOI_origin = sum(element_MOIs(:,3) + (sum(element_CoMs(:,[1 2]).^2,2) .* element_masses));

% calculate POIs about origin
body_XY_POI_origin = sum((element_masses .* element_CoMs(:,1) .* element_CoMs(:,2))*-1); 
body_XZ_POI_origin = sum((element_masses .* element_CoMs(:,1) .* element_CoMs(:,3))*-1); 
body_YZ_POI_origin = sum((element_masses .* element_CoMs(:,2) .* element_CoMs(:,3))*-1); 

% build origin inertial tensor
body_origin_tensor = ...
    [body_X_MOI_origin body_XY_POI_origin body_XZ_POI_origin; ...
    body_XY_POI_origin body_Y_MOI_origin body_YZ_POI_origin; ...
    body_XZ_POI_origin body_YZ_POI_origin body_Z_MOI_origin];

% calculate MOIs about CoM
body_X_MOI_CoM = body_X_MOI_origin - (sum(body_CoM([2 3]).^2,2).* body_mass);
body_Y_MOI_CoM = body_Y_MOI_origin - (sum(body_CoM([1 3]).^2,2).* body_mass);
body_Z_MOI_CoM = body_Z_MOI_origin - (sum(body_CoM([1 2]).^2,2).* body_mass);

% calculate POIs about CoM
body_XY_POI_CoM = body_XY_POI_origin + (body_mass * body_CoM(1) * body_CoM(2));
body_XZ_POI_CoM = body_XZ_POI_origin + (body_mass * body_CoM(1) * body_CoM(3));
body_YZ_POI_CoM = body_YZ_POI_origin + (body_mass * body_CoM(2) * body_CoM(3));

% build CoM inertial tensor
body_CoM_tensor = ...
    [body_X_MOI_CoM body_XY_POI_CoM body_XZ_POI_CoM; ...
    body_XY_POI_CoM body_Y_MOI_CoM body_YZ_POI_CoM; ...
    body_XZ_POI_CoM body_YZ_POI_CoM body_Z_MOI_CoM];

%% SET AND DISPLAY RESULTS

mass = body_mass;
CoM = body_CoM;
tensor_com = body_CoM_tensor;
tensor_origin = body_origin_tensor;