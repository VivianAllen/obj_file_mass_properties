%% SCRIPT DESCRIPTION
% control script for running mass properties analysis on .obj files, either
% singly or as grouped segments

% bug fixes 06/04/17 - issue with blank world origin fixed


% you can set if you want it draw the output here. I'm assuming you do.
drawOut = 'y';

%% PREMADE CONTROLS

useControlFile = input('DO YOU WANT TO USE A PREMADE CONTROLLER FILE? (y/n) ','s');
doManual = 1;

if useControlFile == 'y'
    
    doManual = 0;
    
    % get list of all csv files in current folder
    folderStruct = dir('*.csv');
    
    % convert struct to list
    csvList = cell(size(folderStruct,1),1);
    for loop1 = 1:size(folderStruct,1)
        csvList(loop1) = {folderStruct(loop1).name};
    end
    
    clear folderStruct
    
    csvDispList = csvList;
    
    for loop1 = 1:length(csvDispList);
        csvDispList(loop1) = {[num2str(loop1), '. ', cell2mat(csvList(loop1))]};
    end
    
    display ('CSV FILES CURRENTLY IN FOLDER')
    display (csvDispList)
    disp ' '
    
    csvChoice = input('WHICH FILE NUMBER DO YOU WANT TO USE? ');
    
    % load control file
    controlData = csv2cell(cell2mat(csvList(csvChoice)),'fromfile');
    
    % edit name
    analName = cell2mat(csvList(csvChoice));
    analName = analName(1:end-4);
    analName(strfind(analName,'CONTROLLER'):end) = [];
    
    % parse control file
    for loop1 = 1:length(controlData)
        if strmatch('DENSITY IN kg m3:',controlData(loop1,1));
            density = str2num(cell2mat(controlData(loop1+1,1)));
        elseif strmatch('UNITS (1 = m 0.01 = cm etc):',controlData(loop1,1));
            units = str2num(cell2mat(controlData(loop1+1,1)));
        elseif strmatch('FIDELITY:',controlData(loop1,1));
            fidelity = str2num(cell2mat(controlData(loop1+1,1)));
        elseif strmatch('MANUALLY ENTER SYSTEM ORIGIN BELOW (OR LEAVE BLANK):',controlData(loop1,1));
            system_origin = [str2num(cell2mat(controlData(loop1+2,1)))...
                str2num(cell2mat(controlData(loop1+2,2))) str2num(cell2mat(controlData(loop1+2,3)))];
        elseif strmatch('SOLID OBJECTS:',controlData(loop1,1));
            solLoadList = controlData(loop1+1:end,1);
            cavLoadList = controlData(loop1+1:end,2);
            originLoad = controlData(loop1+1:end,3);
        end
    end
    
    % if manually entered system origin is blank, clear it
    if isempty(system_origin)
        clear system_origin
    end
    
    % check both lists for entries that don't end in .obj
    solDelDex = logical(length(solLoadList));
    cavDelDex = logical(length(cavLoadList));
    originDelDex = logical(length(originLoad));
    
    for loop1 = 1:length(solLoadList)
        solDelDex(loop1) = isempty(cell2mat(strfind(solLoadList(loop1),'obj')));
        cavDelDex(loop1) = isempty(cell2mat(strfind(cavLoadList(loop1),'obj')));
        originDelDex(loop1) = isempty(cell2mat(strfind(originLoad(loop1),'obj')));
    end
    
    % delete bad entries
    solLoadList(solDelDex) = [];
    cavLoadList(cavDelDex) = [];
    originLoad(originDelDex) = [];
    
    % delete cavLoadList if blank
    if isempty(cavLoadList)
        clear cavLoadList
    end
    if isempty(originLoad)
        clear originLoad
    end
           
    disp 'RUNNING ANALYSIS WITH THESE PARAMETERS:'
    disp(['NAME: ', analName])
    disp(['DENSITY: ', num2str(density)])
    disp(['UNITS: ', num2str(units)])
    disp(['FIDELITY: ', num2str(fidelity)])
    disp ' '
    display 'SOLID LIST:'
    display(solLoadList)
    display 'CAVITY LIST:'
    if ~exist('cavLoadList','var')
        disp ' '
    else
        display(cavLoadList)
    end
    disp ' '
    if exist('system_origin','var')
        disp 'USING MANUALLY ENTERED SYSTEM ORIGIN:'
        system_origin
    else
        if ~exist('originLoad','var')
            disp 'USING WORLD ZERO AS SYSTEM ORIGIN'
        else
            disp(['USING CENTRE OF BOUNDING BOX OF ', cell2mat(originLoad), ' AS SYSTEM ORIGIN']);
        end
    end
    
end
%% MANUAL INPUT

if doManual==1;
    
    % name of analysis
    analName = input('NAME ANALYSIS (SETS OUTPUT FILENAME): ','s');
    
    % density
    density = input('SET DENSITY FOR SOLID OBJECTS (kg m3): ');
    
    % units of meshes to be analysed
    units = input('SET UNITS (AS A NUMBER) OF INPUT MESH FILES (1 = m, 0.01 = cm, etc): ');
    
    % fidelity of mass properties analysis
    display('SET FIDELITY OF MASS PROPERTIES ANALYSIS (granularity setting,')
    fidelity = ...
        input('higher is better but takes longer, 50 is a good value): ');
    
    
    disp 'RUNNING ANALYSIS WITH THESE PARAMETERS:'
    disp(['NAME: ', analName])
    disp(['DENSITY: ', num2str(density)])
    disp(['UNITS: ', num2str(units)])
    disp(['FIDELITY: ', num2str(fidelity)])
    disp ' '
    
    % GET AND DISPLAY LIST ALL OBJS IN CURRENT FOLDER
    
    folderStruct = dir('*.obj');
    
    % convert struct to list
    objList = cell(size(folderStruct,1),1);
    for loop1 = 1:size(folderStruct,1)
        objList(loop1) = {folderStruct(loop1).name};
    end
    
    clear folderStruct
    
    dispList = objList;
    
    for loop1 = 1:length(dispList);
        dispList(loop1) = {[num2str(loop1), '. ', cell2mat(objList(loop1))]};
    end
    
    % CHOOSE SOLID OBJECTS AND CAVITIES
    proceed = 'n';
    
    while proceed=='n';
        
        % display list of .obj files
        display ('OBJ FILES CURRENTLY IN FOLDER')
        display (dispList)
        disp ' '
        % get list as string (simplest way of doing it)
        solListStr = input('LIST FILES (BY NUMBER) TO BE TREATED AS SOLID OBJECTS: ','s');
        
        % convert string to list of numbers
        solList =  cell2mat(textscan(solListStr,'%f'));
        solLoadList = objList(solList);
        
        % CHOOSE CAVITY OBJECTS LIST
        % get list as string (simplest way of doing it)
        display('LIST FILES (BY NUMBER) TO BE TREATED AS CAVITIES ')
        cavListStr = input('(Just press enter if there are no cavities):  ','s');
        
        % convert string to list of numbers
        if ~isempty(cavListStr)
            cavList =  cell2mat(textscan(cavListStr,'%f'));
            cavLoadList = objList(cavList);
        end
        
        % CHOOSE SYSTEM ORIGIN
        % get list as string (simplest way of doing it)
        display('CHOOSE FILE (BY NUMBER) TO BE USED FOR SYSTEM ORIGIN ')
        display('ALTERNATELY JUST PRESS ENTER TO USE WORLD ORIGIN  ')
        originListStr = input('OR ENTER m TO MANUALLY ENTER ORIGIN:  ','s');
        
        % if manual entry, ask for manual entry as string and convert
        if ~isempty(originListStr) && originListStr=='m'
            display('MANUAL INPUT SELECTED')
            originListStr2 = input('ENTER SYSTEM ORIGIN: ','s');
            system_origin =  cell2mat(textscan(originListStr2,'%f'))';
        end
        % if not manual entry, convert string to list of numbers
        if ~isempty(originListStr) && originListStr~='m'
            originList =  cell2mat(textscan(originListStr,'%f'));
            originLoad = objList(originList);
        end
        
        
        disp ' '
        display 'SOLID LIST:'
        display(solLoadList)
        display 'CAVITY LIST:'
        if ~exist('cavLoadList','var')
            disp ' '
        else
            display(cavLoadList)
        end
        disp ' '        
        if exist('system_origin','var')
        disp 'USING MANUALLY ENTERED SYSTEM ORIGIN:'
        system_origin
            else
            if ~exist('originLoad','var')
                disp 'USING WORLD ZERO AS SYSTEM ORIGIN'
            else
                disp(['USING CENTRE OF BOUNDING BOX OF ', cell2mat(originLoad), ' AS SYSTEM ORIGIN']);
            end
        end
        disp ' '
        proceed = input('PROCEED WITH ANALYSIS USING ABOVE SETTINGS? (y/n) ','s');
        
        % build cell files
        outData1 = cell(6,3);
        outData1(:,1) = {'DENSITY IN kg m3:'; density;...
            'UNITS (1 = m 0.01 = cm etc):'; units;...
            'FIDELITY:'; fidelity};
        outData2 = cell(3,3);
        outData2(1,1) = {'MANUALLY ENTER SYSTEM ORIGIN BELOW (OR LEAVE BLANK):'};
        outData2(2,:)= {'SYS ORIGIN X' 'SYS ORIGIN Y' 'SYS ORIGIN Z'};
        if exist('system_origin','var')
            outData2(3,:) = num2cell(system_origin);
        end
        if ~exist('cavLoadList','var')
            maxLength = length(solLoadList);
        else
            maxLength = max([length(solLoadList) length(solLoadList)]);
        end
        outData3 = cell(maxLength+1,3);
        outData3(1:length(solLoadList)+1,1) = [{'SOLID OBJECTS:'}; solLoadList];
        if ~exist('cavLoadList','var')
            outData3(1,2) = {'CAVITIES:'};
        else
            outData3(1:length(cavLoadList)+1,2) = [{'CAVITIES:'}; cavLoadList];
        end
        if ~exist('originLoad','var')
            outData3(1,3) = {'SYSTEM ORIGIN OBJ:'};
        else
            outData3(1:2,3) = [{'SYSTEM ORIGIN OBJ:'}; originLoad];
        end
        
        % write controller data as csv
        controller = [outData1; outData2; outData3];
        [rows, cols] = size(controller);
        fid = fopen([analName, 'CONTROLLER.csv'], 'w');
        for i_row = 1:rows
            file_line = '';
            for i_col = 1:cols
                contents = controller{i_row, i_col};
                if isnumeric(contents)
                    contents = num2str(contents);
                elseif isempty(contents)
                    contents = '';
                end
                if i_col < cols
                    file_line = [file_line, contents, ','];
                else
                    file_line = [file_line, contents];
                end
            end
            count = fprintf(fid, '%s\n', file_line);
        end
        st = fclose(fid);
        
        disp(['CONTROLLER WRITTEN TO FILE ', analName,'CONTROLLER.csv'])
    end
    
end

%% DO ANALYSIS AND GET DATA
if ~exist('cavLoadList','var')
    [mass, CoM, ~, inertia_tensor_CoM] = getMassPropSystem...
        (solLoadList, density, fidelity, units);
else
    [mass, CoM, ~, inertia_tensor_CoM] = getMassPropSystem...
        (solLoadList, cavLoadList, density, fidelity, units);
end

%% ZERO IF NEEDED
if exist('originLoad','var') && ~exist('system_origin','var')
    %load vertices from file
    v = zeros(0,3);
    vertex_index = 1;
    fid = fopen(cell2mat(originLoad),'rt');
    line = fgets(fid);
    while ischar(line)
        vertex = sscanf(line,'v %f %f %f');
        % see if line is vertex command if so add to vertices
        if(size(vertex)>0)
            v(vertex_index,:) = vertex;
            vertex_index = vertex_index+1;
            % see if line is simple face command if so add to faces
        end
        line = fgets(fid);
    end
    fclose(fid);
    % get zeroing vector (centre of v cloud bounding box)
    zeroBox = [max(v(:,1)) max(v(:,2)) max(v(:,3)); ...
     min(v(:,1)) min(v(:,2)) min(v(:,3))];
    system_origin = mean(zeroBox);
end
 % make and apply zeroing vector if wanted
 if exist('system_origin','var')     
    zeroVec = system_origin*-1;
    CoM = CoM+zeroVec;
 else
     system_origin = [0 0 0];
 end

%% CALCULATE PRINCIPLE AXES

% HAH TURNS OUT IT'S EASY: calculate the axes and the principle moments
% using the eig function!
[pa_cols,principle_moments]=eig(inertia_tensor_CoM);
% eig outputs COLUMN vectors whereas I think ROWS are nicer, so lets
% reorganise...
principle_axes = [pa_cols(:,1)';pa_cols(:,2)';pa_cols(:,3)'];
% also might as well get a single matrix that has orientation and moments
% sensibly aligned
axes_norm = [repmat(norm(principle_axes(1,:)),1,3);
    repmat(norm(principle_axes(2,:)),1,3);
    repmat(norm(principle_axes(3,:)),1,3)]...
    ./[repmat(norm(principle_moments(1,:)),1,3);
    repmat(norm(principle_moments(2,:)),1,3);
    repmat(norm(principle_moments(3,:)),1,3)];

oriented_principle_axes = principle_axes./axes_norm;

%% DISPLAY DATA
system_origin
mass
CoM
inertia_tensor_CoM
oriented_principle_axes

%% SAVE DATA TO CSV
outData3 = [{'SYSTEM ORIGIN:' '' ''}; num2cell(system_origin)];
outData4 = {'MASS IN KG:' '' ''; mass '' ''};
outData5 = [{'CoM IN m FROM SYSTEM ORIGIN:' '' ''}; num2cell(CoM)];
outData6 = [{'INERTIAL TENSOR ABOUT CoM:'  '' ''}; num2cell(inertia_tensor_CoM)];
outData7 = [{'PRINCIPLE AXES OF INERTIA:'  '' ''}; num2cell(principle_axes)];
outData8 = [{'PRINCIPLE MOMENTS OF INERTIA:'  '' ''}; num2cell(principle_moments)];

% write results data as csv
results = [outData3; outData4; outData5; outData6; outData7; outData8];
[rows, cols] = size(results);
fid = fopen([analName, 'RESULTS.csv'], 'w');
for i_row = 1:rows
    file_line = '';
    for i_col = 1:cols
        contents = results{i_row, i_col};
        if isnumeric(contents)
            contents = num2str(contents);
        elseif isempty(contents)
            contents = '';
        end
        if i_col < cols
            file_line = [file_line, contents, ','];
        else
            file_line = [file_line, contents];
        end
    end
    count = fprintf(fid, '%s\n', file_line);
end
st = fclose(fid);

disp(['RESULTS WRITTEN TO FILE ', analName,'RESULTS.csv'])


%% DRAW RESULTS IF WANTED

if drawOut == 'y';
    
    figure(1)
    clf
    hold on
    
    % draw solids in red then cavities in blue
    if  ~exist('cavLoadList','var')
        plotList = solLoadList;
        plotColorDex = ones(length(plotList),1);
    else
        plotList = [solLoadList; cavLoadList];
        plotColorDex = [ones(length(solLoadList),1); zeros(length(solLoadList),1)];
    end
    plotColors = cell(length(plotColorDex),1);
    plotColors(plotColorDex==1) = {'m'};
    plotColors(plotColorDex==0) = {'c'};
    
    % plot meshes
    for loop1 = 1:length(plotList)
        
        currObj = cell2mat(plotList(loop1));
        
        % LOAD WOBJ FILE
        v = zeros(0,3);
        f = zeros(0,4); %ADDED AN EXTRA VALUE HERE (WORKS WITH QUADS)
        vertex_index = 1;
        face_index = 1;
        fid = fopen(currObj,'rt');
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
        
        % ZERO V IF NEEDED
        if exist('zeroVec','var') 
            v = v + repmat(zeroVec,size(v,1),1);
        end
        
        % IF TRIS NOT QUADS, STRIPS LAST COLUMN (WILL BE ALL ZEROS) FROM FACE
        % LIST, f
        if sum(f(:,4))==0
            f(:,4)=[];
        end
        
        % plot obj file
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
        hold on
        if compoundMesh ==0;
            trimesh(f,v(:,1), v(:,2), v(:,3),'edgecolor',cell2mat(plotColors(loop1)), 'facecolor','none')
        else
            trimesh(quads,v(:,1), v(:,2), v(:,3),'edgecolor',cell2mat(plotColors(loop1)), 'facecolor','none')
            trimesh(tris,v(:,1), v(:,2), v(:,3),'edgecolor',cell2mat(plotColors(loop1)), 'facecolor','none')
        end
        view (3)
        axis('equal')
    end
    
    % get tensor axes normalised to a proportion of the max dimension of
    % the model
    plotTensor = zeros(3,3);
    plotTensor(1,1) = inertia_tensor_CoM(1,1);
    plotTensor(2,2) = inertia_tensor_CoM(2,2);
    plotTensor(3,3) = inertia_tensor_CoM(3,3);
    
    axesNormed = plotTensor./max(max(plotTensor));
    scalVals = [xlim' ylim' zlim'];
    scalFac = max(scalVals(2,:) - scalVals(1,:));
    axesScaled = axesNormed.*(scalFac*0.5);
    
    % draw CoM and inertia
    p1_plot = [CoM; CoM+axesScaled(1,:)];
    p2_plot = [CoM; CoM+axesScaled(2,:)];
    p3_plot = [CoM; CoM+axesScaled(3,:)];
    
    plot3(CoM(1),CoM(2),CoM(3),'o','color','k','MarkerFaceColor','r','MarkerSize',15)
    plot3(p1_plot(:,1),p1_plot(:,2),p1_plot(:,3),'color','r','LineWidth',3)
    plot3(p2_plot(:,1),p2_plot(:,2),p2_plot(:,3),'color','b','LineWidth',3)
    plot3(p3_plot(:,1),p3_plot(:,2),p3_plot(:,3),'color','g','LineWidth',3)
    
    % plot system origin
    plot3(0,0,0,'o','color','k','MarkerFaceColor','k','MarkerSize',6)
    
    % draw principle axes
    % scale oriented principle moments to a sensible display value
    axes_norm = oriented_principle_axes...
        ./max(max([repmat(norm(oriented_principle_axes(1,:)),1,3);
        repmat(norm(oriented_principle_axes(2,:)),1,3);
        repmat(norm(oriented_principle_axes(3,:)),1,3)]));
    normed_axes = axes_norm.*(scalFac*0.6);
    
    % draw princple axes
    pa1_plot = [CoM; CoM+normed_axes(1,:)];
    pa2_plot = [CoM; CoM+normed_axes(2,:)];
    pa3_plot = [CoM; CoM+normed_axes(3,:)];
    
    plot3(CoM(1),CoM(2),CoM(3),'o','color','k','MarkerFaceColor','r','MarkerSize',15)
    plot3(pa1_plot(:,1),pa1_plot(:,2),pa1_plot(:,3),'--r','LineWidth',2)
    plot3(pa2_plot(:,1),pa2_plot(:,2),pa2_plot(:,3),'--b','LineWidth',2)
    plot3(pa3_plot(:,1),pa3_plot(:,2),pa3_plot(:,3),'--g','LineWidth',2)
    
    title('System Origin = black dot, CoM = red dot, Solid Lines = inertia about global XYZ, Dotted Lines = inertia about principle axes')
    
end
%% HOUSEKEEPING
clear all

