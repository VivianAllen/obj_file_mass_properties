function [] = makeControlBlank()
% generate a csv file that can be filled in with user data to control
% doMassPropertiesAnalysis.m

% build cell file
outData1 = cell(6,3);
outData1(:,1) = {'DENSITY IN kg m3:';'';...
    'UNITS (1 = m 0.01 = cm etc):';'';...
    'FIDELITY:';''};
outData2 = cell(3,3);
outData2(1,1) = {'MANUALLY ENTER SYSTEM ORIGIN BELOW (OR LEAVE BLANK):'};
outData2(2,:)= {'SYS ORIGIN X' 'SYS ORIGIN Y' 'SYS ORIGIN Z'};
outData3 = cell(2,3);
outData3(1,1) = {'SOLID OBJECTS:'};
outData3(1,2) = {'CAVITIES:'};
outData3(1,3) = {'SYSTEM ORIGIN OBJ:'};
outData = [outData1; outData2; outData3;];

% write cell data as csv
[rows, cols] = size(outData);
fid = fopen('BlankController.csv', 'w');
for i_row = 1:rows
    file_line = '';
    for i_col = 1:cols
        contents = outData{i_row, i_col};
        if isnumeric(contents)
            contents = num2str(contents);
        elseif isempty(contents)
            contents = '';
        end
        if i_col < cols
            file_line = [file_line, contents, ', '];
        else
            file_line = [file_line, contents];
        end
    end
    count = fprintf(fid, '%s\n', file_line);
end
st = fclose(fid);

display('BLANK CONTROLLER CSV FILE CREATED. FILL IN DENSITY, UNITS AND ')
display('FIDELITY (granularity setting, higher is better but takes longer, ')
display('50 is a good value), THEN LIST SOLID OBJECTS AND CAVITIES UNDER ')
display('THE RELEVANT HEADINGS. LEAVE CAVITIES COLUMN BLANK IF YOU HAVE NO ')
display('CAVITIES. IF YOU ENTER A SYSTEM ORIGIN MANUALLY LEAVE SYSTEM ORIGIN')
display(' OBJ COLUMN BLANK. LEAVE BOTH MANUAL SYSTEM ORIGIN ENTRY AND SYSTEM ')
display('ORIGIN COLUMN BLANK IF YOU JUST WANT TO USE THE WORLD ORIGIN')
end

