%% Header
%
% Function creating all the permutations of names and relevant variables 
% for all the files requested based on the input file. It also checks for
% free inputs varrying without assigned names. Since the code works on a name
% basis, if a free variable has more than one variation it should be included
% inthe name.
%
% V.Pettas/F.Costa 
% University of Stuttgart, Stuttgart Wind Energy (SWE) 2019

function perm_cell = getNamesFromInputs(input)

nameBase = input.nameBase; 
freeInp  = input.freeInp;
fixedInp = input.fixedInp;
AllFixed = input.AllFixed;
HiddenIndex  = find(ismember(AllFixed,fixedInp(:,1))==0);
Hidden.Names = AllFixed(HiddenIndex);
%assign values to the fixed inputs cell vector
for iPat=1:numel(input.PatternNames)
    for i = 1:size(fixedInp,1)
        curNam = fixedInp{i};
        if strcmp(curNam,'Pat')
            fixedInp{i,2} = [iPat];
        elseif any(strcmp(Hidden.Names,'Pat')) && size(input.PatternNames,2) ==1
            index= find(strcmp(Hidden.Names,'Pat'));
            Hidden.Val{index,iPat} = input.PatternNames; %#ok<*FNDSB>
        else
            if ~isempty(HiddenIndex) && any(strcmp(Hidden.Names,'Pat'))
                error('There is an error with the input ''PatternNames'' ')
            end
        end
        if strcmp(curNam,'Tp')
            fixedInp{i,2} = [input.timestep_pat_vec{iPat}];
        elseif any(strcmp(Hidden.Names,'Tp')) && size(input.timestep_pat_vec{1},2) ==1
            index= find(strcmp(Hidden.Names,'Tp'));
            Hidden.Val{index,iPat} = input.timestep_pat_vec{iPat}; %#ok<*FNDSB>
        else
            if ~isempty(HiddenIndex)  && any(strcmp(Hidden.Names,'Tp'))
                error('There is an error with one of input ''timestep_pat_vec'' ')
            end
        end
        if strcmp(curNam,'Tm')
            fixedInp{i,2} = [input.timeStep_Measurements{iPat}];
        elseif any(strcmp(Hidden.Names,'Tm')) && length(input.timeStep_Measurements{1}) ==1
            index= find(strcmp(Hidden.Names,'Tm'));
            Hidden.Val{index,iPat} = input.timeStep_Measurements{iPat}; %#ok<*FNDSB>
        else
            if ~isempty(HiddenIndex) && any(strcmp(Hidden.Names,'Tm'))
                error('There is an error with the input ''timeStep_Measurements'' ')
            end
        end
        if strcmp(curNam,'Pos')
            fixedInp{i,2} = [input.Pos_LiDAR(1)]; %#ok<*NBRAK>
        elseif any(strcmp(Hidden.Names,'Pos')) && length(input.Pos_LiDAR) ==2
            index= find(strcmp(Hidden.Names,'Pos'));
            Hidden.Val{index,iPat} = input.Pos_LiDAR; %#ok<*FNDSB>
        else
            if ~isempty(HiddenIndex) && any(strcmp(Hidden.Names,'Pos'))
                error('There is an error with the input ''Pos_LiDAR'' ')
            end
        end
        if strcmp(curNam,'Fd')
            fixedInp{i,2} = [input.ref_plane_dist];
        elseif any(strcmp(Hidden.Names,'Fd')) && length(input.ref_plane_dist) ==1
            index= find(strcmp(Hidden.Names,'Fd'));
            Hidden.Val{index,iPat} = input.ref_plane_dist; %#ok<*FNDSB>
        else
            if ~isempty(HiddenIndex)  && any(strcmp(Hidden.Names,'Fd'))
                error('There is an error with the input ''ref_plane_dist'' ')
            end
        end
        if strcmp(curNam,'DAv')
            fixedInp{i,2} = [input.distance_av_space];
        elseif any(strcmp(Hidden.Names,'DAv')) && length(input.distance_av_space) ==1
            index= find(strcmp(Hidden.Names,'DAv'));
            Hidden.Val{index,iPat} = input.distance_av_space; %#ok<*FNDSB>
        else
            if ~isempty(HiddenIndex)  && any(strcmp(Hidden.Names,'DAv'))
                error('There is an error with the input ''distance_av_space'' ')
            end
        end
        if strcmp(curNam,'SlAv')
            fixedInp{i,2} = [input.points_av_slice];
        elseif any(strcmp(Hidden.Names,'SlAv')) && length(input.points_av_slice) ==1
            index= find(strcmp(Hidden.Names,'SlAv'));
            Hidden.Val{index,iPat} = input.points_av_slice; %#ok<*FNDSB>
        else
            if ~isempty(HiddenIndex)  && any(strcmp(Hidden.Names,'SlAv'))
                error('There is an error with the input ''points_av_slice'' ')
            end
        end
        if strcmp(curNam,'SmpR')
            fixedInp{i,2} = [input.sample_rate];
        elseif any(strcmp(Hidden.Names,'SmpR')) && length(input.sample_rate) ==1
            index= find(strcmp(Hidden.Names,'SmpR'));
            Hidden.Val{index,iPat} = input.sample_rate; %#ok<*FNDSB>
        else
            if ~isempty(HiddenIndex)  && any(strcmp(Hidden.Names,'SmpR'))
                error('There is an error with the input ''sample_rate'' ')
            end
        end
        if strcmp(curNam,'Ns')
            fixedInp{i,2} = [input.noise_U];
        elseif any(strcmp(Hidden.Names,'Ns')) && length(input.noise_U) ==1
            index= find(strcmp(Hidden.Names,'Ns'));
            Hidden.Val{index,iPat} = input.noise_U; %#ok<*FNDSB>
        else
            if ~isempty(HiddenIndex) && any(strcmp(Hidden.Names,'Ns'))
                error('There is an error with the input ''noise_U'' ')
            end
        end
    end
    % check if something is commented out from fixed inputs and add the value
    UserFixed{iPat} =  fixedInp; %#ok<*AGROW>
end

%Create permutation names for all cases by adding the inputs together
%Create names for lidar and processing
for iPat = 1:numel(input.PatternNames)
    totalMat = [freeInp; UserFixed{iPat}];
    [ArrayVar,~] = GetVariationArray(totalMat); % Sorting them and finding all the unique combinations
    TotalArrayVar{iPat} = ArrayVar;
    for iVary = 1:size(ArrayVar,1)
        if ~isempty(HiddenIndex)
            TotalArrayVarVal{iVary,iPat} = [num2cell(ArrayVar{iVary}) Hidden.Val(:,iPat)'];
        else
            TotalArrayVarVal{iVary,iPat} = [num2cell(ArrayVar{iVary})];
        end
    end
end

%retrieve names of original windfields
totalMatWF = [freeInp];
[ArrayVarWF,~] = GetVariationArray(totalMatWF);

%Create the total array for processing:
FinalArrayVar = {};
for i = 1:size( TotalArrayVar,2)
    FinalArrayVar = [FinalArrayVar ; TotalArrayVar{i}];  
end

% Create the unique names concatanating the variables and assign the values
for iAr = 1:size(FinalArrayVar,1)
    CurArray = FinalArrayVar{iAr};
    CurNam = nameBase;
    for  iVar = 1:length(CurArray)  % in pattern we want tochange the numerical value to string
        if strcmp(totalMat{iVar},'Pat')
            CurNam = [CurNam '_'  input.PatternNames{CurArray(iVar)} ];
        else
            if mod(CurArray(iVar),1) == 0  % check if the number is integer and change the notation accordingly
                CurNam = [CurNam '_' totalMat{iVar} num2str(CurArray(iVar),'%02.f') ];                
            else
                CurNam = [CurNam '_' totalMat{iVar} num2str(CurArray(iVar),'%02.1f') ];
            end
            CurNam = strrep(CurNam,'.','d');
        end
    end
    perm_cell.OutNames{iAr,1} = CurNam; 
end
TotalArrayVarVal = reshape(TotalArrayVarVal,numel(TotalArrayVarVal),1);
for i = 1:size(FinalArrayVar,1)
    perm_cell.values(i,:) = TotalArrayVarVal(i,:); 
end

if ~isempty(HiddenIndex)
    perm_cell.variables{:,1} = [totalMat(:,1); Hidden.Names];
else
    perm_cell.variables{:,1} = totalMat(:,1);    
end

%Create the total array for creating original WF names:
for iAroWF = 1:length(ArrayVarWF)
    CurArrayOWF = ArrayVarWF{iAroWF};
    CurNamWF = nameBase;
    for  iVar = 1:length(CurArrayOWF)  % in pattern we want to change the numerical value to string
        CurNamWF = [CurNamWF '_' totalMatWF{iVar} num2str(CurArrayOWF(iVar),'%02.f') ];
        CurNamWF = strrep(CurNamWF,'.','d');
    end
    perm_cell.namesOWF{iAroWF,1} = CurNamWF;
end

