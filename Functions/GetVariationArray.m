%% Header
%
% function creating the total permutations ofnames and relevant variables 
% for all the files requested based on the input file
%
% V.Pettas/F.Costa 
% University of Stuttgart, Stuttgart Wind Energy (SWE) 2019


function [ArrayVar,VariationPerVariable] = GetVariationArray(totalMat)

ArrayVar = {};
%     totalMat= [freeInp; UserFixed{iPat}];
Variables = size(totalMat,1);
for i = 1:Variables
    VariationPerVariable(i) = length(totalMat{i, 2});
end

veccur = 1; % initialize vector of variations of current variable

% the logic is to combine all variations of each variable in an a array with
% the next one. Then keep this matrix and combine it with the next one, and
% so on, until the one before the last.
for i = 1:Variables-1
    MatrixVar = []; % bad programming but we have to make sure that the dimensions are kept to Variations(i)*Variations(i-1)
    veccur = veccur*VariationPerVariable(i) ; % size of vector of variations of current variable
    if i == 1              % if is the first one there is no array from previous to take
        for ii = 1: veccur % for each variation take all combinations with variations of next variable
            for iii = 1:size(totalMat{i+1, 2},2) % loop over variations of next variables
                MatrixVar{iii,ii} = [totalMat{i,2}(ii) totalMat{i+1, 2}(iii) ]; % creating matrix of variations
            end
        end
        ArrayVar = reshape(MatrixVar,[numel(MatrixVar),1]);% make the matrix a vextor again to continue to the next multiplication
    else
        for ii = 1: veccur % from 2nd step there is a previous array already created
            for iii = 1:size(totalMat{i+1, 2},2)
                MatrixVar{iii,ii} = [ArrayVar{ii} totalMat{i+1, 2}(iii) ]; %#ok<*AGROW> % creating matrix of variations
            end
        end
        ArrayVar = reshape(MatrixVar,[numel(MatrixVar),1])   ; % make the matrix a vextor again to continue to the next multiplication
    end
end
