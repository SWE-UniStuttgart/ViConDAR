%% Header
%
% Takes as input a numerical vector array and returns the indices of the matching values 
% in a cell. Each element of the cell are the indices that have the same
% values. If there are no repeated values the output cell is empty

%Example 
% A=[15  3 3 0 4 6 5 4 8 10];
% Match_index=FindSameValuesAndAverageV(A):
%
% V.Pettas/F.Costa 
% University of Stuttgart, Stuttgart Wind Energy (SWE) 2019

function Match_index=FindSameValuesAndAverageV(A)

Av_ind={};
xx = 0;
for i = 1:length(A)   
    %Compare one by one element and indices
    ind1 = A-A(i);    
    zer_ind = find(ind1==0); % find the matching values indices
    if length(zer_ind) == 1  %assign the matching indices for each of the repeated values
    else
        xx = xx+1;
        Av_ind{xx} = zer_ind ; %#ok<*AGROW>
    end  
end

% Now we have to keep only the unique vectors 
if ~isempty(Av_ind)
    for i = 1:length(Av_ind)
        vec_test = Av_ind{i};
        for j = 1: length(Av_ind)
            if length(Av_ind{i}) == length(Av_ind{i})
                if isequal(Av_ind{j},vec_test) && i~=j && ~isempty(vec_test) % check every vector with all the others
                    Av_ind2{j} = [];
                else
                    Av_ind2{j} = Av_ind{j}; %if not keep it
                end
            end
        end
        Av_ind = Av_ind2 ; % update the  matching indices cell with the remove values for the next iteration
        clear Av_ind2
    end
    Match_index = Av_ind(~cellfun(@isempty,Av_ind));
else
    Match_index = {};  % if no matching values exist
end
