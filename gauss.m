% Philipp Jost
% Linalg S2

function x = gauss( A, b )

% A is a n*n matrix and b a n*1 vector

    % get length of vector b
	n = length(b);
	
	% If there is no unique solution then "false" is returned
	x = false;
	
	% We use the Gauss algorithm to transform B into 
	% row echelon form
	B = [A, b];	
	
	% Transform B into upper triangular form
    % for n rows there will hopefully be n pivots
	for k = 1:n
        
		% Check for Zero column
		col = B(:,1);
        zeroCol = zeros(length(col),1);
        if col == zeroCol
            break % terminate execution of the for loop k = 1:n
        end
		
		% Find pivot element
        % [Y,I] = max(X) returns the indices of the maximum values in vector I.
        % searching the untematrix from k to n rows in the n-th column
		[pivot, pivotInd] = max(abs(B(k:n, k))); %in which case is pivot as return value of max ddifferent from manual retreiving with index
		% Index of pivot element
		pivotInd = pivotInd + (k-1); %add k-1 because were only searching max from k:n
 		pivot = B(pivotInd, k); %i think this is unnecessary => definitely not unessecary .. but why?
        
        % not quite sure if this is needed but it will break execution
        % if no pivot can be choosen which is not zero
        if pivot == 0
            x = false;
            return
        end
        
		% Swap pivot row
		temp = B(k, :);
		B(k, :) = B(pivotInd,:);
		B(pivotInd, :) = temp;
		
		% Normalize row
        % pivots will become 1
		B(k, :) = B(k, :) / pivot;
		
		% Eliminate rest of column
		for l = k+1:n
            % subtract x times the row with the pivot from the rows
            % underneath the pivot, where x is the first element in the row
            % to be eliminated
			B(l,:) = B(l,:) - B(l,k)*B(k,:);
		end
    end
        
	% Back substitution to solve for x
	% Transform B into the form [I, x]
    
    
    % for each row except for the last (since this is already 1x = _)
	for k = (n-1):-1:1
        % subtract x times the row that eliminated everything except 1x
        % to bring B in reduced row echeleon form        
		for l = k+1:1:n %!! changed the loop
            B(k,:) = B(k,:) - B(k,l) * B(l,:);
		end
	end	

	% x is now the last column of B
	x = B(:, end);
end