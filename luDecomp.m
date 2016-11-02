% LU Decomposition of regular square matrices
% @author Philipp Jost
% returns Permutationmatrix, Lower Triangle and Upper Triangle
function [ L U P ] = luDecomp( A )
[m n] = size(A);
if m ~= n
    throw(MException('luDecomp:noSquareMatrix','please enter a square matrix'))
end

% init with zeroes for L so we can swap rows without hurting the diagonal
% elements
L = zeros(n);
P = eye(n);
U = A;

% for each column
% i-th step of the algorithm
for i = 1:n
    % find absolute max of the column i
    [pivot pivotIndex] = max(abs(U(:,i)));
    
    % if pivot is not in top row
    if pivotIndex ~= i
        % swap first with pivot row
        U = swapRows(U,i,pivotIndex);
        % do the same for permutation matrix
        P = swapRows(P,i,pivotIndex);

        % #TODO! swap lambdas in L!
        
    end
    
    % for each row underneath pivot
    for j = (i+1):n
        
        % calculate ratio which is needed to achieve zeroes underneath
        % the pivot
        try
            lambda = U(j,i)/pivot;
        catch
                throw(MException('luDecomp:divisionByZero','Pivot cannot be 0'))
        end
        
        % save lambda to elementary matrix
        L(j,i) = lambda;
        % could also just write 0
        U(j,:) = U(j,:) - lambda * U(i,:);
    end
end
    % take an Matrix and swap row i with toSwapWith
    function A = swapRows(A, i, toSwapWith)
        tempRow = A(i,:);
        A(i,:) = A(toSwapWith,:);
        A(toSwapWith,:) = tempRow;
    end
end

