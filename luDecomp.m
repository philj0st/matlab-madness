% LU Decomposition of regular square matrices
% @author Philipp Jost
% returns Permutationmatrix, Lower Triangle and Upper Triangle
function [L,U,P] = luDecomp( A )
[m,n] = size(A);
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
    % find absolute max of the column i within rows i to n
    [pivot, pivotIndex] = max(abs(U(i:end,i)));
    %pivotIndex must be respective to whole matrix
    pivotIndex = pivotIndex + i -1;
    % if pivot is not in top row
    if pivotIndex ~= i
        % swap top with pivot row
        U([i pivotIndex],:) = U([pivotIndex i],:);
        % do the same for permutation matrix
        P([i pivotIndex],:) = P([pivotIndex i],:);
        % and the factors in L
        L([i pivotIndex],:) = L([pivotIndex i],:);
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
% add identity to L so it will be able to be multiplied with
L = L + eye(n);
end

