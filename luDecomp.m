% LU Decomposition of regular square matrices
% @author Philipp Jost
function [ L U ] = luDecomp( A )
[m n] = size(A);
if m ~= n
    throw(MException('luDecomp:noSquareMatrix','please enter a square matrix'))
end

L = eye(n);
U = A;

%for each column
for i = 1:n
    pivot = U(i,i);
    
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
        U(j,:) = U(j,:) - lambda * U(i,:);
    end
end
end

