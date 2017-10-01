%% This function performs NMF via least square gradient descent 
function [W, H] = gdcls(V, k, maxiter, lambda, options)
    myeps = 10^-9;
    if strcmp(options, 'nonneg');
        neg = 1;
    else
        neg = 0;
    end
    
    [m,n] = size(V);
    W = rand(m, k);
    H = zeros(k, n);
    for j = 1 : maxiter
        A = W' * W + lambda * eye(k);
        for i = 1 : n
            b = W' * V(:, i);
            H(:,i) = (b'/(A'))';
%             H(:,i) = b / A;
        end
        
        % Removing any negative elements
        if neg == 1
            H = H .* (H > 0);
%             H(H<0) = 0;
        end
        
        % Updating W
        W = W .* (V * H') ./ (W * (H * H') + myeps);
        
        for qq = 1 : size(W, 2)
            W(:, qq) = W(:, qq) / sum(W(:,qq));
        end
    end
end