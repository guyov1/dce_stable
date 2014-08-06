function D = distfun(X, C, dist, iter,rep, reps)
%DISTFUN Calculate point to cluster centroid distances.
[n,p] = size(X);
D = zeros(n,size(C,1));
nclusts = size(C,1);

switch dist
    case 'sqeuclidean'
        for i = 1:nclusts
            D(:,i) = (X(:,1) - C(i,1)).^2;
            for j = 2:p
                D(:,i) = D(:,i) + (X(:,j) - C(i,j)).^2;
            end
            % D(:,i) = sum((X - C(repmat(i,n,1),:)).^2, 2);
        end
    case 'cityblock'
        for i = 1:nclusts
            D(:,i) = abs(X(:,1) - C(i,1));
            for j = 2:p
                D(:,i) = D(:,i) + abs(X(:,j) - C(i,j));
            end
            % D(:,i) = sum(abs(X - C(repmat(i,n,1),:)), 2);
        end
    case {'cosine','correlation'}
        % The points are normalized, centroids are not, so normalize them
        normC = sqrt(sum(C.^2, 2));
        if any(normC < eps(class(normC))) % small relative to unit-length data points
            if reps==1
                error(message('stats:kmeans:ZeroCentroid', iter));
            else
                error(message('stats:kmeans:ZeroCentroidRep', iter, rep));
            end
            
        end
        
        for i = 1:nclusts
            D(:,i) = max(1 - X * (C(i,:)./normC(i))', 0);
        end
    case 'hamming'
        for i = 1:nclusts
            D(:,i) = abs(X(:,1) - C(i,1));
            for j = 2:p
                D(:,i) = D(:,i) + abs(X(:,j) - C(i,j));
            end
            D(:,i) = D(:,i) / p;
            % D(:,i) = sum(abs(X - C(repmat(i,n,1),:)), 2) / p;
        end
end
end % function