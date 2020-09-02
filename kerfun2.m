function distmat = kerfun2(data1, data2, sigmaList, dim)
if isequal(data1, data2)
    flag = 1; m1 = length(data1); m2 = m1;
else
    flag = 0; m1 = length(data1); m2 = length(data2);
end
% to generate all pairings of rows of fmat1 rows of fmat2
distmat = zeros(m1, m2, length(sigmaList));
[y, x] = meshgrid(1:m1, 1:m2);
idx = [y(:), x(:)];
if flag == 1
    idx(idx(:,1) < idx(:,2), :) = [];
end
for k = 1:size(idx, 1)
    disp(k);
    % pick vector from each feature matrix
    iy = idx(k, 1);     ix = idx(k, 2);
    diag1 = data1(iy).pers_data;
    diag2 = data2(ix).pers_data;
    % only pick points of given dimension
    diag1(diag1(:,1) ~= dim, :) = [];
    diag2(diag2(:,1) ~= dim, :) = [];
    % remove the dimension column - not needed anymore
    diag1(:,1) = [];    diag2(:,1) = [];
    % compute the inner product between two diagrams
    distmat(iy, ix, :) = ksigma2(diag1, diag2, sigmaList);
    if flag == 1
        distmat(ix,iy,:) = distmat(iy,ix,:);
    end
end

end
