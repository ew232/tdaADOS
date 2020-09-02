
function dist = ksigma2(diag1, diag2, sigmaList)
sz1 = size(diag1, 1);
sz2 = size(diag2, 1);
% create all possible pairings of points
[i, j] = meshgrid(1:sz1, 1:sz2);
diag1 = diag1(i(:), :);
diag2 = diag2(j(:), :);
% compute exponential terms in summation
dist = zeros(length(sigmaList),1);
for si = 1:length(sigmaList)
    s = sigmaList(si);
    t1 = exp(-1 * sum((diag1 - diag2).^2, 2)/(8*s));
    t2 = exp(-1 * sum((diag1 - diag2(:,[2,1])).^2, 2)/(8*s));
    % compute inner product
    dist(si) = (1/(8*pi*s)) * (sum(t1 - t2));
end
end