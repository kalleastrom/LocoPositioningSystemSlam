function [rsols,ids] = find_realsol(sols)

ids = find(sum(imag(sols)) == 0);
rsols = sols(:,ids);
