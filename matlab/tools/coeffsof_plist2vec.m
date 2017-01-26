function [coefficient_vector] = coeffsof_plist2vec(eqs)

coefficient_vector = [];
neq  = length(eqs);

for i = 1:neq
    eqs_mon = monomials(eqs(i));
    eqs_mon = multipol(ones(1,size(eqs_mon,2)),eqs_mon);    
    coefficient_vector = [coefficient_vector;coeffsof_p2m(eqs_mon,eqs(i))'];
end

% coefficient_vector (abs(coefficient_vector) < 1e-13) = [];