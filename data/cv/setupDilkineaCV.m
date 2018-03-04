X = zeros(100, 10);
for i = 1:10
    c = cvpartition(y, 'kfold', 10);
    for j = 1:10
        z = test(c,j);
        X(z,i) = j;
    end
end