stationary_test = zeros(4,12,10,'single');

for i=1:10
    for j=1:12
        stationary_test(:,j,i) = isstationary(preprocessed_data{i}(:,j),0.9);
    end
end

indices = find(stationary_test==1);