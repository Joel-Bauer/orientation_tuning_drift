function [r_mean, numb_of_non_nan_elements] = cocoe_mean(r_vector) % correlation coefficient mean calculator, with graphic representation of transformation
r_vector(isnan(r_vector)) = []; % remove nan elements
numb_of_non_nan_elements = length(r_vector);

if numb_of_non_nan_elements == 0
    r_mean = nan;
elseif numb_of_non_nan_elements == 1
    r_mean = r_vector;
elseif numb_of_non_nan_elements > 1
    z_vector = arrayfun(@(r)atanh(r), r_vector);
    z_mean = mean(z_vector);
    r_mean = tanh(z_mean);
end

end