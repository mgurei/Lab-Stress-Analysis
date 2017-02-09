function [ new_qrs ] = transforMat( qrs_mat,matLength )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

new_qrs = zeros(matLength,1);

for i = 1:length(qrs_mat)
    new_qrs(qrs_mat(i)) = 1;
end

end

