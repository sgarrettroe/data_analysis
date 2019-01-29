function out = cumkronSparse(input_cell)
% calculate the cumulative tensor (kronecker) product of a cell array of
% operators

temp1 = input_cell(2:end);%round braces are important so temp1 stays a cell
temp2 = input_cell{1}; %curly braces here make temp2 a matrix
for ii = 1:length(temp1)
    temp2 = kron(temp2,temp1{ii});%the curly braces on temp1 take out the matrix
end

out = temp2;
end
