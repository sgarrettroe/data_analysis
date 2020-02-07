function displayCoeffMatrix(PSI,inputspace)

c = calculateCoeffMatrix(PSI,inputspace);
printCoeffMatrix(c,inputspace);

end

% function c = calculateCoeffMatrix(PSI,inputspace)
%
% u1 = inputspace(1).u;
% u2 = inputspace(2).u;
% nstates1 = inputspace(1).nstates;
% nstates2 = inputspace(2).nstates;
%
% %c will be a matrix of coefficients for the basis states
% c = zeros([inputspace.nstates]);
% for ii = 1:nstates1
%     for jj = 1:nstates2
%         c(ii,jj) = kron(u1(ii),u2(jj))'*PSI;
%     end
% end
%
% end

function c = calculateCoeffMatrix(PSI,inputspace)

nstates = [inputspace.nstates];
ndims = length(inputspace)+1;
nrows = prod(nstates);

%SGR: something is wrong here, argh
basis = permute(reshape(eye(nrows),[nrows nstates(end:-1:1)]),[1 ndims:-1:2]);
%basis = reshape(eye(nrows),[nrows nstates]);

temp = bsxfun(@times,basis,PSI);
c = squeeze(sum(temp)); %I don't think there should be a ' but my size is wrong!!!
%c will be a matrix of coefficients for the basis states

end

function printCoeffMatrix(c,inputspace)
% plot the matrix in a humanly readable way (?)

switch  length(inputspace)
    case 2

    nstates1 = inputspace(1).nstates;
    nstates2 = inputspace(2).nstates;
    
    %col labels
    colHeadersString = '\t';
    for jj = 1:nstates2,
        colHeadersString = strcat(colHeadersString,'|:,',num2str(jj-1),'>\t\t');
    end
    colHeadersString = strcat(colHeadersString,'\n');
    
    % row labels
    rowHeadersCell = cell(nstates1,1);
    for ii = 1:nstates1,
        rowHeadersCell(ii) = {strcat('|',num2str(ii-1),',:>\t')};
    end
    
    
    % print matrix and labels
    fprintf(1,colHeadersString);
    for ii = 1:nstates1
        fprintf(1,rowHeadersCell{ii});
        fprintf(1,'%5.2f + %5.2fi\t',[real(c(ii,:));imag(c(ii,:))]);
        fprintf(1,'\n');
    end
    fprintf(1,'\n');
    
    otherwise
        %punt for now...
        disp(c)
end

end
