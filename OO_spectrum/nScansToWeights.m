function weightMatrix = nScansToWeights(dataMatrix,nScans)
weightMatrix = bsxfun(@times,ones(size(dataMatrix)),reshape(nScans,1,1,length(nScans)));