function dataStruct = compressStruct(dataIn)
empty_elems = arrayfun(@(s) all(structfun(@isempty,s)),dataIn);
dataStruct = dataIn(~empty_elems);