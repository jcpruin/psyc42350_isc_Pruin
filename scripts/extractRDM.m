%%extract RDM funtion
function flatRDM = extractRDM(inputmatrix)

if size(inputmatrix,1)==size(inputmatrix,2)
    flatRDM =inputmatrix(tril(true(size(inputmatrix)),-1));
    end
end


