function [coefs] = decomposeEigenFunction(eigenFunctions,data)

PSIMAT = eigenFunctions;
BCOEF = (PSIMAT.'*PSIMAT)\(PSIMAT.'*data);
coefs = BCOEF;

end
