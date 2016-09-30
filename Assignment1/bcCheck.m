%bcCheck
function boundary=bcCheck(psi,e)
boundary = psi(1)*beta(e) + psi(2);
end