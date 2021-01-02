% calculates lower bound on key-rate
function x=keyRate(SB,SE)
    x=-1*BinEntropy(SB)-log2(SE)