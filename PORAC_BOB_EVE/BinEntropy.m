% calculate binary entropy for a probability distribution [x,1-x]
function h=BinEntropy(x)
    if ( x ~= 1 && x ~= 0 )
        h= - x * log2(x) - ( 1 - x ) * log2( 1 - x );
    else
        h = 0;
    end

