


function fct(a,b,c)

    @show a,b,c

    function fct1(v)
        a,b,c  = (a,b,c) .+ v
    end


    fct1(4)

    @show a, b, c
end



fct(1,2,3)
