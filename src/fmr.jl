function fct2()
    println("undef")
end
function fct1()
    fct2()
end
function fct2()
    println("fct2")
end

fct1()

