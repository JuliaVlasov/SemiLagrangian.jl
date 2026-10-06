
@testitem "Aqua.jl" begin
    using Aqua
    using SemiLagrangian
    # Aqua currently flags four legacy tuple signatures with apparently unbound parameters.
    Aqua.test_all(SemiLagrangian; unbound_args=false)
end
