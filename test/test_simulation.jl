


struct SLModel1D1V

    fn :: Array{Float64, 2}
    ft :: Array{Float64, 2}
    nt :: Int
    dt :: Float64
    interp :: InterpolationType
    splitting :: SplittingType

    function SLModel1D1V( df, interp, tspan, nt )

        fn = zero(df)
        ft = collect(transpose(fn))
        dt = tf / nt

        new( fn, ft, nt, dt, interp, splitting )

    end


end
