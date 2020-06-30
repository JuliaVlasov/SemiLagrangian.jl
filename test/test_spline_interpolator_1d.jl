

@testset " Spline Interpolator 1D cos " begin

    function f( x; d=0 ) 
    
        k = 2π
    
        k^d .* cos.( 0.5*π*d .+ k*x )
    
    end

    xmin, xmax = 0.0, 1.0
    ncells = 10

    for degree = 1:9

        obj =  InterpolatorSpline1D( degree, xmin, xmax, ncells )

        if isodd(degree)
            @test obj.tau ≈ [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
        else
            @test obj.tau ≈ [0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95]
        end

        gtau = map( f, obj.tau ) 

        derivs_left = zeros(Float64, degree÷2)
        s = Int64(iseven(degree))
        for j = 1:degree÷2
            derivs_left[j] = f( xmin; d=j-s )
        end

        derivs_right = zeros(Float64, degree÷2)
        for j = 1:degree÷2
            derivs_right[j] = f( xmax; d=j-s )
        end

        compute_interpolant!( obj, gtau, derivs_left, derivs_right )

        error = 0.0
        for (i, tau) in enumerate(obj.tau)
            delta = gtau[i] - eval_value( obj.bspl, tau )
            error += abs( delta )
        end

        @test error ≈ 0.0 atol = 1e-14

        nx = 20
        error = 0.0
        for x in range(xmin, stop=xmax, length=nx) 
            y = f(x)
            delta = y - eval_value( obj.bspl, x )
            error += abs( delta)
        end

        @test error/nx < obj.bspl.step^degree

        error = 0.0
        for x in obj.tau
            y = f(x; d=1)
            delta = y - eval_deriv( obj.bspl, x )
            error += abs( delta)
        end

        @test error/nx < obj.bspl.step^(degree-1)

    end 


end

@testset " Spline Interpolator 1D polynomial " begin

    function falling_factorial( x :: Int64, n :: Int64 ) 
        c = 1
        for k = 0:n-1
            c = c * (x-k)
        end
        c
    end

    xmin, xmax, ncells = -1.0, 1.0, 23

    for degree = 1:9

        coeffs = 1.0 .- rand(degree+1)

        function f( x :: Float64; d = 0 :: Int64 )

            y = 0.0
            for i = d:degree
               c = falling_factorial( i, d ) * x^(i-d)
               y = y + coeffs[i+1] * c
            end
            y

        end

        obj =  InterpolatorSpline1D( degree, xmin, xmax, ncells )

        gtau = map( f, obj.tau ) 

        derivs_left = zeros(Float64, degree÷2)
        s = Int64(iseven(degree))
        for j = 1:degree÷2
            derivs_left[j] = f( xmin; d=j-s )
        end

        derivs_right = zeros(Float64, degree÷2)
        for j = 1:degree÷2
            derivs_right[j] = f( xmax; d=j-s )
        end

        compute_interpolant!( obj, gtau, derivs_left, derivs_right )

        error = 0.0
        for (i, tau) in enumerate(obj.tau)
            delta = gtau[i] - eval_value( obj.bspl, tau )
            error += abs( delta )
        end

        @test error ≈ 0.0 atol = 1e-14

        nx = 20
        error = 0.0
        for x in range(xmin, stop=xmax, length=nx) 
            y = f(x)
            delta = y - eval_value( obj.bspl, x )
            error += abs( delta)
        end

        @test error/nx ≈ 0.0 atol = 1e-14

        error = 0.0
        for x in obj.tau
            y = f(x; d=1)
            delta = y - eval_deriv( obj.bspl, x )
            error += abs( delta)
        end

        @test error/nx ≈ 0.0 atol = 1e-13

    end 


end
