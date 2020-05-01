@testset " Lagrange 1D " begin

  import SemiLagrangian:lagrange_interpolation_1d_disp_fixed_no_bc

  num_points = 100
  alpha = 0.2
  
  xi = zeros(num_points+1)
  fi = zeros(num_points+1)

  # x values with offset
  xp = zeros(num_points+1)

  # data initialization

  f(x, num_points) = cos(2π * x / num_points)

  xmin = 0.0
  xmax = num_points-1.0
  l = xmax - xmin
  for i = 1:num_points+1
    xi[i] = i - 1
    fi[i] = f(xi[i], num_points)
    xp[i] = xi[i] + alpha
  end

  function test_interpolation(num_points, fi, alpha, xp, order, type)

    fp = zeros(num_points)

    if type > 1
       num_cells = num_points-1
    else
       num_cells = num_points
    end

    if type == 0 
       @info "Test fixed_no_bc with order $order "
       lagrange_interpolation_1d_disp_fixed_no_bc(fi, fp, alpha, order)
    #elseif type == 1 # periodic
    #   println("Test fixed_periodic with order $order")
    #   lagrange_interpolation_1d_disp_fixed_periodic(fi, fp, alpha, order)
    #elseif type == 2 # periodic with last value
    #   println("Test fixed_periodic_last with order $order")
    #   lagrange_interpolation_1d_disp_fixed_periodicl(fi, fp, alpha, order)
    #elseif type == 3 # periodic centered
    #   println("Test centered_periodic_last with order $order")
    #   lagrange_interpolation_1d_disp_centered_periodicl(fi, fp, alpha, order)
    #else
    #   @error "Interpolation type not implemented"
    end
       

    diff = maximum(abs.(f.(xp, num_points) .- fp))

    println("error = $diff ")

    return diff

  end 


  @test test_interpolation(num_points, fi[1:num_points], alpha, xp[1:num_points], 5, 0) < 1e-8
# @test test_interpolation(num_points, fi[1:num_points], alpha, xp[1:num_points], 3, 1) < 8e-6
# @test test_interpolation(num_points+1, fi, alpha, xp, 5, 2) < 7e-9
# @test test_interpolation(num_points+1, fi, alpha, xp, 4, 3) < 3e-7
# @test test_interpolation(num_points+1, fi, alpha, xp, 6, 3) < 2e-10

end
