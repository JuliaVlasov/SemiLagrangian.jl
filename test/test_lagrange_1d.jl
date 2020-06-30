@testset " Lagrange 1D " begin

  # import SemiLagrangian:lagrange_interpolation_1d_disp_fixed_no_bc
  # import SemiLagrangian:lagrange_interpolation_1d_disp_fixed_periodic
  # import SemiLagrangian:lagrange_interpolation_1d_disp_fixed_periodicl
  # import SemiLagrangian:lagrange_interpolation_1d_disp_centered_periodicl

  n = 100
  alpha = 0.2

  fi = zeros(n+1)
  xi = zeros(n+1)
  fp = zeros(n+1)
  xp = zeros(n+1)
  
  f(x, n) = cos(2π * x / n)

  for i = 1:n+1
    xi[i] = i - 1
    fi[i] = f(xi[i], n)
    xp[i] = xi[i] + alpha
  end

   
  lagrange_interpolation_1d_disp_fixed_no_bc(view(fi,1:n), view(fp,1:n), alpha, 5)

  @test maximum(abs.(fp[1:n] .- f.(xp[1:n],n))) < 1e-8

  lagrange_interpolation_1d_disp_fixed_periodic(view(fi,1:n), view(fp,1:n), alpha, 3)

  @test maximum(abs.(fp[1:n] .- f.(xp[1:n],n))) < 8e-6

  lagrange_interpolation_1d_disp_fixed_periodicl(fi, fp, alpha, 5)

  @test maximum(abs.(fp .- f.(xp,n))) < 7e-9

  lagrange_interpolation_1d_disp_centered_periodicl(fi, fp, alpha, 4)

  @test maximum(abs.(fp .- f.(xp,n))) < 3e-7

  lagrange_interpolation_1d_disp_centered_periodicl(fi, fp, alpha, 6)
       
  @test maximum(abs.(fp .- f.(xp,n))) < 2e-10

end
