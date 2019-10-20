struct VlasovPoisson

    interp_x1 :: InterpolatorSpline1D
    interp_x2 :: InterpolatorSpline1D
    interp_x3 :: InterpolatorSpline1D
    interp_x4 :: InterpolatorSpline1D

    function VlasovPoisson( degree, xmin, xmax, ncells )

        inter_x1 = InterpolatorSpline1D( degree, xmin[1], xmax[1], ncells[1])
        inter_x2 = InterpolatorSpline1D( degree, xmin[2], xmax[2], ncells[2])
        inter_x3 = InterpolatorSpline1D( degree, xmin[3], xmax[3], ncells[3])
        inter_x4 = InterpolatorSpline1D( degree, xmin[4], xmax[4], ncells[4])

        ex  = zeros(Float64, (ncells[1], ncells[2]))
        ey  = zeros(Float64, (ncells[1], ncells[2]))
        rho = zeros(Float64, (ncells[1], ncells[2]))

    end

end

function advection_x1(interp_x1, f, dt)

    for l=1,ncells[4], k=1,ncells[3]
        alpha = (xmin[3] +(k-1)*delta[3])*dt
        for j=1,ncells[2]
            interpolate!(interp_x1, f[:,j,k,l], alpha)
        end
    end

end

function advection_x2(interp_x2, f, dt)

    for l=1,ncells[4]
        alpha = (xmin[4] +(l-1)*delta[4])*dt
        for k=1,ncells[3], i=1,ncells[1]
            interpolate!(interp_x2, f[i,:,k,l], alpha)
        end
    end

end 

function advection_x3(interp_x3, ft, dt)

  for l=1,ncells[4], j=1,ncells[2], i=1,ncells[1]
     alpha = ex[i,j] * dt
     interpolate!(interp_x3, ft[:,l,i,j], alpha)
  end

end 

function advection_x4(interp_x4, ft, dt)

  for k=1:ncells[3], j=1:ncells[2], i=1:ncells[1]
      alpha = ey[i,j] * dt
      interpolate!(interp_x4, ft[k,:,i,j], alpha)
  end

end
