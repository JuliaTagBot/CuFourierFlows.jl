struct CuTwoDGrid{T} <: AbstractTwoDGrid
  nx::Int
  ny::Int
  nk::Int
  nl::Int
  nkr::Int
  Lx::T
  Ly::T
  dx::T
  dy::T
  x::CuArray{T,2}
  y::CuArray{T,2}
  X::CuArray{T,2}
  Y::CuArray{T,2}
  k::CuArray{T,2}
  l::CuArray{T,2}
  kr::CuArray{T,2}
  Ksq::CuArray{T,2}      # K^2 + L^2
  invKsq::CuArray{T,2}   # 1/KKsq, invKKsq[1, 1]=0
  Krsq::CuArray{T,2}     # Kr^2 + Lr^2
  invKrsq::CuArray{T,2}  # 1/KKrsq, invKKrsq[1, 1]=0

  fftplan::CuArrays.FFT.cCuFFTPlan{Complex{T},-1,false,2}
  rfftplan::CuArrays.FFT.rCuFFTPlan{T,-1,false,2}

  # Range objects that access the aliased part of the wavenumber range
  ialias::UnitRange{Int}
  iralias::UnitRange{Int}
  jalias::UnitRange{Int}
end

"""
    CuTwoDGrid(nx, Lx, ny=nx, Ly=Lx;  x0=-Lx/2, y0=-Ly/2)

Construct a CuTwoDGrid object. The two-dimensional domain has size `(Lx, Ly)`,
resolution `(nx, ny)` and bottom left corner at `(x0, y0)`. The type of `CuTwoDGrid` is
inferred from the type of `Lx`.
"""
function CuTwoDGrid(nx, Lx, ny=nx, Ly=Lx; x0=-0.5*Lx, y0=-0.5*Ly)
  T = typeof(Lx)
  dx = Lx/nx
  dy = Ly/ny
  nk = nx
  nl = ny
  nkr = Int(nx/2+1)

  # Physical grid
  x = Array{T}(reshape(linspace(x0, x0+Lx-dx, nx), (nx, 1)))
  y = Array{T}(reshape(linspace(y0, y0+Ly-dy, ny), (1, ny)))
  X = [ x[i] for i = 1:nx, j = 1:ny]
  Y = [ y[j] for i = 1:nx, j = 1:ny]

  # Wavenubmer grid
  i1 = 0:Int(nx/2)
  i2 = Int(-nx/2+1):-1
  j1 = 0:Int(ny/2)
  j2 = Int(-ny/2+1):-1

  k  = reshape(2π/Lx*cat(1, i1, i2), (nk, 1))
  kr = reshape(2π/Lx*cat(1, i1), (nkr, 1))
  l  = reshape(2π/Ly*cat(1, j1, j2), (1, nl))

  Ksq  = @. k^2 + l^2
  invKsq = 1 ./ Ksq
  invKsq[1, 1] = 0

  Krsq = @. k^2 + lr^2
  invKrsq = 1 ./ Krsq
  invKrsq[1, 1] = 0

  # Convert Arrays to CuArrays
  @cuconvertarrays x y X Y k kr l Ksq invKsq Krsq invKrsq

  # FFT plans
    fftplan = plan_fft(CuArray{Complex{T},2}(nx, ny))
   rfftplan = plan_rfft(CuArray{T,2}(nx, ny))

  # Index endpoints for aliased i, j wavenumbers
  iaL, iaR = Int(floor(nk/3))+1, 2*Int(ceil(nk/3))-1
  jaL, jaR = Int(floor(nl/3))+1, 2*Int(ceil(nl/3))-1

  ialias  = iaL:iaR
  iralias = iaL:nkr
  jalias  = iaL:iaR

  CuTwoDGrid(nx, ny, nk, nl, nkr, Lx, Ly, dx, dy, x, y, X, Y,
             k, l, kr, Ksq, invKsq, Krsq, invKrsq,
             fftplan, rfftplan, ialias, iralias, jalias)
end

"""
    CuTwoDGrid(T, args...; kwargs...)

Construct a CuTwoDGrid object with type T. The other args and kwargs are identical
to the constructor for CuTwoDGrid with no explicit type.
"""
CuTwoDGrid(T::DataType, nx, Lx, args...; kwargs...) = CuTwoDGrid(nx, T(Lx), args...; kwargs...)

"""
    CuTwoDGrid(g)

Construct a CuTwoDGrid object with the same type and size as the TwoDGrid object `g`.
"""
CuTwoDGrid(g::FourierFlows.TwoDGrid) = CuTwoDGrid(g.nx, g.Lx, g.ny, g.Ly; x0=g.x[1], y0=g.y[1])

#=
"""
    makefilter(K; order=4, innerK=0.65, outerK=1)

Returns a filter acting on the non-dimensional wavenumber K that decays exponentially
for K>innerK, thus removing high-wavenumber content from a spectrum it is multiplied with.
The decay rate is determined by order and outerK determines the outer wavenumber at which
the filter is smaller than machine precision.
"""
function makefilter(K; order=4, innerK=0.65, outerK=1)
  TK = typeof(K)
  K = Array(K)
  decay = 15*log(10) / (outerK-innerK)^order # decay rate for filtering function
  filt = @. exp( -decay*(K-innerK)^order )
  filt[real.(K) .< innerK] .= 1
  TK(filt)
end

function makefilter(g::AbstractTwoDGrid; realvars=true, kwargs...)
  K = realvars ?
      @.(sqrt((g.Kr*g.dx/π)^2 + (g.Lr*g.dy/π)^2)) : @.(sqrt((g.K*g.dx/π)^2 + (g.L*g.dy/π)^2))
  makefilter(K; kwargs...)
end

function makefilter(g::AbstractOneDGrid; realvars=true, kwargs...)
  K = realvars ? g.kr*g.dx/π : @.(abs(g.k*g.dx/π))
  makefilter(K; kwargs...)
end

makefilter(g, T, sz; kwargs...) = ones(T, sz).*makefilter(g; realvars=sz[1]==g.nkr, kwargs...)
=#

makefilter(g::CuTwoDGrid, T, sz; kwargs...) = cuones(T, sz).*makefilter(g; realvars=sz[1]==g.nkr, 
                                                                                     kwargs...)
