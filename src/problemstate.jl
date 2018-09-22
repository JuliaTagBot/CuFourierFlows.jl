function CuProblem(g, v, p, eq, ts)
  st = CuState(cxeltype(eq.LC), size(eq.LC), ts.dt)
  Problem(g, v, p, eq, ts, st)
end

mutable struct CuState{T,dim} <: AbstractState
  t::T
  step::Int
  dt::T
  sol::CuArray{T,dim}
end

CuState(T::DataType, sz::Tuple, dt) = CuState(T(0), 0, T(dt), cu(zeros(T, sz)))

mutable struct CuDualState{T,dimc,dimr} <: AbstractState
  t::T
  step::Int
  dt::T
  solc::CuArray{T,dimc}
  solr::CuArray{T,dimr}
end

function CuDualState(T::DataType, sizec, sizer, dt)
  solc = CuArray(zeros(T, sizec))
  solr = CuArray(zeros(T, sizer))
  CuDualState(T(0), 0, T(dt), solc, solr)
end

# Equation Composite Type
"""
This type defines the linear implicit and explicit components of an equation.
The linear implicit part of an equation is defined by an array of coefficients
which multiply the solution. The explicit part of an equation is calculated
by a function that may define linear and nonlinear parts.
"""
struct CuEquation{T,dim} <: AbstractEquation
  LC::CuArray{T,dim} # Coeffs of the eqn's implicit linear part
  calcN!::Function # Function that calcs linear & nonlinear parts
end

struct DualCuEquation{T,dimc,dimr} <: AbstractEquation
  LCc::CuArray{Complex{T},dimc}
  LCr::CuArray{Complex{T},dimr}
  calcN!::Function
end

CuEquation(eq::AbstractEquation) = CuEquation(cu(eq.LC), eq.calcN!)
