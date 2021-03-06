module CuFourierFlows

export
  CuTwoDGrid,

  CuProblem,
  CuEquation,
  CuState,
  CuDualState,

  CuForwardEulerTimeStepper,
  CuFilteredForwardEulerTimeStepper,
  CuRK4TimeStepper,
  CuFilteredRK4TimeStepper,

  @cuconvertarrays,
  @createcuarrays

using 
  CuArrays,
  FFTW,
  FourierFlows

using FourierFlows: 
  AbstractState,

  AbstractTwoDGrid,
  AbstractOneDGrid,

  AbstractForwardEulerTimeStepper,
  AbstractFilteredForwardEulerTimeStepper,
  AbstractRK4TimeStepper,
  AbstractFilteredRK4TimeStepper,

  makefilter,
  cxeltype

import FourierFlows: makefilter

include("utils.jl")
include("domains.jl")
include("problemstate.jl")
include("timesteppers.jl")

end # module
