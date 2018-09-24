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
  @cucreatearrays

using 
  CuArrays,
  FFTW,
  Reexport

import FourierFlows

import FourierFlows: 
  AbstractState,

  AbstractTwoDGrid,
  AbstractOneDGrid,

  AbstractForwardEulerTimeStepper,
  AbstractFilteredForwardEulerTimeStepper,
  AbstractRK4TimeStepper,
  AbstractFilteredRK4TimeStepper,

  makefilter
  
@reexport using FourierFlows

include("utils.jl")
include("domains.jl")
include("problemstate.jl")
include("timesteppers.jl")

end # module
