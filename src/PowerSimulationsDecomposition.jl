module PowerSimulationsDecomposition

export MultiRegionProblem
export MultiProblemTemplate

export SplitAreaPTDFPowerModel

export Relief

import PowerSimulations
import PowerNetworkMatrices
import PowerSystems
import InfrastructureSystems
import InfrastructureSystems: @assert_op
import JuMP
import Dates
import MPI
import MathOptInterface
import DataStructures: OrderedDict, SortedDict

const PSI = PowerSimulations
const PSY = PowerSystems
const PNM = PowerNetworkMatrices
const IS = InfrastructureSystems
const ISOPT = InfrastructureSystems.Optimization
const ISSIM = InfrastructureSystems.Simulation
const PM = PSI.PM
const MOI = MathOptInterface

using DocStringExtensions
@template (FUNCTIONS, METHODS) = """
                                 $(TYPEDSIGNATURES)
                                 $(DOCSTRING)
                                 """

include("core/algorithms.jl")
include("core/definitions.jl")
include("core/formulations.jl")
include("core/mpi_info.jl")
include("core/parameters.jl")
include("core/auxiliary_variables.jl")

include("multiproblem_template.jl")
include("multi_optimization_container.jl")
include("algorithms/sequential_algorithm.jl")
include("algorithms/mpi_parallel_algorithm.jl")
include("problems/multi_region_problem.jl")
include("models/network_models.jl")
include("print.jl")

end
