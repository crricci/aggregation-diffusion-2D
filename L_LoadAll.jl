# progress for DifferentialEquations
using Logging: global_logger
using TerminalLoggers: TerminalLogger
global_logger(TerminalLogger())

using LinearAlgebra
using DifferentialEquations
using ImageFiltering
using Parameters
using PyPlot

include("L_parameters.jl")
include("L_diff.jl")
include("L_plotLib.jl")
