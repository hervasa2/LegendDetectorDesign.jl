# This file is a part of LegendDetectorDesign.jl, licensed under the MIT License (MIT).

"""
    LegendDetectorDesign

Detector design tools for the LEGEND experiment.
"""
module LegendDetectorDesign

using LinearAlgebra
using PropDicts
using Unitful
using LegendDataManagement
using SolidStateDetectors

import SolidStateDetectors: 
            SSDFloat, AbstractImpurityDensity, Simulation, AbstractCoordinatePoint
import Base: show, print, println

include("Units.jl")
include("Geometry/Geometry.jl")
include("DetectorDesign.jl")
include("ImpurityDensities.jl")
include("Characterize.jl")

export DetectorDesign, InvertedCoaxDesign, InvertedCoaxGeometry, SegregationImpurityDensity, ValidGeometry, InvalidGeometry, characterize!

end # module
