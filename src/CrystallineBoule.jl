# This file is a part of LegendDetectorDesign.jl, licensed under the MIT License (MIT).

abstract type AbstractCrystallineBoule{T <: SSDFloat} end

mutable struct CrystallineBoule{T} <: AbstractCrystallineBoule{T}
    name::AbstractString
    order::AbstractString
    impurity_model::Function
    geometry::Union{BouleGeometry{<:Any,T}, Missing}
    impurity_resistivity::Union{Vector{T}, Missing}
    z_resistivity::Union{Vector{T}, Missing}
    impurity_hall::Union{Vector{T}, Missing}
    z_hall::Union{Vector{T}, Missing}
    impurity_model_parameters::Union{Vector{T}, Missing}
end

function CrystallineBoule(::Type{T};
        name::AbstractString,
        order::AbstractString,
        impurity_model::Symbol,
        geometry::Union{BouleGeometry{<:Any, T}, Missing} = missing,
        impurity_resistivity::Union{Vector{<:Number}, Missing} = missing,
        z_resistivity::Union{Vector{<:Number}, Missing} = missing,
        impurity_hall::Union{Vector{<:Number}, Missing} = missing,
        z_hall::Union{Vector{<:Number}, Missing} = missing,
        impurity_model_parameters::Union{Vector{<:Number}, Missing} = missing
    ) where {T <: SSDFloat}
    CrystallineBoule{T}(name, order, fit_function(impurity_model), geometry, to_internal_units.(impurity_resistivity), to_internal_units.(z_resistivity), to_internal_units.(impurity_hall), to_internal_units.(z_hall), to_internal_units.(impurity_model_parameters))
end

get_unitful_property(boule::CrystallineBoule, prop::Symbol) = get_unitful_property(boule, Val(prop))

get_unitful_property(boule::CrystallineBoule, ::Val{:impurity_model_parameters}) = boule.impurity_model_parameters .* fit_parameter_units(nameof(boule.impurity_model))

get_impurity_density(boule::CrystallineBoule, z::Number) = boule.impurity_model(to_internal_units(z), boule.impurity_model_parameters) * internal_impurity_quantity

#need to be consistent with get_ functions (do they always return unitfill values?)
get_impurity_density(boule, z::AbstractArray) = broadcast(z -> get_impurity_density(boule, z), z)