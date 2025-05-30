# This file is a part of LegendDetectorDesign.jl, licensed under the MIT License (MIT).

struct BouleGeometry{N, T} <: AbstractBouleGeometry{N, T}
    distance_from_seed_end::SVector{N, T}
    radius::SVector{N, T}
    spline::Interpolations.Extrapolation{T}
end

function BouleGeometry{T}(distance_from_seed_end::Vector{<:Number}, radius::Vector{<:Number}) where {T}
    @assert length(distance_from_seed_end) == length(radius) "Vectors must be the same length"
    N = length(distance_from_seed_end)

    d = SVector{N, T}(to_internal_length_units.(distance_from_seed_end))
    r = SVector{N, T}(to_internal_length_units.(radius))

    new_d, new_r = resample_with_min_step(d, r)
    spline = cubic_spline_interpolation(new_d, new_r)
    BouleGeometry{N, T}(d, r, spline)
end

function resample_with_min_step(x::AbstractVector{<:Real}, y::AbstractVector{<:Real})
    @assert length(x) == length(y) "x and y must have the same length"
    @assert issorted(x) "x must be sorted in increasing order"

    # Compute the minimum step
    dxs = diff(x)
    min_step = minimum(dxs)

    # Construct StepRangeLen from min step
    x_min, x_max = first(x), last(x)
    n = Int(cld(x_max - x_min, min_step)) + 1  # cld = ceil division
    new_x = StepRangeLen(x_min, min_step, n)

    # Interpolate with flat extrapolation
    itp = LinearInterpolation(x, y, extrapolation_bc=Flat())
    new_y = SVector{n}(itp.(new_x))

    new_x, new_y
end