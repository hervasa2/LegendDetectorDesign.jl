# This file is a part of LegendDetectorDesign.jl, licensed under the MIT License (MIT).

struct BouleGeometry{N, T} <: AbstractBouleGeometry{N, T}
    z::SVector{N, T}
    radius::SVector{N, T}
    spline::Interpolations.Extrapolation{T}
end

function BouleGeometry(::Type{T};
        z::Union{Number, Vector{<:Number}},
        radius::Union{Number, Vector{<:Number}}
    ) where {T}

    @assert length(z) == length(radius) "Vectors must be the same length"
    

    if length(z) == 1 z, radius = [0*unit(z), z], [radius, radius] end
    
    N = length(z)

    idx = sortperm(z)
    d = SVector{N, T}(to_internal_length_units.(z[idx]))
    r = SVector{N, T}(to_internal_length_units.(radius[idx]))

    new_d, new_r = resample_with_min_step(d, r)
    spline = cubic_spline_interpolation(new_d, new_r)
    BouleGeometry{N, T}(d, r, spline)
end

get_boule_radius(geo::BouleGeometry, z::Number) = geo.spline(to_internal_length_units(z))

function get_physical_volume(geo::BouleGeometry{N,T})::T where {N,T}
    n_slices = 1000000
    z = range(geo.z[1], geo.z[end], n_slices)
    step = z[2] - z[1]
    r = geo.spline(z)
    ustrip(u"cm^3", sum(r.^2 * π * step)*internal_length_unit^3)
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

function get_unicode_rep(::BouleGeometry; cut::Bool = false)
    if cut 
        " ╭──╮╭──────────╮╭───╮  ", 
        "/   ││          ││    \\ ",
        "│   ││          ││     )",
        "\\   ││          ││    / ",
        " ╰──╯╰──────────╯╰───╯  "
    else
        " ╭───────────────────╮  ", 
        "/                     \\ ",
        "│                      )",
        "\\                     / ",
        " ╰───────────────────╯  "
    end
end

function print(io::IO, geo::BouleGeometry)
    g1,g2,g3,g4,g5 = get_unicode_rep(geo)
    r, l = Int.(round.((maximum(geo.radius), geo.z[end]-geo.z[1])))
    println(io, "$g1  ╮    $(typeof(geo))")
    println(io, "$g2$(lpad(string(r), 3, ' '))ₘₘ  ╰─Radius: Maximum, Minimum")
    println(io, "$g3  ╯      ╰─$(maximum(geo.radius)), $(minimum(geo.radius)) mm")
    println(io, "$g4       ╰─Length: $(geo.z[end]-geo.z[1]) mm")
    println(io, "$g5       ╰─Mass: $(Int(round(get_physical_volume(geo)*ge_76_density))) g ")
    println(io, "╰────────$(lpad(string(l), 2, ' '))ₘₘ─────────╯ ")
end