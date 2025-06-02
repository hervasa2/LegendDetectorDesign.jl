# This file is a part of LegendDetectorDesign.jl, licensed under the MIT License (MIT).

function get_default_xtal_meta(det::DetectorDesign)
    PropDict(
        :impurity_curve => PropDict(
                :model => "constant_boule",
                :parameters => PropDict("value" => 0,)
                ),
        :slices => PropDict(Symbol(det.name[end]) => PropDict("detector_offset_in_mm" => 0))
    )
end

fit_function(model::Symbol) = fit_function(Val(model))
fit_parameter_units(model::Symbol) =  fit_parameter_units(Val(model))

@. linear_boule(z,p) = p[1] + p[2]*z
fit_function(::Val{:linear_boule}) = linear_boule

function fit_parameter_units(::Val{:linear_boule})
    [
        internal_impurity_quantity,
        internal_impurity_quantity / internal_length_unit
    ]
end

@. parabolic_boule(z,p) = p[1] + p[2]*z + p[3]*z^2
fit_function(::Val{:parabolic_boule}) = parabolic_boule

function fit_parameter_units(::Val{:parabolic_boule})
    [
        internal_impurity_quantity,
        internal_impurity_quantity / internal_length_unit,
        internal_impurity_quantity / internal_length_unit^2
    ]
end

@. linear_exponential_boule(z,p) = p[1] + p[2]*z + p[3]*exp((z-p[4])/p[5])
fit_function(::Val{:linear_exponential_boule}) = linear_exponential_boule

function fit_parameter_units(::Val{:linear_exponential_boule})
    [
        internal_impurity_quantity,
        internal_impurity_quantity / internal_length_unit,
        internal_impurity_quantity,
        internal_length_unit,
        internal_length_unit
    ]
end

@. parabolic_exponential_boule(z,p) = p[1] + p[2]*z + p[3]*z^2 + p[4]*exp((z-p[5])/p[6])
fit_function(::Val{:parabolic_exponential_boule}) = parabolic_exponential_boule

function fit_parameter_units(::Val{:parabolic_exponential_boule})
    [
        internal_impurity_quantity,
        internal_impurity_quantity / internal_length_unit,
        internal_impurity_quantity / internal_length_unit^2,
        internal_impurity_quantity,
        internal_length_unit,
        internal_length_unit
    ]
end