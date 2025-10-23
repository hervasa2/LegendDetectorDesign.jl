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
impurity_density_model(model::Symbol) = impurity_density_model(Val(model))
fit_parameter_names(model::Symbol) =  fit_parameter_names(Val(model))
fit_parameter_units(model::Symbol) =  fit_parameter_units(Val(model))

@. linear_boule(z,p) = p[1] + p[2]*z
fit_function(::Val{:linear_boule}) = linear_boule
impurity_density_model(::Val{:linear_boule}) = LinBouleImpurityDensity
fit_parameter_names(::Val{:linear_boule}) = Symbol.(["a", "b"])

function fit_parameter_units(::Val{:linear_boule})
    [
        internal_impurity_quantity,
        internal_impurity_quantity / internal_length_unit
    ]
end

@. parabolic_boule(z,p) = p[1] + p[2]*z + p[3]*z^2
fit_function(::Val{:parabolic_boule}) = parabolic_boule
impurity_density_model(::Val{:parabolic_boule}) = ParBouleImpurityDensity
fit_parameter_names(::Val{:parabolic_boule}) = Symbol.(["a", "b", "c"])

function fit_parameter_units(::Val{:parabolic_boule})
    [
        internal_impurity_quantity,
        internal_impurity_quantity / internal_length_unit,
        internal_impurity_quantity / internal_length_unit^2
    ]
end

@. linear_exponential_boule(z,p) = p[1] + p[2]*z + p[3]*exp((z-p[4])/p[5])
fit_function(::Val{:linear_exponential_boule}) = linear_exponential_boule
impurity_density_model(::Val{:linear_exponential_boule}) = LinExpBouleImpurityDensity
fit_parameter_names(::Val{:linear_exponential_boule}) = Symbol.(["a", "b", "n", "l", "m"])

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
impurity_density_model(::Val{:parabolic_exponential_boule}) = ParExpBouleImpurityDensity
fit_parameter_names(::Val{:parabolic_exponential_boule}) = Symbol.(["a", "b", "c", "n", "l", "m"])

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