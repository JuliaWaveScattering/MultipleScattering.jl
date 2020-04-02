"""
Holds information about the physical properties of the medium, the dimension of
the field and the number of dimensions it is a function of.
"""
abstract type PhysicalMedium{T<:AbstractFloat,Dim,FieldDim} end

"Extract the dimension of the field of this physical property"
field_dimension(p::PhysicalMedium{T,Dim,FieldDim}) where {T,Dim,FieldDim} = FieldDim

"Extract the dimension of the field of this type of physical property"
field_dimension(p::Type{P}) where {Dim,FieldDim,T,P<:PhysicalMedium{T,Dim,FieldDim}} = FieldDim

"Extract the dimension of the space that this physical property lives in"
spatial_dimension(p::PhysicalMedium{T,Dim,FieldDim}) where {Dim,FieldDim,T} = Dim

"Extract the dimension of the space that this type of physical property lives in"
spatial_dimension(p::Type{P}) where {Dim,FieldDim,T,P<:PhysicalMedium{T,Dim,FieldDim}} = Dim

"""
Basis functions in a specific dimension for specific physcial properties..
"""
function regular_basis_function(medium::P, ω::T) where {P<:PhysicalMedium,T}
    error("No basis function implmented for this physics type.")
end

"""
Basis of outgoing wave functions in a specific dimension for specific physcial properties..
"""
function outgoing_basis_function(medium::P, ω::T) where {P<:PhysicalMedium,T}
    error("No basis function implmented for this physics type.")
end

"""
the field inside an AbstractParticle a some given point x.
"""
internal_field

"""
A tuples of vectors of the field close to the boundary of the shape. The field is calculated from sim::FrequencySimulation, but the PhysicalMedium inside and outside of the shape are assumed to be given by inside_medium and outside_medium.
"""
boundary_data
