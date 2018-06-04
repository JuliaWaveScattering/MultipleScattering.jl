"""
Holds information about the physical properties of the medium, the dimension of
the field and the number of dimensions it is a function of.
"""
abstract type PhysicalProperties{T<:AbstractFloat,Dim,FieldDim} end

"Extract the dimension of the field of this physical property"
field_dim(p::PhysicalProperties{T,Dim,FieldDim}) where {T,Dim,FieldDim} = FieldDim

"Extract the dimension of the field of this type of physical property"
field_dim(p::Type{P}) where {Dim,FieldDim,T,P<:PhysicalProperties{T,Dim,FieldDim}} = FieldDim

"Extract the dimension of the space that this physical property lives in"
dim(p::PhysicalProperties{T,Dim,FieldDim}) where {Dim,FieldDim,T} = Dim

"Extract the dimension of the space that this type of physical property lives in"
dim(p::Type{P}) where {Dim,FieldDim,T,P<:PhysicalProperties{T,Dim,FieldDim}} = Dim

"""
Basis functions in a specific dimension for a specific physics type.
"""
function basis_function(medium::PhysicalProperties, ω::T) where {T}
    error("No basis function implmented for this physics type.")
end

"""
the field inside an AbstractParticle a some given point x.  
"""
internal_field

"""
Returns a function that gives the value of the besselj expansion centred at centre
"""
besselj_field
