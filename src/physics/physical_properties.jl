"""
Holds information about the physical properties of the medium, the dimension of
the field and the number of dimensions it is a function of.
"""
abstract type PhysicalProperties{Dim,FieldDim,T<:AbstractFloat} end

"Extract the dimension of the field of this physical property"
field_dim(p::PhysicalProperties{Dim,FieldDim,T}) where {Dim,FieldDim,T} = FieldDim

"Extract the dimension of the field of this type of physical property"
field_dim(p::Type{P}) where {Dim,FieldDim,T,P<:PhysicalProperties{Dim,FieldDim,T}} = FieldDim

"Extract the dimension of the space that this physical property lives in"
dim(p::PhysicalProperties{Dim,FieldDim,T}) where {Dim,FieldDim,T} = Dim

"Extract the dimension of the space that this type of physical property lives in"
dim(p::Type{P}) where {Dim,FieldDim,T,P<:PhysicalProperties{Dim,FieldDim,T}} = Dim

"""
Basis functions in a specific dimension for a specific physics type.
"""
function get_basis_function(medium::PhysicalProperties, ω::T) where {T}
    error("Get basis function not implmented for this physics type.")
end
