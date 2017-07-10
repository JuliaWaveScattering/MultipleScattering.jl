# UNFINISHED
# If you are feeling very adventurous, you can define your own shape
import MultipleScattering

type MyShape <: MultipleScattering.Shape
end

MultipleScattering.volume(shape::MyShape) = 0.0

MultipleScattering.name(shape::MyShape) = "MyShape"

MultipleScattering.bounding_box(shape::MyShape) = MultipleScattering.Rectangle()