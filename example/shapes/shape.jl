# # Shapes

# ## Existing shapes
# The package provides 3 built in basic shapes to put your random particles in, 
# you can plot them using:

using MultipleScattering

rectangle = Rectangle([0.0,0.0],[2.0,3.0])
plot(rectangle)

circle = Circle(3.0,[-1.0,2.0])
plot!(circle)

timeofflight = TimeOfFlight([-1.0,0.0],3.0)
plot!(timeofflight)

# Time of flight is a shape which contains shapes from a half space which take at
# most `t` time to reach from the listener.

# ## New shape
# If you are feeling very adventurous, you can define your own shape
# First you must import the package in order to add to existing functions

import MultipleScattering

type MyShape <: MultipleScattering.Shape
end


# To describe the characteristics and behaviour of the function you must define
# the following functions:

MultipleScattering.volume(shape::MyShape) = 0.0

MultipleScattering.name(shape::MyShape) = "MyShape"

MultipleScattering.bounding_box(shape::MyShape) = MultipleScattering.Rectangle()

# When you have this, you can make use of your shape to generate particles in it




