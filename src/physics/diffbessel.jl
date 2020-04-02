# "Derivative of Hankel function of zeroth kind"
# function diffhankelh1(n, z)
#     if n != 0
#         (hankelh1(-1 + n, z) - hankelh1(1 + n, z)) / 2
#     else
#         -hankelh1(one(n), z)
#     end
# end
#
# "Derivative of Bessel function of first kind"
# function diffbesselj(n, z)
#     if n != 0
#         (besselj(-1 + n, z) - besselj(1 + n, z)) / 2
#     else
#         -besselj(one(n), z)
#     end
# end
