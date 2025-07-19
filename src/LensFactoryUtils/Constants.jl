module Constants

export RV, ROA
export CONST_C, CONST_G
export MASS_SUN, RADIUS_SUN, LUMINOSITY_SUN
export DIST_PC, DIST_KPC, DIST_MPC, DIST_GPC
export ANGLE_DEGREE, ANGLE_ARCMIN, ANGLE_ARCSEC, ANGLE_MILLIARCSEC, ANGLE_MICROARCSEC
export YEAR2SECOND, YEAR2HOUR, YEAR2DAY, DAY2HOUR, DAY2SECOND

# Manual definition of types that LensFactory can work with
const  RV = Union{Int64, Float64}
const ROA = Union{Vector{<:RV}, Matrix{<:RV}}

# Natural constants
const CONST_C = 299792458.0                         # Speed of light value (in m/s)
const CONST_G = 6.6743E-11                          # Gravitational constant (in m^3 kg/s^2)


# Solar constants
const MASS_SUN = 1.98855E+30                        # Mass of the Sun (in kg)
const RADIUS_SUN = 696340E+03                       # Radius of the Sun (in meters)
const LUMINOSITY_SUN = 3.828E26                     # Solar luminosity (in Watts = Joule/s)


# Distance units
const DIST_PC = 3.0856775714409184E+16              # One parsec (in meters)
const DIST_KPC = 3.0856775714409184E+19             # One Kilo-parsec (in meters)
const DIST_MPC = 3.0856775714409184E+22             # One Mega-parsec (in meters)
const DIST_GPC = 3.0856775714409184E+25             # One Giga-parsec (in meters)


# Angular units
const ANGLE_MICROARCSEC= 4.848136811095361E-12      # One micro-arcsec (in radian)
const ANGLE_MILLIARCSEC= 4.84813681109536E-09       # One milli-arcsec (in radian)
const ANGLE_ARCSEC = 4.84813681109536E-06           # One arcsec value (in radian)
const ANGLE_ARCMIN = 0.0002908882086657216          # One arcmin value (in radian)
const ANGLE_DEGREE = 0.017453292519943295           # One degree value (in radian)


# Temporal units
const YEAR2SECOND = 31557600.0                      # Number of seconds in one year
const YEAR2HOUR = 8766.0                            # Number of hours in one year
const YEAR2DAY = 365.25                             # Number of days in one year
const DAY2HOUR = 24.0                               # Number of hours in one day
const DAY2SECOND = 86400.0                          # Number of seconds in one day

end