import numpy as np

# constant coordinate format (alt,az) with np arrays

# initializing tables
sun_coords_ephimeris_previous = np.array([0,0])
sun_coords_ephimeris_now = np.array([0,0])
d_coords_movement = np.array([0,0])

# defining basic parameters
altCount = 10 # rows count
azCount = 10 # columns count
d_coords_scan = np.array([0.05]) # degrees

# matrix that will contain all the data
matrix = np.array([]) # will be altCount x azCount

def getSunCoordsNow():
    return np.array([1,1])

def telescopeMoveAltAz(dAltAz):
    moveAlt(dAltAz(0))
    moveAz(dAltAz(1))

    return

def getSdrSignalStrength():
    # ...
    return

# get initial (false) coordinates from IMU and point radiotelescope based on them
# SPEED = MINIMUM

initialCoords = getImuCoords()

dCoordsInitial = getSunCoordsNow() - initialCoords

telescopeMoveAltAz(dCoordsInitial)

# scanning
# SPEED = MAXIMUM

sun_coords_ephimeris_previous = getSunCoordsNow()

for i in range(altCount):
    for j in range(azCount):

        sun_coords_ephimeris_now = getSunCoordsNow() # the real coords of the sun

        d_coords_movement = sun_coords_ephimeris_now - sun_coords_ephimeris_previous # coords change due to sun's movement

        # d_coords_scan = constant, defined on the beginning

        d_coords_total = d_coords_movement + d_coords_scan # total coords change
        
        # Each time the coordinates will be incremented by d_coords_movement and by d_coords_scan

        telescopeMoveAltAz(d_coords_total)

        signalStrength = getSdrSignalStrength()

        theoreticalCoords = theoreticalCoords + d_coords_total # based on the initial coords form the IMU

        # t=time_now_ACCURATE # perhaps not needed

        matrix_to_append = np.array([sun_coords_ephimeris_now, theoreticalCoords, signalStrength])

        matrix.append(matrix_to_append)


# finding maximum signal strength

max_strength_index = matrix.find_max_strength_i

# finding the corrections !!!

corrections = sun_real_coords(max_strength_index) - theoreticalCoords(max_strength_index)

# will be 1 x 2 matrix

# so whenever we will need to find a point we will use the initial IMU coord, but corrected with the corrections
# keeping track of the movements of the radiotelescope, we will know where it points