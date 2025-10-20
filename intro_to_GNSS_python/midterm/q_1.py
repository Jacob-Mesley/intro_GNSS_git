# work from question one of the midterm

# imports
import numpy as np
import convert
import calculate

# given constants to use
c = 3*10**8  # [m/s]
mu_e = 398600.5 * 1000**3  # [m^2/s^3]

# QUESTION 1 -------------------------------------------------------------------------


# PART A ----------------------
print("\n-------------- PART A --------------")

# given sat positions in ECEF
sat1 = [1193698, -4454942, 26156494]
sat2 = [832273, -4720052, 4792867]
sat3 = [-3994192, 22652188, 13280000]
sat4 = [21189069, -36700543, 0]
sat5 = [7867022, -21614466, 13280000]
sats = [sat1, sat2, sat3, sat4, sat5]

# find the magnitude of their position vector, and magnitude of velocity
sats_pos_mag = [0, 0, 0, 0, 0]
sats_velo_mag = [0, 0, 0, 0, 0]
print("Velocity Magnitudes for each sat")
for i in range(len(sats)):
    sats_pos_mag[i] = np.sqrt(sats[i][0]**2 + sats[i][1]**2 + sats[i][2]**2)
    sats_velo_mag[i] = np.sqrt(mu_e / sats_pos_mag[i])
    # print("sat"+str(i+1)+"_pos_mag -> " + str(sats_pos_mag[i]) + " [m]")
    print("sat"+str(i+1)+" -> " + str(sats_velo_mag[i]) + " [m/s]")


# PART B -----------------------
print("\n-------------- PART B --------------")

# given observer position in lat, long, altitude
observer_LLA = [26, 80, 0]

# convert the lat long altitude (LLA) to ECEF using function from hw 2
observer_ECEF = convert.lla2ecef(observer_LLA)
print("Observer ECEF Position [m]")
print(str(np.array(observer_ECEF)))


# PART C ------------------------
print("\n-------------- PART C --------------")

# calculate the al_el_range for each position
print("Az [deg], El [deg] and Range [m] for each sat, DO NOT account for visibility")
for i in range(len(sats)):
    az_el_range = calculate.az_el_range(np.array([observer_ECEF]), np.array([sats[i]]), only_when_visible=False)
    print("sat" + str(i + 1) + " -> " + str(az_el_range))

# PART D --------------------------------
print("\n-------------- PART D --------------")

print("\nSat Position Magnitude (GPS has an altitude of aprox ~20,200 km, so plus earth radius ~26578 km)")
for i in range(len(sats)):
    print("sat"+str(i+1)+" -> "+str(sats_pos_mag[i]/1000)+" [km]")
