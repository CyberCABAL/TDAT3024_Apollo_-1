import math
import numpy as np
import matplotlib.pyplot as plt

G = 6.67408 * 10**(-11)


class CelestialObject:
    def __init__(self,
                 position = [0, 0, 0],
                 dirVec = [0, 0, 0],
                 mass = 1, r = 1, name = ""):
        self.mass = mass
        self.gm = G * mass
        self.r = r
        self.position = np.array(position)
        self.dirVec = np.array(dirVec)
        self.name = name
        self.stop = False;

class Rocket(CelestialObject):
    def __init__(self,
                 position = [0, 0, 0],
                 dirVec = [0, 0, 0],
                 mass = 1,
                 A = 1,
                 name = "Saturn V",
                 #system_init = [0, 0, 0],
                 force = 1,
                 C_d = 0.5,
                 stop = False,
                 saturn_v = None):
        CelestialObject.__init__(self, position, dirVec, mass, 1, name);
        #self.system_init = system_init;
        self.A = A;
        self.s_V = saturn_v;
        self.force = force;
        self.C_d = C_d;
        self.stop = stop;

    def update(self, t):
        if (not self.stop):
            self.mass = self.s_V.get_mass(t);
            self.A = self.s_V.get_area(t);
            self.force = self.s_V.get_force(t);
        else:
            #self.mass = 0;
            self.force = 0;
            self.dirVec = [0] * len(self.dirVec);

    def a_R(self):
        if (self.stop):
            return np.array([0] * len(self.dirVec));
        l = np.linalg.norm(self.dirVec, 2);
        if (l == 0 or self.mass == 0):
            return np.array([0] * len(self.dirVec));
        return (self.dirVec / l) * (self.force / self.mass);

    def get_h(self, planet):  #Height
        return np.linalg.norm(planet.position - self.position, 2) - planet.r;

    def a_Atmos(self, planet):  #Resistance
        if (self.stop):
            return np.array([0] * len(self.dirVec));
        l = np.linalg.norm(np.array(self.dirVec), 2);
        if (l == 0 or self.mass == 0):
            return np.array([0] * len(self.dirVec));
        return -(self.dirVec / l) * (F_d_h(self.C_d, self.get_h(planet), self.A, l) / self.mass);
    

# Total force of gravity
def F_G(x_0, x_1):
    r = dist(x_0, x_1)
    return G * x_0.mass * x_1.mass / (r * r)


# Distance between two objects
def dist(x_0, x_1):
    p0 = x_0.position
    p1 = x_1.position
    delta = [p0[i] - p1[i] for i in range(len(p0))]
    distSum = 0
    for i in range(len(p0)):
        distSum += delta[i] * delta[i]
    return math.sqrt(distSum)


# Force of drag
def F_d(C_d, ρ_air, A, v):
    return (C_d * ρ_air * A * v * v) / 2

# Force of drag
def F_d_h(C_d, h, A, v):
    return (C_d * ρ_atmos_h(h) * A * v * v) / 2

# Density of atmosphere
def ρ_atmos(p_air, T):
    return (p_air * 3.4855) / (T * 1000)

# Density of atmosphere
def ρ_atmos_h(h):
    return (p_air(h) * 3.4855) / (T(h) * 1000);


# Temperature at height h
def T(h):
    T_a = 1
    T_b = 1
    if h < 11000:
        T_a = 288.19
        T_b = -0.00649
    elif h < 25000:
        T_a = 216.69
        T_b = 0
    elif h >= 25000:
        T_a = 141.94
        T_b = 0.00299
    return T_a + T_b * h


# Pressure at height h
def p_air(h):
    p_a = 1
    p_b = 1
    exp = 1
    if h < 11000:
        p_a = 101290
        p_b = 288.08
        exp = 5.256
    elif h < 25000:
        return 101290 * math.pow(math.e, -0.000157 * h)
    else:
        p_a = 2488
        p_b = 216.6
        exp = -11.388
    return p_a * math.pow(T(h) / p_b, exp)


if __name__ == "__main__":
    # execute only if run as a script
    x = np.arange(0, 50000, 5000)
    pressure, temp, density = np.ones(10), np.ones(10), np.ones(10)
    for i in range(10):
        temp[i] = T(x[i])/100  # unit = 100K
        pressure[i] = p_air(x[i])/10000  # unit = 10kPa
        density[i] = ρ_atmos(pressure[i], temp[i])
    plt.plot(x, pressure, 'b-', x, temp, 'r-', x, density, 'y-')
    plt.ylabel('Temp (100K) - Trykk (10kPa) - Tetthet(100g/m3)')
    plt.xlabel('Høyde (m)')
    plt.show()
