import math
import numpy as np

G = 6.67408 * 10**(-11);

class CelestialObject:
    def __init__(self,
                 position = [0, 0, 0],
                 dirVec = [0, 0, 0],
                 mass = 1, r = 1, name = ""):
        self.mass = mass;
        self.gm = G * mass;
        self.r = r;
        self.position = np.array(position);
        self.dirVec = np.array(dirVec);
        self.name = name;

# Total force of gravity
def F_G(x_0, x_1):
    r = dist(x_0, x_1);
    return G * x_0.mass * x_1.mass / (r * r);

# Distance between two objects
def dist(x_0, x_1):
    p0 = x_0.position;
    p1 = x_1.position;
    delta = [p0[i] - p1[i] for i in range(len(p0))];
    distSum = 0;
    for i in range(len(p0)):
        distSum += delta[i] * delta[i];
    return math.sqrt(distSum);

# Force of drag
def F_d(C_d, p_air, A, v):
    return (C_d * p_air * A * v * v) / 2;

# Density of atmosphere
def ρ_atmos(p_air, T):
    return p_air * 3.4855 / T;

# Temperature at height h
def T(h):
    T_a = 1;
    T_b = 1;
    if (h < 11000):
        T_a = 288.19;
        T_b = -0.00649;
    elif (h < 25000):
        return 216.69;
    else:
        T_a = 141,94;
        T_b = 0.00299;
    return T_a + T_b * h;

# Pressure at height h
def p_air(h):
    p_a = 1;
    p_b = 1;
    exp = 1;
    if (h < 11000):
        p_a = 101290;
        p_b = 288.08;
        exp = 5.256;
    elif (h < 25000):
        return 101290 * math.pow(math.e, -0.000157 * h);
    else:
        p_a = 2488;
        p_b = 216.6;
        exp = -11.388;
    return p_a * math.pow(T(h) / p_b, exp);
