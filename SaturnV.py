'''
Saturn V info: 
Totalvekt:  2970000  kg

Steg:  1 
Totalvekt:  2290000  kg 
Vekt uten drivstoff:  130000  kg 
Drivstoff i kg:  2160000  kg
Varighet i s:  168  s 
Masseendring:  -12857.14 kg/s 
Hastighet:  2730.0  m/s
Skyvekraft:  35100000  N

Steg:  2 
Totalvekt:  496200.0  kg 
Vekt uten drivstoff:  40100.0  kg 
Drivstoff i kg:  456100.0  kg
Varighet i s:  360  s 
Masseendring:  -1266.94 kg/s 
Hastighet:  4057.79  m/s
Skyvekraft:  5141000  N

Steg:  3 
Totalvekt:  123000  kg
Vekt uten drivstoff:  13500.0  kg 
Drivstoff i kg:  109500.0  kg
Varighet i s:  500  s 
Masseendring:  -219.0  kg/s 
Hastighet:  4566.21  m/s
Skyvekraft:  1000000  N
'''


class SaturnV(object):
    '''
    class SaturnV(self, totmasse, stegtotmasse, stegtorrmasse, stegskyvekraft, stegtid)
        total_mass: Den totale massen til raketten
        step_total_mass: En tabell med total masse for hvert steg
        step_mass_without_fuel: En tabell med massen til hvert steg uten drivstoff
        step_force: En tabell med skyvekraften i hvert steg
        step_time: En tabell med tiden hvert steg tar

    Standardverdiene som er lagt inn i konstruktøren er hentet fra https://en.wikipedia.org/wiki/Saturn_V
    '''
    def __init__(self,
                 total_mass=2970 * 10 ** 3,
                 step_total_mass=[2290 * 10 ** 3, 496.2 * 10 ** 3, 123 * 10 ** 3],
                 step_mass_without_fuel=[130 * 10 ** 3, 40.1 * 10 ** 3, 13.5 * 10 ** 3],
                 step_force=[35100 * 10 ** 3, 5141 * 10 ** 3, 1000 * 10 ** 3],
                 step_time=[168, 360, 500]):
        self.mass = total_mass
        self.step_mass = step_total_mass
        self.step_mass_without_fuel = step_mass_without_fuel
        self.step_force = step_force
        self.step_time = step_time

    def get_fuel(self, x):
        # Finner hvor mye drivstoff et steg har i kg
        return self.step_mass[x] - self.step_mass_without_fuel[x]

    def get_mass_difference(self, x):
        # Finner masseendringen til et steg i kg/s
        return -(self.get_fuel(x) / self.step_time[x])

    def get_velocity(self, x):
        # Estimerer hastigheten til eksosgassene i hvert trinn
        return self.step_force[x] / -(self.get_mass_difference(x))

    def step_info(self, step_number):
        steg = step_number - 1
        print("Steg: ", step_number, "\nTotalvekt: ", self.step_mass[steg], " kg",
              "\nVekt uten drivstoff: ", self.step_mass_without_fuel[steg], " kg",
              "\nDrivstoff i kg: ", self.get_fuel(steg), " kg",
              "\nVarighet i s: ", self.step_time[steg], " s",
              "\nMasseendring: ", round(self.get_mass_difference(steg),2), " kg/s",
              "\nHastighet: ", round(self.get_velocity(steg),2), " m/s",
              "\nSkyvekraft: ", self.step_force[steg], " N")

    def saturn_v_info(self):
        print("Saturn V info:", "\nTotalvekt: ", self.mass, " kg")
        for i in range(3):
            print("\n")
            self.step_info(i + 1)

    def mass_table(self):
        print("Steg 1:")
        for i in range(4):
            print("t:", i*50, "s", "\tmasse:", round(self.get_mass(i*50),0), "kg")
        print("t: 168 s", "\tmasse:", round(self.get_mass(168),0), "kg", "\n\nSteg 2:")
        print("t: 168.1 s", "\tmasse:", round(self.get_mass(168.1),0), "kg")
        for i in range(4,11):
              print("t:", i*50, "s", "\tmasse:", round(self.get_mass(i*50),0), "kg")
        print("t: 528 s", "\tmasse:", round(self.get_mass(528),0), "kg", "\n\nSteg 3:")
        print("t: 528.1 s", "\tmasse:", round(self.get_mass(528.1),0), "kg")
        for i in range (11, 21):
            print("t:", i*50, "s", "\tmasse:", round(self.get_mass(i*50),0), "kg")
        print("t: 1028 s", "\tmasse:", round(self.get_mass(1028),0), "kg")
        print("t: 1028.1 s", "\tmasse:", round(self.get_mass(1028.1),0), "kg")

    def get_step_start_mass(self, x):
        # Finner totalmassen til raketten rett før hver trinn tennes
        if 0 <= x < 3:
            if x == 0:
                return self.mass
            if x == 1:
                return self.mass - self.step_mass[0]
            else:
                return self.mass - self.step_mass[0] - self.step_mass[1]
            
    def get_mass(self, t):
        # En funksjon m(t) som gir massen til raketten ved tiden t sekunder etter oppskytingen
        if 0 <= t <= (self.step_time[0]+self.step_time[1]+self.step_time[2]):
            if t <= self.step_time[0]:
                return self.get_step_start_mass(0) + self.get_mass_difference(0) * t
            if self.step_time[0] < t <= (self.step_time[0]+self.step_time[1]):
                st = t-self.step_time[0]
                return self.get_step_start_mass(1) + self.get_mass_difference(1) * st
            else:
                st = t-self.step_time[0]-self.step_time[1]
                return self.get_step_start_mass(2) + self.get_mass_difference(2) * st
        else:
            if t > (self.step_time[0]+self.step_time[1]+self.step_time[2]):
                return self.mass - self.step_mass[0] - self.step_mass[1] - self.step_mass[2]
            else:
                raise ValueError('velg en t i intervallet [0,->]')

    def get_force(self, t):
        # En funksjon som gir skyvekraften til raketten ved tiden t sekunder etter oppskytingen.
        # Siden vi antar at skyvekraften er konstant i hvert steg, så sjekker metoden bare hvilket steg det er snakk om
        if 0 <= t <= (self.step_time[0]+self.step_time[1]+self.step_time[2]):
            if t <= self.step_time[0]:
                return self.step_force[0]
            if self.step_time[0] < t <= (self.step_time[0]+self.step_time[1]):
                return self.step_force[1]
            else:
                return self.step_force[2]
        else:
            if t > (self.step_time[0]+self.step_time[1]+self.step_time[2]):
                return 0
            else:
                raise ValueError('Velg en t i intervallet [0, ->]')

def main():
    saturn_v = SaturnV()
    saturn_v.saturn_v_info()

    print("\nTotalmasse i starten av steg 1: ",saturn_v.get_step_start_mass(0), "kg")
    print("\nTotalmasse i starten av steg 2: ",saturn_v.get_step_start_mass(1), "kg")
    print("\nTotalmasse i starten av steg 3: ",saturn_v.get_step_start_mass(2), "kg")

    print("\nMasse ved ulike t-verdier")
    saturn_v.mass_table()

if __name__ == "__main__":
    main();
        



