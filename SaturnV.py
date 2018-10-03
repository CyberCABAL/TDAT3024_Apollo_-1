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

import math

class SaturnV(object):
    '''
    class SaturnV(self, totmasse, stegtotmasse, stegtorrmasse, stegskyvekraft, stegtid)
        total_mass: Den totale massen til raketten
        step_total_mass: En tabell med total masse for hvert steg
        step_mass_without_fuel: En tabell med massen til hvert steg uten drivstoff
        step_force: En tabell med skyvekraften i hvert steg
        step_time: En tabell med tiden hvert steg tar
        step_diameter: En tabell med diameteren i hvert steg. Den siste verdien er diametern i kommandomodulen

    Standardverdiene som er lagt inn i konstruktøren er hentet fra https://en.wikipedia.org/wiki/Saturn_V
    '''
    def __init__(self,
                 total_mass=2970 * 10 ** 3,
                 step_total_mass=[2290 * 10 ** 3, 496.2 * 10 ** 3, 123 * 10 ** 3],
                 step_mass_without_fuel=[130 * 10 ** 3, 40.1 * 10 ** 3, 13.5 * 10 ** 3],
                 step_force=[35100 * 10 ** 3, 5141 * 10 ** 3, 1000 * 10 ** 3],
                 step_time=[168, 360, 500],
                 step_diameter=[10.1, 10.1, 6.6]):
        self.mass = total_mass
        self.step_mass = step_total_mass
        self.step_mass_without_fuel = step_mass_without_fuel
        self.step_force = step_force
        self.step_time = step_time
        self.step_diameter = step_diameter

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
        tempmass=self.mass
        for step in range (len(self.step_time)):
            if x==step:
                return tempmass
            tempmass -= self.step_mass[step]
            
    def get_mass(self, t):
        # En funksjon m(t) som gir massen til raketten ved tiden t sekunder etter oppskytingen
        st = t;
        steptime=self.step_time[0]
        for step in range(len(self.step_time)):
            if t <= steptime:
                return self.get_step_start_mass(step) + self.get_mass_difference(step) * st
            st -= self.step_time[step]
            steptime += self.step_time[step+1] if len(self.step_time) > step + 1 else 0
        return self.mass - sum(self.step_mass)
        
    def get_force(self, t):
        # En funksjon som gir skyvekraften til raketten ved tiden t sekunder etter oppskytingen.
        # Siden vi antar at skyvekraften er konstant i hvert steg, så sjekker metoden bare hvilket steg det er snakk om
        steptime=self.step_time[0]
        for step in range(len(self.step_time)):
            if t <= steptime:
                return self.step_force[step]
            steptime += self.step_time[step+1] if len(self.step_time) > step + 1 else 0
        return 0

    def get_area(self, t):
        # Denne funksjonen ble lagt til med tanke på atmosfæren i oppgave 6.
        # Den returnerer arealet til raketten etter t sekunder.
        steptime=self.step_time[0]
        for step in range(len(self.step_time)):
            if t <= steptime:
                return (((self.step_diameter[step]/2)**2)*math.pi)
            steptime += self.step_time[step+1] if len(self.step_time) > step + 1 else 0
        return (((3.9/2)**2)*math.pi) #area of comando module
    
def main():
    saturn_v = SaturnV()
    saturn_v.saturn_v_info()
    #print(saturn_v.get_area(1100))

    print("\nTotalmasse i starten av steg 1: ",saturn_v.get_step_start_mass(0), "kg")
    print("\nTotalmasse i starten av steg 2: ",saturn_v.get_step_start_mass(1), "kg")
    print("\nTotalmasse i starten av steg 3: ",saturn_v.get_step_start_mass(2), "kg")

    print("\nMasse ved ulike t-verdier")
    saturn_v.mass_table()

    

if __name__ == "__main__":
    main();
        



