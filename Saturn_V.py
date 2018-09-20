'''
Saturn V info: 
Totalvekt:  2970000  kg

Steg:  1 
Totalvekt:  2290000  kg 
Vekt uten drivstoff:  130000  kg 
Drivstoff i kg:  2160000  kg 
Masseendring:  -12857.142857142857  kg/s 
Hastighet:  2730.0  m/s
Skyvekraft:  35100000  N

Steg:  2 
Totalvekt:  496200.0  kg 
Vekt uten drivstoff:  40100.0  kg 
Drivstoff i kg:  456100.0  kg 
Masseendring:  -1266.9444444444443  kg/s 
Hastighet:  4057.7943433457576  m/s
Skyvekraft:  5141000  N

Steg:  3 
Totalvekt:  123000  kg
Vekt uten drivstoff:  13500.0  kg 
Drivstoff i kg:  109500.0  kg 
Masseendring:  -219.0  kg/s 
Hastighet:  4566.2100456621  m/s
Skyvekraft:  1000000  N
'''

#BTW så e d bare å endre variabelnavn og metodenavn
class Saturn_V:
    '''
    class Saturn_V(self, totmasse, stegtotmasse, stegtorrmasse,
        stegskyvekraft, stegtid)
    totmasse: Den totale massen til raketten
    stegtotmasse: En tabell med total masse for hvert steg
    stegtorrmasse: En tabell med massen til hvert steg uten drivstoff
    stegskyvekraft: En tabell med skyvekraften i hvert steg
    stegtid: En tabell med tiden hvert steg tar

    standardverdiene som er lagt inn i konstruktøren er hentet fra
    https://en.wikipedia.org/wiki/Saturn_V
    '''
    def __init__(self,
               totmasse = 2970*10**3,
               stegtotmasse = [2290*10**3, 496.2*10**3, 123*10**3],
               stegtorrmasse = [130*10**3, 40.1*10**3, 13.5*10**3],
               stegskyvekraft = [35100*10**3, 5141*10**3, 1000*10**3],
               stegtid = [168, 360, 500]):
        self.m = totmasse;
        self.stegM = stegtotmasse;
        self.stegTorrM = stegtorrmasse;
        self.stegF = stegskyvekraft;
        self.stegT = stegtid;

    def getDrivstoff(self, x):
        '''Finner hvor mye drivstoff et steg har i kg'''
        return self.stegM[x]-self.stegTorrM[x]

    def getMasseEndring(self, x):
        '''Finner masseendringen til et steg i kg/s'''
        return -(self.getDrivstoff(x)/self.stegT[x])

    def getV(self,x):
        '''Estimerer hastigheten til eksosgassene i hvert trinn'''
        return self.stegF[x]/-(self.getMasseEndring(x))

    def stegInfo(self, stegNr):
         steg=stegNr-1
         print("Steg: ", stegNr, "\nTotalvekt: ", self.stegM[steg], " kg",
               "\nVekt uten drivstoff: ", self.stegTorrM[steg], " kg",
               "\nDrivstoff i kg: ", self.getDrivstoff(steg), " kg",
               "\nMasseendring: ", self.getMasseEndring(steg), " kg/s",
               "\nHastighet: ", self.getV(steg), " m/s",
               "\nSkyvekraft: ", self.stegF[steg], " N")

    def saturn_vInfo(self):
        print("Saturn V info:", "\nTotalvekt: ", self.m, " kg")
        for i in range(3):
            print("\n")
            self.stegInfo(i+1)

    def getStegStartMasse(self,x):
        '''finner totalmassen til raketten rett før hver trinn tennes'''
        if(x<3 and x>=0):
            if(x==0):
                return self.m
            if(x==1):
                return self.m-self.stegM[0]
            else:
                return self.m-self.stegM[0]-self.stegM[1]
            
    def mt(self, t):
        '''En funksjon m(t) som gir massen til raketten ved
        tiden t sekunder etter oppskytingen'''
        if(t>=0 and t<=1028):
            
            if(t<=168):
                return self.getStegStartMasse(0)+self.getMasseEndring(0)*t
            if(t>168 and t<=528):
                st=t-168
                return self.getStegStartMasse(1)+self.getMasseEndring(1)*st
            else:
                st=t-168-360
                return self.getStegStartMasse(2)+self.getMasseEndring(2)*st
        else:
            if(t>1028):
                return self.m-self.stegM[0]-self.stegM[1]-self.stegM[2]
            else:
                raise ValueError('velg en t i intervallet [0,->]')

    def skyvekraft(self, t):
        '''En funksjon som gir skyvekraften til raketten ved
        tiden t sekunder etter oppskytingen, siden vi antar
        at skyvekraften er konstant i hvert steg, så
        sjekker metoden bare hvilket steg det er snakk om'''
        if(t>=0 and t<=1028):
            if(t<=168):
                return self.stegF[0]
            if(t>168 and t<=528):
                return self.stegF[1]
            else:
                return self.stegF[2]
        else:
            if(t>1028):
                return 0
            else:
                raise ValueError('velg en t i intervallet [0, ->]')
        

saturnv = Saturn_V();
#saturnv.saturn_vInfo()
print(saturnv.mt(300))
print(saturnv.skyvekraft(600))

