import os

def makespec(timeSlots=100, timeStep=0.1, horizon=1, gridSize=200, smallNum=10000,
        smallMass=0.0001, smallRadius=0.0001, largeNum=0, largePtc=None):
    spec = '''TimeSlots: {TimeSlots}
TimeStep: {TimeStep}
Horizon: {Horizon}
GridSize: {GridSize}
NumberOfSmallParticles: {NumberOfSmallParticles}
SmallParticleMass: {SmallParticleMass}
SmallParticleRadius: {SmallParticleRadius}
NumberOfLargeParticles: {NumberOfLargeParticles}
'''.format(TimeSlots=timeSlots, TimeStep=timeStep, Horizon=horizon,
            GridSize=gridSize, NumberOfSmallParticles=smallNum, 
            SmallParticleMass=smallMass, SmallParticleRadius=smallRadius,
            NumberOfLargeParticles=largeNum)
    if largeNum == 0:
        return spec
    else:
        ptcs = map(lambda x : [str(y) for y in x], largePtc)
        return spec + '\n'.join([' '.join(ptc) for ptc in ptcs])

if __name__ == "__main__":
    os.system("make clean && make")
    largePtc = [[8,1.5,56.5,70.21],[6,1,156,20.1],[12,3.1,20,30],[8,1.7,180,90]]
    os.system("cp spec.txt spec_backup.txt")
    try:
        with open("speedup.txt", "w+", 1) as resFile:
            for horizon in [0,1,2,3]:
                for smallNum in [50, 100, 1000, 5000, 10000]:
                    for np in [1, 4, 9, 16, 25, 36, 49, 64]:
                        with open("spec.txt", "w+", 1) as spec:
                            spec.write(makespec(horizon=horizon, smallNum=smallNum, 
                                    largeNum=4, largePtc=largePtc))
                        excuTime = os.popen('./run.sh {}'.format(np)).read().strip()
                        res = "{}, {}, {}, {}\n".format(horizon, smallNum, np, excuTime)
                        resFile.write(res)
                        print(res)
            for horizon in [0,1,2,3]:
                np = 16
                for smallNum in [20000,50000]:
                    with open("spec.txt", "w+", 1) as spec:
                        spec.write(makespec(horizon=horizon, smallNum=smallNum, 
                                largeNum=4, largePtc=largePtc))
                    excuTime = os.popen('./run.sh {}'.format(np)).read().strip()
                    res = "{}, {}, {}, {}\n".format(horizon, smallNum, np, excuTime)
                    resFile.write(res)
                    print(res)
    except KeyboardInterrupt:
        os.system("mv ppmresults results")
        os.system("make clean")
        os.system("mv spec_backup.txt spec.txt")