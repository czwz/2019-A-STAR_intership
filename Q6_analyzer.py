###########################################################################################
####This file is specifically designed for the LAMMPS trajectory file to calculate ALBO####
####(Averaged Local Bond Order parameter). It requires the output format to accuire Q4-####
####Q12 plus every Ylm-s for l=6.                                                      ####
###########################################################################################
import numpy as np

def distance(x0, x1, dimensions):
    delta = np.abs(x0 - x1)
    delta = np.where(delta >= 0.5 * dimensions, delta - dimensions, delta)
    return np.sqrt((delta ** 2).sum(axis=-1))

class trj_analyzer: 
    
    def __init__(self, filename):
        self.boxes = []
        self.positions = []
        self.Ylm = []
        self.Q6 = []
        self.atoms = 0
        self.timesteps = 0
        self.cutoff = 5.32
                
        with open(str(filename)) as fileobject:
            line_number = 1
    
            for line in fileobject:                    
                line = line.split()
                if line_number%(self.atoms+9) == 1:
                    box, position, Y, Q = [], [], [], []
                    
                if line_number==4:
                    self.atoms = int(line[0]) 
                    
                if len(line) == 2:
                    try:
                        box.append(float(line[1]) - float(line[0]))
                    except:
                        pass
            
                if len(line) > 5:
                    try:
                        line[0] = int(line[0])
                        position.append(np.array([float(line[3]),float(line[4]),float(line[5])]))
                        Q.append(np.array([float(line[7])]))
                        Y.append(np.array([float(line[11]),float(line[12]),float(line[13]),float(line[14]),float(line[15]),float(line[16]),
                                           float(line[17]),float(line[18]),float(line[19]),float(line[20]),float(line[21]),float(line[22]),
                                           float(line[23]),float(line[24]),float(line[25]),float(line[26]),float(line[27]),float(line[28]),
                                           float(line[29]),float(line[30]),float(line[31]),float(line[32]),float(line[33]),float(line[34]),
                                           float(line[35]),float(line[36])]
                                         )
                                )
                    except Exception:
                        pass
            
                if line_number%(self.atoms+9) == 0:
                    self.positions.append(np.array(position))
                    self.Ylm.append(np.array(Y))
                    self.Q6.append(np.array(Q))
                    self.boxes.append(np.array(box))
                line_number = line_number + 1
            self.timesteps = len(self.boxes)
            
    def ALBO(self):
        new_Q6 = np.zeros((self.timesteps,self.atoms))

        for t in range(self.timesteps):
            Xlm = np.zeros_like(self.Ylm[t])
            Ylm_denormalized = np.zeros_like(self.Ylm[t][0])
            Nb = np.zeros(self.atoms)
            dist_matrix = []
    
            for i in range(self.atoms):
                dist_matrix = distance(self.positions[t][i], self.positions[t], self.boxes[t])
                for j in range(self.atoms):
                    if dist_matrix[j] <= self.cutoff:
                        Xlm[i] = Xlm[i] + self.Ylm[t][j]
                        Nb[i] = Nb[i] + 1
            
            for i in range(self.atoms):
                Ylm_denormalized = Xlm[i]*( self.Q6[t][i]**2 / (4*np.pi*(1/13)) )**0.5
                new_Q6[t][i] = ( 4*np.pi*(1/13)*np.sum( (Ylm_denormalized/Nb[i])**2 ) )**0.5
                
            return new_Q6
        
    def OUTPUT_ALBO_trj(self, New_Q6, old_file, new_file):
        f = open(str(new_file), "w")
        with open(str(old_file)) as fileobject:
            index = 0    
            for wholeline in fileobject:
                line = wholeline.split()     
                if len(line) > 6:
                    try:
                        step = index//self.atoms
                        atom = index%self.atoms
                        int( line[0] )
                        print(line[0], line[1], line[2], 
                              self.positions[step][atom][0], self.positions[step][atom][1], self.positions[step][atom][2], 
                              New_Q6[step][atom], float(line[37]), file=f)
                        index += 1
                    except Exception:
                        print(line[0], line[1], line[2], line[3], line[4], line[5], line[6], line[7], line[8], line[39], file=f)
                else:
                    print(wholeline, file=f, end = '')
        f.close()
