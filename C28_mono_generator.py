###########################################################################################
####This file is specifically designed to generate the xyz file for monoclinic crystall####
####structure of C28. It also provide the script for VMD-topo tool to create LAMMPS dat####
####a file.                                                                            ####
###########################################################################################
import numpy as np
from scipy.spatial.transform import Rotation as R

class C28_mono:
    def __init__(self):
        self.backbone = np.array([ [ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00],
                                   [ 8.55493930e-01, -2.42061080e-01,  1.25743253e+00],
                                   [ 0.00000000e+00,  0.00000000e+00,  2.51434564e+00],
                                   [ 8.55493930e-01, -2.42061080e-01,  3.77177817e+00],
                                   [-8.16520000e-04, -2.88740000e-04,  5.02919121e+00],
                                   [ 8.54677410e-01, -2.42349830e-01,  6.28662374e+00],
                                   [-2.61000000e-04, -4.45920000e-04,  7.54435337e+00],
                                   [ 8.55232930e-01, -2.42507010e-01,  8.80178590e+00],
                                   [-9.20340000e-04, -1.67908000e-03,  1.00589102e+01],
                                   [ 8.54573590e-01, -2.43740170e-01,  1.13163427e+01],
                                   [-9.20340000e-04, -1.67908000e-03,  1.25732558e+01],
                                   [ 8.54573590e-01, -2.43740170e-01,  1.38306884e+01],
                                   [-1.18133000e-03, -2.12501000e-03,  1.50889179e+01],
                                   [ 8.54312600e-01, -2.44186090e-01,  1.63463504e+01],
                                   [-1.18133000e-03, -2.12501000e-03,  1.76032636e+01],
                                   [ 8.54312600e-01, -2.44186090e-01,  1.88606961e+01],
                                   [-1.99785000e-03, -2.41375000e-03,  2.01181091e+01],
                                   [ 8.53496080e-01, -2.44474840e-01,  2.13755417e+01],
                                   [-1.44233000e-03, -2.57093000e-03,  2.26332713e+01],
                                   [ 8.54208780e-01, -2.45576430e-01,  2.38904151e+01],
                                   [-2.10167000e-03, -3.80409000e-03,  2.51478281e+01],
                                   [ 8.53392260e-01, -2.45865180e-01,  2.64052606e+01],
                                   [-2.10167000e-03, -3.80409000e-03,  2.76621738e+01],
                                   [ 8.53392260e-01, -2.45865180e-01,  2.89196063e+01],
                                   [-2.36266000e-03, -4.25002000e-03,  3.01778358e+01],
                                   [ 8.53131260e-01, -2.46311100e-01,  3.14352684e+01],
                                   [-2.36266000e-03, -4.25002000e-03,  3.26921815e+01],
                                   [ 8.53131260e-01, -2.46311100e-01,  3.39496140e+01]
                                ])
        self.bond = 5
        self.bond_z = 3.8
        backbone_axis = ( self.backbone[2] - self.backbone[0] ) / np.linalg.norm( self.backbone[2] - self.backbone[0] )
        rotation_axis = np.array( [1.0, 0.0, 0.0] )
        rotation_degree_right = -(30/180)*np.pi
        rotation_degree_left = +(30/180)*np.pi
        rotate_right = R.from_rotvec(rotation_degree_right*rotation_axis)
        rotate_left = R.from_rotvec(rotation_degree_left*rotation_axis)
        self.backbone_down = rotate_right.apply( self.backbone )
        self.backbone_up = rotate_left.apply( self.backbone ) + np.array( [0, self.backbone_down[-1][1], self.backbone_down[-1][2]] ) + np.array( [0,0,3.8] )
        self.trans_x = self.bond
        self.trans_y = self.bond*(3**0.5)
        self.trans_z = self.backbone_up[27][2] - self.backbone[0][2] + self.bond_z
        self.trans_x_num = 10
        self.trans_y_num = 20
        self.trans_z_num = 1
        self.box = []
        
    def OUTPUT_2BEADS_TYPE(self, filename):
        out = open(str(filename),"w")
        print(self.trans_x_num*self.trans_y_num*self.trans_z_num*4*len(self.backbone), file=out)
        print("Required", file=out)
        for j in range(self.trans_y_num):
            for i in range(self.trans_x_num):
                for k in range(self.trans_z_num):
                    for atom in range(len(self.backbone_down)):
                        if atom==0 or atom==27:
                            print( "CH3", self.backbone_down[atom][0] + i*self.trans_x, 
                                          self.backbone_down[atom][1] + j*self.trans_y,  
                                          self.backbone_down[atom][2] + k*self.trans_z, file=out)
                        else:
                            print( "CH2", self.backbone_down[atom][0] + i*self.trans_x, 
                                          self.backbone_down[atom][1] + j*self.trans_y,  
                                          self.backbone_down[atom][2] + k*self.trans_z, file=out)
                    for atom in range(len(self.backbone_up)):
                        if atom==0 or atom==27:
                            print( "CH3", self.backbone_up[atom][0] + i*self.trans_x, 
                                          self.backbone_up[atom][1] + j*self.trans_y,  
                                          self.backbone_up[atom][2] + k*self.trans_z, file=out)
                        else:
                            print( "CH2", self.backbone_up[atom][0] + i*self.trans_x, 
                                          self.backbone_up[atom][1] + j*self.trans_y,  
                                          self.backbone_up[atom][2] + k*self.trans_z, file=out)

                    for atom in range(len(self.backbone_down)):
                        if atom==0 or atom==27:
                            print( "CH3", self.backbone_down[atom][0] + i*self.trans_x + 0.5*self.bond, 
                                          self.backbone_down[atom][1] + j*self.trans_y+0.5*self.bond*(3**0.5), 
                                          self.backbone_down[atom][2] + k*self.trans_z, file=out)
                        else:
                            print( "CH2", self.backbone_down[atom][0] + i*self.trans_x + 0.5*self.bond, 
                                          self.backbone_down[atom][1] + j*self.trans_y+0.5*self.bond*(3**0.5), 
                                          self.backbone_down[atom][2] + k*self.trans_z, file=out)
                    for atom in range(len(self.backbone_up)):
                        if atom==0 or atom==27:
                            print( "CH3", self.backbone_up[atom][0] + i*self.trans_x + 0.5*self.bond, 
                                          self.backbone_up[atom][1] + j*self.trans_y+0.5*self.bond*(3**0.5), 
                                          self.backbone_up[atom][2] + k*self.trans_z, file=out)
                        else:
                            print( "CH2", self.backbone_up[atom][0] + i*self.trans_x + 0.5*self.bond, 
                                          self.backbone_up[atom][1] + j*self.trans_y+0.5*self.bond*(3**0.5), 
                                          self.backbone_up[atom][2] + k*self.trans_z, file=out)
                        
        self.box = [[ min(self.backbone_up[27][0],self.backbone_down[0][0]), max(self.backbone_up[27][0],self.backbone_down[0][0]) + self.trans_x*(self.trans_x_num-1) + self.bond ],
                    [ min(self.backbone_up[27][1],self.backbone_down[0][1]), max(self.backbone_up[27][1],self.backbone_down[0][1]) + self.trans_y*(self.trans_y_num-1) + self.bond ],
                    [ min(self.backbone_up[27][2],self.backbone_down[0][2]), max(self.backbone_up[27][2],self.backbone_down[0][2]) + self.trans_z*(self.trans_z_num-1) + self.bond_z ]]
        return self.box
    
    def OUTPUT_1BEADS_TYPE(self, filename):
        out = open(str(filename),"w")
        print(self.trans_x_num*self.trans_y_num*self.trans_z_num*4*len(self.backbone), file=out)
        print("Required", file=out)
        for j in range(self.trans_y_num):
            for i in range(self.trans_x_num):
                for k in range(self.trans_z_num):
                    for atom in range(len(self.backbone_down)):
                        print( "CH2", self.backbone_down[atom][0] + i*self.trans_x, 
                                      self.backbone_down[atom][1] + j*self.trans_y,  
                                      self.backbone_down[atom][2] + k*self.trans_z, file=out)
                    for atom in range(len(self.backbone_up)):
                        print( "CH2", self.backbone_up[atom][0] + i*self.trans_x, 
                                      self.backbone_up[atom][1] + j*self.trans_y,  
                                      self.backbone_up[atom][2] + k*self.trans_z, file=out)
                        
                    for atom in range(len(self.backbone_down)):
                        print( "CH2", self.backbone_down[atom][0] + i*self.trans_x + 0.5*self.bond, 
                                      self.backbone_down[atom][1] + j*self.trans_y+0.5*self.bond*(3**0.5), 
                                      self.backbone_down[atom][2] + k*self.trans_z, file=out)
                    for atom in range(len(self.backbone_up)):
                        print( "CH2", self.backbone_up[atom][0] + i*self.trans_x + 0.5*self.bond, 
                                      self.backbone_up[atom][1] + j*self.trans_y+0.5*self.bond*(3**0.5), 
                                      self.backbone_up[atom][2] + k*self.trans_z, file=out)                        
        self.box = [[ min(self.backbone_up[27][0],self.backbone_down[0][0]), max(self.backbone_up[27][0],self.backbone_down[0][0]) + self.trans_x*(self.trans_x_num-1) + self.bond ],
                    [ min(self.backbone_up[27][1],self.backbone_down[0][1]), max(self.backbone_up[27][1],self.backbone_down[0][1]) + self.trans_y*(self.trans_y_num-1) + self.bond ],
                    [ min(self.backbone_up[27][2],self.backbone_down[0][2]), max(self.backbone_up[27][2],self.backbone_down[0][2]) + self.trans_z*(self.trans_z_num-1) + self.bond_z ]]
        return self.box        
    
    def VMD_topo_script_2beads(self,filename):
        print("topo" + ";" + 
              "topo guessangles" + ";" + 
              "topo guessdihedrals" + ";" + 
              "topo guessimpropers" + ";" +
              "set sel [atomselect top \"within 1.6 of type CH3\"]" + ";" +
              "set id [$sel get index]" + ";"+
              "foreach {id1 id2} $id {topo addbond $id1 $id2 -bondtype CH3-CH2}" + ";"+
              "mol reanalyze top" + ";"+
              "topo writelammpsdata " + str(filename)
             )
        
    def VMD_topo_script_1beads(self,filename):
        print("topo" + ";" + 
              "topo guessangles" + ";" + 
              "topo guessdihedrals" + ";" + 
              "topo guessimpropers" + ";" +
              "topo writelammpsdata " + str(filename)
             )
        