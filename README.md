# A*START: Oil-field-Chemistry Project

Requirements :
```diff
- an up-to-date pythob version (>= python2.0)
- basic library NumPy 
```

- This repository is the guidance written especially for those who take the project "Oil-Field Chemistry". The goal is to reduce the time and efforts one needs to catch up. The content includes the input fille preparation and ALBO computation.

- LAMMPS requires input file for the calculation and there is a certain format for the input file. The input file roughly contains three categories: simulation box, simulation set up, and the simulation constraint. It should be noted that the preparation of LAMMPS input file needs one to create manually. To provide how an input file looks like, here an example file C28_CG is included in this repo.

## Data File

- We choose VMD TopoTool because it is a very common and free MD analysis software. VMD TopoTool requires configuration fiie such as xyz file for topology. In the following it is shown how to create the data file by C28_mono_generator.py. First xyz file for monoclinic crystal of C28 (octacosane) is created by "C28_mono_generator.py". The output The can be seen by VMD.

```diff
$ python
>>> import C28_mono_generator
>>> out = C28 mono generator.C28 mono()
>>> out.OUTPUT 1BEADS TYPE("DATA FILE NAME.xyz")
```

- Second is to create data file by VMD TopoTool.

```diff
$ python
>>> import C28 mono generator
>>> out = C28 mono generator.C28 mono()
>>> out.VMD topo script 1beads("DATA_FILE_NAME")
```

- Copy the output. This is basically the script for VMD to identify, define bonds and create data file. The script has to be modified depending on different configuration. Details can be found here: https://www.ks.uiuc.edu/Training/Tutorials/vmd/tutorial-html/node4.html. Go to the console of VMD and type it.

```diff
vmd > topo; topo guessangles; topo guessdihedrals; topo guessimpropers; topo writelammpsdata DATA_FILE_NAME
```

- For two bead type, that is CH3 and CH2, simply change command -1BEADS to -2BEAD. Now the data file "DATA_FILE_NAME" has been create in the directory of VMD. However, one needs to manually type the corresponding coeffcients, box size and masses themselves.

## ALBO Post Calculation (Q6)

- During the project, ALBO calculation (Q6) is used for the indirect study of crystal growth mechanism. In this section, we will show how to use ALBO on LAMMPS and Python script "Q6_analyzer.py". There are seven columns in the output. They are id, type, molID, x, y, z, and Q6 respectively. To visualize, the last column Q6 can be viewed as color.

```diff
$ python
>>> import Q6_analyzer
>>> data = Q6_analyzer.trj analyzer("Q6.lammpstrj")
>>> q6 = data.ALBO()
>>> out = data.OUTPUT ALBO trj(q6,"Q6.lammpstrj","new Q6.lammpstrj")
```
