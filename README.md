# FEMpy_3D
Finite Element Method (FEM) Solver in 3D

Hi and welcome to the FEMpy_3D project! FEMpy_3D is a fully object-oriented FEM solver in 3D-space, which can use both the direct stiffness method as well as explicit time integration (Forward-Euler discretisation) to solve the given set of equations. The solver is written in Python and all necessary files are located in the source folder. The folder "Verification_documents" provides pdf-files which compare the results of the solver against the analytical solution. Examples are given in the "Examples" Folder. These examples cover everything of what the solver can do and should provide a good insight to get started. 

# Thank you!
A big thank you goes to the team behind [PyVista](https://github.com/pyvista/pyvista) for their great work, which made the visualisation possible in 3D-Space!
Furthermore, [NumPy](https://github.com/numpy/numpy) provided all necessary algebraic tools to realise the calculations which take place.

# Get started
To install FEMpy_3D, download the source code and unzip the folder. Open a PowerShell window in the extracted folder (shift + rightclick on windows) and enter the following into the PowerShell window:
```bash
pip install .
```
Once the installation was successful, the scripts provided in "Examples" can be run. The "Example_Demonstrator.py" functions mainly as demonstrator for how all implemented element types work together. Its deformed structure is shown below. 
![image](https://github.com/user-attachments/assets/e55496f8-a603-4558-8551-c4af205bbbee)

# Meshing:
Automatic meshing is provided for simple structures. Either a 2D-Mesh or a 3D-Mesh can be constructed automatically using triangle or hexahedral elements. Other than that, more complex meshs can be constructed manually as the "Example_Demonstrator.py" shows. An example for the automatic meshing is shown below:
![image](https://github.com/user-attachments/assets/1002d953-7de6-42f6-b884-76011fc82a7e)
Alternatively, it is also possible to read .k files from LS-DYNA to import a mesh modelled in LS-PrePost:
![image](https://github.com/user-attachments/assets/93b9f8d2-2d40-4260-899d-e86545647654)


