# INTRODUCTION

A set of functions for estimated the mass properties (mass, center of mass and
inertia) for a 3D object, input as a [Wavefront OBJ file](https://en.wikipedia.org/wiki/Wavefront_.obj_file). Developed at the
[Royal Veterinary College](https://www.rvc.ac.uk/) [Structure and Motion Lab](https://www.rvc.ac.uk/research/research-centres-and-facilities/structure-and-motion)
to estimate the inertial properties of extinct dinosaurs, as part of my PhD in evolutionary biomechanics under [Professor John Hutchinson](https://www.rvc.ac.uk/about/our-people/john-hutchinson). Written in [Matlab](https://www.mathworks.com/products/matlab.html) using standard libraries.


# INSTRUCTIONS FOR USE

Clone or download / unzip this repository into a local directory, then follow the
walkthrough below. This was originally written when I was using 3DSmax for mesh editing
but the procedure is essentially the same for any 3D CGI or mesh editing app.
For a more in-depth guide to the modelling process, see the mass modelling for blender.docx
file in this repo.

## 1. TEST .OBJ FILE

The ducky.obj file is provided for you to test things out.

## 2. THINGS TO DO BEFORE YOU START

Check that all normals are the right way out:
	•	In 3DSmax, select your mesh, go to the ‘modify’ tab [little box on the right dialogue, looks like a blue rainbow], then go to the ‘face’ sub-object level [solid red triangle], and check the ‘show normals’ box. If all of the little blue sticks are not facing outwards on the mesh, click the ‘unify’ button in the ‘surface properties’ box until they all face outwards on the mesh).
	•	In Blender, select your me sh, go to ‘edit’ mode (tab key), open the ‘transform’ window (n key) and click the third option in the ‘normals’ box (looks like a cube with one yellow side). If all of the little orange sticks are not facing outwards on the mesh, press CTRL+N to recalculate normal.
	•	Check there are no holes. Use either 3DSmax or blender to patch up any holes you find.

## 3. TO ESTIMATE MASS PROPERTIES FOR A SINGLE .OBJ FILE

- Place the .obj file in the same folder as the mass properties code

- In matlab, run the ‘doMassPropertiesAnalysis.m’ script (either right click on it and select ‘run’, left click on it and press ‘F9’, open it and press ‘F5’ or open it and press the big green ‘run’ arrow in the top menu of the editor).

- Follow the instructions printed in the Command Window. It will ask you for the following:
- If you want to use a premade controller file (see below). If this is the first time you have analysed this object, the answer to this is probably ‘n’
- Analysis name – this is what the results will be saved as.
- Density - the density of the object in kg/m3. 1060 is a standard value for vertebrate muscle.
- Fidelity - a controller for the amount of parts the body will be split into to estimate its mass properties – 50 is a good value for this
- Units = the units of the input file, 1 for metres, 0.01 for centimetres, etc. The output will be scaled to metres based on this.
- the .obj file to analyse as a solid object, which you pick from a numbered list (ignore the cavities list for now).
- the .obj file to use as the basis for the system origin (where zero is). Leave blank if you just want to use the world origin, or chose a file to use. The centroid of this object will then be used to zero your CoM position.

- The code will either then output an error (in which case you need to check your meshes – the figure that appears should help you work out where on the mesh the problem is) or will run to completion

- On completion, two files will be output – a CONTROLLER csv file and a RESULTS csv file. The results file has the estimated mass properties for your body. The controller file has all of the information about the analysis – what body was processed with what density, what units and what fidelity.

- To redo the analysis a second time, simply answer ‘y’ to the question about premade controller files, and then choose your CONTROLLER csv file from the displayed list.

## 4. TO ESTIMATE MASS PROPERTIES FOR A MULTIPLE BODY SYSTEM (SUCH AS A WHOLE BODY MODEL)

- place all .obj files in the same folder as the mass properties code

- There are now two options. Firstly, you can make your own CONTROLLER .csv file by running the ‘makeBlankController.m’ file. This makes a csv file that you can then fill in with the relevant info (see PART 1 above for details. NB, density, units and fidelity values should be entered BELOW their headings, not to the side). An additional feature of multiple body analysis is that some bodies can be treated as cavities within other bodies (lungs, nasal cavities etc). There are two columns at the bottom of the blank controller file, one for listing .obj files that will be treated as solid bodies (the SOLID OBJECTS column) and one for .obj files that will be treated as cavities (the  CAVITIES column). These should be listed underneath the appropriate heading.

- If you have pre-made a controller file in this way, simply run ‘doMassPropertiesAnalysis.m’ and answer ‘y’ to the question about premade controller files, and then choose your CONTROLLER csv file from the displayed list.

- The second option is to follow the onscreen instructions, as in PART 1 above, except selecting multiple solid objects from the list shown, and (if needed) also selecting cavities and a system origin/zero object.

- The code will then run as in PART 1.
