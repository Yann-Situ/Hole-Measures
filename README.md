# Yann-topo
-----------------

Software for computing the thickness-breadth balls from a 3D volumetric object represented by polygonal mesh of its surface.


## Usage
This software can read a polygonal surface mesh in [OFF](https://en.wikipedia.org/wiki/OFF_(file_format)) format and computes the thickness-breadth balls.

To compile, create a `build` directory in the git initial folder and run `cmake .. ; make` in it.

To execute this software, just do
```
./main_V filename.off
```
For instance, `./main_V ../data/eight.off`

It will compute TB-balls pairs and store them in `filename_V.tb`.

## Output
The thickness-breadth pairs are a pair of balls associated to each hole of the mesh. For each hole of dimension `dim`, the software computes:
- the thickness ball with center at point `(t_x t_y t_z)` with radius `t`
- the breadth ball with center at point `(b_x b_y b_z)` with radius `b`

In the _.tb_ file, every pair of thickness-breadth balls is outputted in a line as
```
dim t t_x t_y t_z b b_x b_y b_z
```
Note that there is always a thickness ball of dimension 0 that has no breadth ball. It appears as
```
0 t t_x t_y t_z inf
```
