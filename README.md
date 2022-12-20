# Hole-Measures
Software for computing the thickness-breadth balls of topological holes from a 3D volumetric object represented by polygonal mesh of its surface.
-----------------

## Requirement
This project requires the [CGAL](https://www.cgal.org/) library with Boost and QT5. The library [PHAT](https://github.com/blazs/phat) is already in the repository (`include/phat/`). The project also requires Cmake with version above 2.8.10.

## Usage
This software can read a polygonal surface mesh in [OFF](https://en.wikipedia.org/wiki/OFF_(file_format)) format and computes the thickness-breadth balls.

To compile, create a `build` directory in the git initial folder and run `cmake .. ; make` in it.

To execute this software, just do
```
./main_voronoi filename.off
```
For instance, `./main_voronoi ../data/eight.off`

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

## Visualization
To visualize the result, the user can open the _.tb_ file and its associated _.off_ file in the GUI `tools/gui/yann-topo-gui.pro`.
`tools/gui/yann-topo-gui.pro` can compiled using `qtcreator`.
