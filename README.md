# Hole-Measures
Software for computing the hole measures (thickness and breadth-balls) of topological holes from a 3D volumetric object represented by polygonal mesh of its surface. The theoretical work linked with this repository has been published in a special issue of _Computer Aided Design_: [Computing Geometrical Measures of Topological Holes](https://doi.org/10.1016/j.cad.2023.103563) ([HAL version](https://hal.science/hal-04322031v1)).
![graphical-abstract](https://github.com/Yann-Situ/Hole-Measures/assets/51318484/8060e38f-420c-435c-89f3-f8a7bcf37fda)

-----------------

## Requirement
This project requires the [CGAL](https://www.cgal.org/) (on Linux distributions, `sudo apt-get install libcgal-dev libcgal-qt5-dev` should do the job).
The library [PHAT](https://github.com/blazs/phat) is already in the repository (`include/phat/`). The project also requires [Cmake](https://cmake.org/download/) with version above 2.8.10.

## Usage
This software can read a polygonal surface mesh in [OFF](https://en.wikipedia.org/wiki/OFF_(file_format)) format and computes the hole-balls. Two different approaches are developed:
- on one hand, the *link approach* pretends to compute every thickness and breadth-balls of an object as well as their homological link (i.e. the pairing between thickness and breadth-balls that corresponds to the same hole).
- On the other hand, the *medial axes approach* is a faster method that computes thickness and breadth-balls without their homological link.

To compile, create a `build` directory in the git initial folder and run `cmake .. ; make` in it. It will compile and create two executable programs: `main_link` and `main_medial`.

## Link approach
The `main_link` executable corresponds to the *link approach*. More precisely, it computes the hole-balls and the *TB*-pairs using the Voronoi filtration.

Usage: `main_link object.off [-E] [-o output_file] [-h]`
Compute hole measures of the 3D object *object.off* using the Voronoi filtration.
By default, store only the present hole measures in the *object.link.tb* file using the *tb* format (see [Output section](#output)).

Some parameters can be used:
- `-E`,` --exhaustive` : store every hole measures (not only the present hole measures).
- `-o output_file`     : write the hole measures in *output_file*.
- `-h`, `--help`       : display the usage message.

## Medial axes approach
The `main_medial` executable corresponds to the *medial axes approach*. More precisely, it computes the hole-balls the inner and outer medial axis filtration. Note that it produces arbitrary *TB*-pairs that do not corresponds to the real ones.

Usage: `main_medial object.off [-I] [-O] [-E] [-c] [-o output_file] [-h]`
Compute hole-balls of the 3D object *object.off* using the medial axes filtration.
By default, store only the present hole measures in the *object.medial.tb* file using the *.tb* format (see [Output section](#output)).

Some parameters can be used:
- `-I`, `--in`         : compute only holes measures of the inner medial axis filtration.
- `-O`, `--out`        : compute only holes measures of the outer medial axis filtration.
- `-E`, `--exhaustive` : store every hole measures (not only the present hole measures).
- `-c`, `--critical`   : take into account Delaunay critical points to correct wrong hole-balls.
- `-o output_file`     : write the hole measures in *output_file*.
- `-h`, `--help`       : display the usage message.

## Output
A hole measure is a pair of two hole-balls (the *T*-ball and the *B*-ball) that can be associated to a topological hole of an volumetric object.
The software computes hole measures and store them in a *.tb* file.
In the *.tb* file, every pair of hole-balls is outputted in a line as
```
dim t t_x t_y t_z b b_x b_y b_z
```
Where:
- `dim` is the dimension of the considered hole (0, 1 or 2 for 3D objects);
- the thickness-ball is a ball of radius `t` centered at point `(t_x t_y t_z)`;
- the breadth-ball is a ball of radius `b` centered at point `(b_x b_y b_z)`;

Note that there is always a thickness-ball of dimension 0 that has no corresponding breadth-ball. It appears as
```
0 t t_x t_y t_z inf inf inf inf
```

## Visualization
To visualize the result, user can use the dedicated python script `hole-ball-viewer.py` with the corresponding *.off* and *.tb* files.
Note that this script requires [pyvista](https://pyvista.org/), which can be installed with command `pip install 'pyvista[all]'`.

For instance run the following command line:
```
python3 hole-ball-viewer.py data/example/eight.off data/example/eight.tb
```
You should obtain the following figure:
<img width="804" height="636" alt="Capture d’écran du 2026-07-06 17-44-40" src="https://github.com/user-attachments/assets/fbd0cf96-3e6e-4365-b307-3eb50632fcf3" />

The *Filter* slider allows to discard *TB*-pairs that have a low persistence value (i.e. noise).

It is also possible to visualize the corresponding persistence diagram with the argument `--diagram`.
For example, the following line display the exhaustive *TB*-pairs and the 3D model *eight.off*, as well as the corresponding persistence diagram.
```
python3 hole-ball-viewer.py data/example/eight.off data/example/eight.exhaustive.tb --diagram
```
