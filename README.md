# mkvol
This is a program for producing test data for [BccFccRaycaster](https://github.com/Elima85/bccfccraycaster), a volume renderer for data sampled on body-centered cubic (BCC) and face-centered cubic (FCC) lattices. 

mkvol can produce the following images sampled on a Cartesian cubic (CC), BCC, or FCC lattice:
* A sphere,
* A cube,
* An octahedron,
* [The Marschner-Lobb phantom](http://dl.acm.org/citation.cfm?id=951109),
* [The Linnér-Strand phantom](http://link.springer.com/chapter/10.1007/978-3-642-29843-1_57).

The output image is cubic, and the number of voxels is as close as possible to the specified value, given the constraints imposed by the sampling lattice.

## Building
### Requirements
* GCC

### Build instructions
Call
```bash
make
```
in the directory containing mkvol.c and Makefile. 

## Usage
Call 
```bash
./mkvol
```
for instructions.
