# Neurotica ####################################################################

Neurotica is an open-source Neuroscience library for Mathematica. This library
is currently in development and is intended to replace the now depricated 
MmaSurfer library.

## Contents ####################################################################

This repository contains a single Mathematica package, Neurotica; this package
provides several basic functionalities:
 * **Cortical Surface Meshes** - Neurotica contains several methods for 
   representing the cortical surface mesh and computing over it. The mesh data
   themselves are treated simultaneously as Mathematica BoundaryMeshRegion 
   objects, Graph objects, and support properties and options.
 * **MR Images** - The Neurotica library includes code for handling and
   representing 3D MR images; the interface is based on that for Image3D.
 * **FreeSurfer Interoperability** - The ability to communicate between 
   FreeSurfer and Mathematica is one of the main goals of this library. Import
   and Export functions for most FreeSurfer file formats (including MGH/MGZ as
   well as surface, w, and label files) are supported, and FreeSurfer subjects
   can be queried much like data structures.
 * **NifTI Support** - Neurotica supports NifTI-1, NifTI-2, CifTI, and GifTI
   file format importing.
 * **Visual Cortex** - Neurotica contains a small set of functions for examining
   fMRI data specifically relating to the occipital pole.
 * **Registration** - Neurotica includes its own registration library, which
   is designed for the registration of cortical surface-data to ideal 2D models.


## Notes #######################################################################

This library is currently undergoing many rapid changes; please be patient and 
expect the contents to change over the next couple of months.

The notebook Tutorial.nb, included in the root directory of this repository, is
the best place to start learning the library. It includes many examples and
documents all currently supported public functions.

## License #####################################################################

This README file is part of the Neurotica library.

The Neurotica library is free software: you can redistribute it and/or
modify it under the terms of the GNU General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
