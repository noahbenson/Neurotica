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

## Installation ################################################################

To install Neurotica, you will need to make sure that Mathematica can find a
copy of it in its library directories. On Linux/Unix systems, the local library
directory is "~/.Mathematica/Applications"; on Mac OSX, it is
"~/Library/Mathematica/Applications"; and in Windows it is
"C:\Users\<username>\AppData\Roaming\Mathematica\Applications". To obtain the
library's code itself, you can download the library from GitHub:
https://github.com/noahbenson/Neurotica. The recommended installation method is
to clone this github repository then to make a symbolic link to Mathematica's
library directory for the file src/Neurotica.m and the directory
src/Neurotica. The following code demonstrates how this installation occurs on
Mac OSX from the Terminal:

    ~$ cd Code
    ~/Code$ git clone https://github.com/noahbenson/Neurotica
    Cloning into 'Neurotica'...
    remote: Counting objects : 726, done.
    remote: Compressing objects : 100% (14/14), done.
    remote: Total 726 (delta 4), reused 0 (delta 0), pack - reused 712
    Receiving objects : 100% (726/726), 422.29 KiB | 0 bytes/s, done.
    Resolving deltas : 100% (427/427), done.
    Checking connectivity ... done.
    ~/Code$ cd Neurotica
    ~/Code/Neurotica$ ls
    LICENSE     README.md   Tutorial.nb src
    ~/Code/Neurotica$ cd src
    ~/Code/Neurotica/src$ 
    Neurotica   Neurotica.m
    ~/Code/Neurotica/src$ cd ~/Library/Mathematica/Applications
    ~/Library/Mathematica/Applications$ ln -s ~/Code/Neurotica/Neurotica.m .
    ~/Library/Mathematica/Applications$ ln -s ~/Code/Neurotica/Neurotica/ .
    ~/Library/Mathematica/Applications$ ls Neurotica/
    Coordinates.m  Global.m       Mesh.m         Registration.m VisualCortex.m
    FreeSurfer.m   MRImage.m      NifTI.m        Util.m
    
    # After this installation is complete, you can include the entire Neurotica library
    # into Mathematica by using the typical <<Neurotica` syntax. Additionally, you can
    # upgrade the library by typing 'git pull' from the ~/Code/Neurotica directory.

### JLink ######################################################################

Mathematica interacts with Java libraries through an interface it calls
JLink. The JLink interface boots up an instance of the Java Virtual Machine
(JVM) whenever the JLink package is loaded (via <<JLink`). By default,
Mathematica only allows the JVM to use a small amount of RAM, and this can cause
Neurotica to encounter problems when interfacing with nben library
(https://github.com/noahbenson/nben) that it uses to perform registrations, due
to the moderate memory requirements of keeping coordinate matrices and meta-data
in memory. There are multiple ways to fix this, listed here in the order of the
Neurotica Author's preference:

1. Edit your init file to include the following lines:

     <<JLink`;
     SetOptions[InstallJava, JVMArguments->"-Xmx2g"];
     SetOptions[ReinstallJava, JVMArguments->"-Xmx2g"];
     ReinstallJava[];

   (Note that this will allow the JVM to take up at most 2 GB of RAM; for a
   different amount you can edit the JVMArguments).  Neurotica provides a
   function that will make this edit for you:
   `NeuroticaFixJLinkMemoryPermanent[amount]` where amount is a string that
   contains the memory allocation, e.g., "2g" for 2 gigabytes. This makes a
   permanent edit to the beginning of your init.m file.

2. In a notebook, load JLink first, setup the memory, then load Neurotica. This
   fix is quick and easy; in the cell of your notebook in which you load Neurotica,
   include these lines, ending with your inclusion of Neurotica:

     <<Jlink`
     SetOptions[InstallJava, JVMArguments->"-Xmx2g"];
     SetOptions[ReinstallJava, JVMArguments->"-Xmx2g"];
     ReinstallJava[];
     <<Neurotica`

3. Reinstall Java yourself, then reload Neurotica. This can be done with the
   following code:

     ReinstallJava[JVMArguments->"Xmx2g"];
     NeuroticaReload[];

   Neurotica provides a function that will perform this fix for you, which is
   identical to the `NeuroticaFixJLinkMemoryPermanent[]` amount except that it
   provides only a temporary fix: `NeuroticaJLinkFixMemory[amount]`,
   e.g. `NeuroticaJLinkFixMemory["2g"]` to allocate JLink a max of 2 GB of
   RAM. The downsides of this method are that it interrupts existing JLink
   connections (if you have any) and causes Neurotica to forget certain cached
   data, the latter of which is not usually a significant problem. As long as
   you run this function early in your initialization, such as immediately after
   loading Neurotica, there shouldn't be any problems.

### FreeSurfer #################################################################

Neurotica needs to be able to find your FreeSurfer Subjects' directory
($SUBJECTS_DIR).  You don't need to have FreeSurfer installed to use Neurotica,
but in order to execute the examples in the tutorial, you'll need to have a
subjects directory containing, at the least, the fsaverage, fsaverage_sym, and
bert subjects. (The subject 'bert' is an example subject that is distributed
with FreeSurfer). If you do not have FreeSurfer or FreeSurfer's default subject
directory, please visit [http://freesurfer.net/](http://freesurfer.net/) to
obtain these.

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
