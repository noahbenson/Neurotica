(* Neurotica.m
 *
 * Core entrance to the Neurotica library. This file primarily just includes all other relevant
 * packages.
 *
 * Copyright (C) 2014-2015 by Noah C. Benson.
 * This file is part of the Neurotica library.
 *
 * This program is free software: you can redistribute it and/or modify it under the terms of the
 * GNU General Public License as published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See
 * the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with this program.
 * If not, see <http://www.gnu.org/licenses/>.
 *)

(**************************************************************************************************)
BeginPackage[
  "Neurotica`", 
  {"JLink`",
   "Neurotica`Global`", 
   "Neurotica`Util`",
   "Neurotica`Coordinates`",
   "Neurotica`Mesh`",
   "Neurotica`Registration`",
   "Neurotica`MRImage`",
   "Neurotica`FreeSurfer`",
   "Neurotica`NifTI`",
   "Neurotica`VisualCortex`"}];
   (*"Neurotica`MEG`",*)
Unprotect["Neurotica`*", "Neurotica`Private`*"];
ClearAll[ "Neurotica`*", "Neurotica`Private`*"];

$NeuroticaMajorVersion::usage = "$NeuroticaMajorVersion yields the major version number of the current Neurotical library.";
$NeuroticaMinorVersion::usage = "$NeuroticaMinorVersion yields the minor version number of the current Neurotical library.";
$NeuroticaVersion::usage = "$NeuroticaVersion yields the Neurotica library version number. This is equivalent to {$NeuroticaMajorVersion, $NeuroticaMinorVersion}.";
$NeuroticaVersionNumber::usage = "$NeuroticaVersionNumber yields a version number, which is guaranteed to increase monotonically with the Neurotica` package version.";

$NeuroticaPath::usage = "$NeuroticaPath yields the path of the Neurotica.m primary file in the Neurotica` package.";

NeuroticaReload::usage = "NeuroticaReload[] yields the Neurotica version number after forcing the re-evaluation of all Neurotica library source code. This clears all Neurotica namespace values, thus any calls to functions such as AddFreeSurferSubjectsDirectory[] must be made again.";

Neurotica::initerr = "Initialization Error: `1`";
Neurotica::initwarn = "Initialization Warning: `1`";
Neurotica::initmsg = "Initialization Message: `1`";

$NeuroticaJLinkDetails::usage = "$NeuroticaJLinkDetails evaluates to a string that explains the Neurotica-JLink status on your system.";

NeuroticaFixJLinkMemory::usage = "NeuroticaFixJLinkMemory[amountString] reinstalls the Java Virtual Machine with a maximum amount of RAM as specified by the amountString. The amountString is the text that would appear after the \"-Xmx\" command-line argument to the java program, and generally should take a form such as 1g or 2g (for 1 or 2 gigabytes of max RAM). This requires both resetting Mathematica's JLink connections and reloading the Neurotica library, which can result in the eventual reloading of some cached data and will interrupt JLink connections maintained by other libraries, which also may need to be reloaded. If you use this solution to the JLink memory probem, it is suggested that you run it early in your notebook's initialization.";
NeuroticaFixJLinkMemoryPermanent::usage = "NeuroticaFixJLinkMemoryPermanent[amountString] adds a few lines to the beginning of your Mathematica Kernel initialization file, creating it for you if it does not yet exist. The code injected should ensure that sufficient RAM is always allocated to the Java Virtual Machine upon Mathematica startup. The maximum amount of RAM for the JVM is specified by the amountString. The amountString is the text that would appear after the \"-Xmx\" command-line argument to the java program, and generally should take a form such as \"1g\" (1 gigabyte), \"2g\" (2 gigabytes), or larger. Subsequent calls to this function should overwrite the old code injected by previous calls rather than appending to the file. The first time this function is called, a backup is made. Note that this function does not reset the currently running JVM; after running this function you will need to either call NeuroticaFixJLinkMemory[amountString] or restart the Mathematica Kernel.";

(* These are included here because they have dependencies from many of Neurotica's modules *)
ImageToCortex::usage = "ImageToCortex[img, mesh] yields a a list of the values in the given MRImage3D object img interpolated to the given mesh vertices.
ImageToCortex[img, sub, hemi] yields the interpolation, averaged across the thickness of the cortex, for the given subject and hemisphere. If the hemisphere is LR or All, the result is given as {ImageToCortex[img, sub, LH], ImageToCortex[img, sub, RH]}.
ImageToCortex[img, sub] is equivalent to ImageToCortex[img, sub, LR].

The additional option Weight may be given to instruct the algorithm to use a particular image as a weights in the interpolation. This must be an MRImage3D, Image3D, 3D array, or 3D SparseArray the same size as img.  If Weights are not None (the default), then the option Indeterminate may also be specified, and its value is used whenever the total weight for a vertex is 0 (default: 0).";
CortexToImage::usage = "CortexToImage[sub, hemi, property] yields an MRImage3D object in which the voxels have been interpolated from the given property values on the cortical surface of the given hemisphere of the given subject; interpolation is done using the VoxelToVertexMap[sub, hemi]. If hemi is LR, then both hemispheres are interpolated.
CortexToImage[sub, property] is equivalent to CortexToImage[sub, LR, property].

The option Weight may be given to specify the property name or values of a weight to use when interpolating. If the Weight is not None, then the option Indeterminate may be used to indicate the value that is filled in when the total weight applied to a voxel is 0 (default: 0).";
ImageToCortex::badarg = "Bad argument given to ImageToCortex: `1`";
CortexToImage::badarg = "Bad argument given to ImageToCortex: `1`";

Begin["`Private`"];

Protect[Neurotica];

(* #$NeuroticaPath ********************************************************************************)
$NeuroticaPath = $InputFileName;
Protect[$NeuroticaPath];

(* #NeuroticaReload *******************************************************************************)
NeuroticaReload[] := With[
  {path = FileNameJoin[Append[Most @ FileNameSplit[$NeuroticaPath], "Neurotica"]]},
  Get[FileNameJoin[{path, #}]]& /@ {
    "Global.m",
    "Util.m",
    "Coordinates.m",
    "Mesh.m",
    "Registration.m",
    "MRImage.m",
    "FreeSurfer.m",
    "NifTI.m",
    "VisualCortex.m"};
  Get[$NeuroticaPath];
  $NeuroticaVersionNumber];
Protect[NeuroticaReload];

(* #NeuroticaVersion ******************************************************************************)
$NeuroticaMajorVersion = 0;
$NeuroticaMinorVersion = 1;
$NeuroticaVersion := {$NeuroticaMajorVersion, $NeuroticaMinorVersion};
$NeuroticaVersionNumber = N[$NeuroticaMajorVersion + $NeuroticaMinorVersion / 100];
Protect[$NeuroticaMajorVersion, $NeuroticaMinorVersion, $NeuroticaVersion, $NeuroticaVersionNumber];

(* #JLink Maintenance *****************************************************************************)
(* We want to make sure that we have enough memory in the JVM to run registration; in my
   experiments, I've found that at 512M is enough for non-intense use... *)
LoadJavaClass["java.lang.Runtime"];
With[
  {maxMem = (java`lang`Runtime`getRuntime[])@maxMemory[]},
  $NeuroticaJLinkDetails = If[maxMem < 450*10^6,
    Message[
      Neurotica::initwarn,
      "Registration may not work due to low JLink memory; for more information, evaluate "
      <> "$NeuroticaJLinkDetails or visit https://github.com/noahbenson/Neurotica and see the "
      <> "section on JLink."];
    StringJoin[
      "Mathematica interacts with Java libraries through an interface it calls JLink. The JLink ",
      "interface boots up an instance of the Java Virtual Machine (JVM) whenever the JLink ",
      "package is loaded (via <<JLink`). By default, Mathematica only allows the JVM to use a ",
      "small amount of RAM, and this can cause Neurotica to encounter problems when interfacing ",
      "with nben library (https://github.com/noahbenson/nben) that it uses to perform ",
      "registrations, due to the moderate memory requirements of keeping coordinate matrices and ",
      "meta-data in memory. There are multiple ways to fix this, listed here in the order of the ",
      "Neurotica Author's preference:\n\n",
      "(1) Edit your init file (", FileNameJoin[{$UserBaseDirectory, "Kernel", "init.m"}], ") ",
      "to include the following lines:\n\n",
      "    <<JLink`;\n",
      "    SetOptions[InstallJava, JVMArguments->\"-Xmx2g\"];\n",
      "    SetOptions[ReinstallJava, JVMArguments->\"-Xmx2g\"];\n",
      "    ReinstallJava[];\n\n",
      "    (Note that this will allow the JVM to take up at most 2 GB of RAM; for a different ",
      "amount you can edit the JVMArguments).\n",
      "    Neurotica provides a function that will make this edit for you: ",
      "NeuroticaFixJLinkMemoryPermanent[amount] where amount is a string that contains the memory ",
      "allocation, e.g., \"2g\" for 2 gigabytes. This makes a permanent edit to the beginning of ",
      "your init.m file.\n\n",
      "(2) In a notebook, load JLink first, setup the memory, then load Neurotica. This fix is ",
      "quick and easy; in the cell of your notebook in which you load Neurotica, include these ",
      "lines, ending with your inclusion of Neurotica:\n\n",
      "    <<Jlink`\n",
      "    SetOptions[InstallJava, JVMArguments->\"-Xmx2g\"];\n",
      "    SetOptions[ReinstallJava, JVMArguments->\"-Xmx2g\"];\n",
      "    ReinstallJava[];\n",
      "    <<Neurotica`\n\n",
      "(3) Reinstall Java yourself, then reload Neurotica. This can be done with the following ",
      "code:\n\n",
      "    ReinstallJava[JVMArguments->\"Xmx2g\"];\n",
      "    NeuroticaReload[];\n\n",
      "    Neurotica provides a function that will perform this fix for you, which is identical ",
      "to the NeuroticaFixJLinkMemoryPermanent[amount] function except that it provides only a ",
      "temporary fix: NeuroticaFixJLinkMemory[amount], e.g. NeuroticaFixJLinkMemory[\"2g\"] to ",
      "allocate JLink a max of 2 GB of RAM. The downsides of this method are that it interrupts ",
      "existing JLink connections (if you have any) and causes Neurotica to forget certain cached ",
      "data, the latter of which is not usually noticeable. As long as you run this function ",
      "early in your initialization, such as immediately after loading Neurotica, there shouldn't ",
      "be any problems."],
    StringJoin[
      "Your JLink configuration has allocated a maximum of ", ToString@N[maxMem/(1024^2)],
      " MB of RAM; this should be sufficient for Neurotica's purposes."]]];
Protect[$NeuroticaJLinkDetails];

(* #NeuroticaFixJLinkMemory ***********************************************************************)
NeuroticaFixJLinkMemory[amount_String:"2g"] := (
  SetOptions[InstallJava, JVMArguments -> "-Xmx"<>amount];
  SetOptions[ReinstallJava, JVMArguments -> "-Xmx"<>amount];
  ReinstallJava[];
  NeuroticaReload[];
  True);
Protect[NeuroticaFixJLinkMemory];

(* #NeuroticaFixJLinkMemoryPermanent **************************************************************)
$NeuroticaFixJLinkMemoryPermanentTags = {
  "(*************NEUROTICA*START****************************)",
  "(*************NEUROTICA*END******************************)"};
NeuroticaFixJLinkMemoryPermanentCode[amount_String] := StringJoin[
  "\n",
  $NeuroticaFixJLinkMemoryPermanentTags[[1]], "\n",
  "(*  This block of code has been added by the Neurotica  *)\n",
  "(*  library in order to ensure that sufficient memory   *)\n",
  "(*  is allocated to the Java Virtual Machine.           *)\n",
  "<<JLink`;\n",
  "SetOptions[InstallJava, JVMArguments->\"-Xmx", amount, "\"];\n",
  "SetOptions[ReinstallJava, JVMArguments->\"-Xmx", amount, "\"];\n",
  "ReinstallJava[];\n",
  $NeuroticaFixJLinkMemoryPermanentTags[[2]], "\n\n"];
NeuroticaFixJLinkMemoryPermanent[amount_String:"2g"] := Catch@With[
  {initFile = FileNameJoin[{$UserBaseDirectory, "Kernel", "init.m"}],
   initDir = FileNameJoin[{$UserBaseDirectory, "Kernel"}]},
  If[!FileExistsQ[initFile],
    (* Just write it *)
    (If[!DirectoryQ[initDir] && $Failed === CreateDirectory[initDir], Throw[$Failed]];
     If[$Failed === Export[initFile, NeuroticaFixJLinkePermanentCode[amount], "Text"],
       $Failed,
       True]),
    (* Otherwise, we need to append to it *)
    With[
      {text = SplitBy[
         Import[initFile, "Lines"],
         MemberQ[$NeuroticaFixJLinkMemoryPermanentTags,#]&]},
      With[
        {start = If[Length[#] == 0, Missing["NotFound"], #[[-1, 1]]]&@Position[
           text,
           $NeuroticaFixJLinkMemoryPermanentTags[[1]]],
         end   = If[Length[#] == 0, Missing["NotFound"], #[[-1, 1]]]&@Position[
           text,
           $NeuroticaFixJLinkMemoryPermanentTags[[2]]]},
        With[
          {valid = !MissingQ[start] && !MissingQ[end] && start+2 == end},
          With[
            {prolog = If[valid, text[[ 1     ;; start-1 ]], text],
             epilog = If[valid, text[[ end+1 ;; All     ]], ""]},
            (* If we don't find the old lines in here, we need to make a backup *)
            If[!valid, 
              Check[
                Replace[
                  CopyFile[initFile, initFile <> ".Neurotica-Backup"],
                  $Failed :> Throw[$Failed]],
                Throw[$Failed]]];
            (* Now, write the edited file *)
            If[StringQ@Export[
                 initFile,
                 StringJoin[prolog, NeuroticaFixJLinkMemoryPermanentCode[amount], epilog],
                 "Text"],
              True,
              $Failed]]]]]]];
Protect[NeuroticaFixJLinkMemoryPermanent];

(* #ImageToCortex *********************************************************************************)
Options[ImageToCortex] = {Weight -> None, Indeterminate -> 0};
ImageToCortex[img_?MRImageQ, mesh_?CorticalMeshQ, opts:OptionsPattern[]] := 0; (* #here *)
ImageToCortex[img_?MRImageQ, sub_?SubjectQ, hemi:LH|RH|LR, opts:OptionsPattern[]] := 0;
ImageToCortex[img_?MRImageQ, sub_?SubjectQ, opts:OptionsPattern[]] := ImageToCortex[
  img, sub, LR, opts];
Protect[ImageToCortex];

(* #CortexToImage *********************************************************************************)
Options[CortexToImage] = {Weight -> None, Indeterminate -> 0};
CortexToImage[sub_?SubjectQ, hemi:LH|RH, prop_, opts:OptionsPattern[]] := 0; (* #here *)
CortexToImage[sub_?SubjectQ, LR, prop_, opts:OptionsPattern[]] := {
  CortexToImage[sub, LH, prop, opts],
  CortexToImage[sub, RH, prop, opts]};
CortexToImage[sub_?SubjectQ, prop_, opts:OptionsPattern[]] := CortexToImage[sub, LR, prop, opts];
Protect[CortexToImage];


End[];
EndPackage[];