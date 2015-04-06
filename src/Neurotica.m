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
  {"Neurotica`Global`", 
   "Neurotica`Util`",
   "Neurotica`Coordinates`",
   "Neurotica`Mesh`",
   "Neurotica`Registration`",
   "Neurotica`MRImage`",
   "Neurotica`FreeSurfer`",
   "Neurotica`NifTI`"}];
   (*"Neurotica`VisualCortex`",*)
   (*"Neurotica`MEG`",*)
Unprotect["Neurotica`*", "Neurotica`Private`*"];
ClearAll[ "Neurotica`*", "Neurotica`Private`*"];

$NeuroticaMajorVersion::usage = "$NeuroticaMajorVersion yields the major version number of the current Neurotical library.";
$NeuroticaMinorVersion::usage = "$NeuroticaMinorVersion yields the minor version number of the current Neurotical library.";
$NeuroticaVersion::usage = "$NeuroticaVersion yields the Neurotica library version number. This is equivalent to {$NeuroticaMajorVersion, $NeuroticaMinorVersion}.";
$NeuroticaVersionNumber::usage = "$NeuroticaVersionNumber yields a version number, which is guaranteed to increase monotonically with the Neurotica` package version.";

$NeuroticaPath::usage = "$NeuroticaPath yields the path of the Neurotica.m primary file in the Neurotica` package.";

NeuroticaReload::usage = "NeuroticaReload[] yields the Neurotica version number after forcing the re-evaluation of all Neurotica library source code. This clears all Neurotica namespace values, thus any calls to functions such as AddFreeSurferSubjectsDirectory[] must be made again.";

Begin["`Private`"];

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
    "NifTI.m"};
  Get[$NeuroticaPath];
  $NeuroticaVersionNumber];
Protect[NeuroticaReload];

(* #NeuroticaVersion ******************************************************************************)
$NeuroticaMajorVersion = 0;
$NeuroticaMinorVersion = 1;
$NeuroticaVersion := {$NeuroticaMajorVersion, $NeuroticaMinorVersion};
$NeuroticaVersionNumber = N[$NeuroticaMajorVersion + $NeuroticaMinorVersion / 100];
Protect[$NeuroticaMajorVersion, $NeuroticaMinorVersion, $NeuroticaVersion, $NeuroticaVersionNumber];

End[];
EndPackage[];