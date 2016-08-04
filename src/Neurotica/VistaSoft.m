(* VistaSoft.m
 *
 * The Neurotica`VistaSoft contains code used for loading and interacting with VistaSoft subjects
 * data.
 *
 * Copyright (C) 2016 by Noah C. Benson.
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
  "Neurotica`VistaSoft`",
  {"Neurotica`Global`",
   "Neurotica`Util`"}];
Unprotect["Neurotica`VistaSoft`*", "Neurotica`VistaSoft`Private`*"];
ClearAll[ "Neurotica`VistaSoft`*", "Neurotica`VistaSoft`Private`*"];

ImportVistaSoftMATFile::usage = "ImportVistaSofrMatFile[filename] yields a structured set of data that approcimately resembles the type read in by Matlab when loading VistaSoft MAT files.";

ImportMatlabDataDirectory::usage = "ImportMatlabDataDirectory[dir] yields an association that forms a tree structure of the directories and matlab files (with their contents) of the given directory dir. The data is loaded lazily as it is requested so load-times should be fast even for very deep directories. There is no particular distinction between files and directories in the loaded association tree save for their loaded contents.";

ImportVistaSoftGrayData::usage = "ImportVistaSoftGrayData[path] yields an association containing the data stored in the \"Gray\" directory given by path.";

VistaSoftPRFModel::usage = "VistaSoftPRFModel[fit, coords] yields an association of PRF model data detived from the given PRF fit data and the given subject coordinate data. If fit and/or coords are filenames fo the PRF result file and the coords.mat file, then these are loaded automatically; otherwise these argument should be the result of loading these data files.";
VistaSoftPRFModel::badarg = "Bad argument given to VistaSoftPRFModel: `1`";

VistaSoftImport::badarg = "Bad argument given to VistaSoft importer `1`: `2`";
VistaSoft::err = "VistaSoft error in `1`: `2`";

Begin["`Private`"];

(*  Some Notes
 *  ================================================================================================
 *   * VistaSoft stores its coords (in Gray/coords.mat /coords field and Volime/coords.mat) in the 
 *     IPR orientation; to convert to RAS, {#[[3]], MAX - #[[2]], MAX - #[[1]]}
 *
 *)

(* #ImportVistaSoftMATFile ************************************************************************)
ImportVistaSoftMATFile[fl_String?FileExistsQ] := With[
  {data = Fold[
     ReplaceAll,
     Quiet[
       Import[fl, {"MAT", "LabeledData"}],
       {Partition::ilsmp, Transpose::nmtx}],
     {HoldPattern[Transpose@Partition[{}, _]] :> {},
      HoldPattern[Partition[{}, _]] :> {},
      HoldPattern[Transpose[{}]] :> {}}]},
  ReplaceRepeated[
    Fold[
      Function@Replace[#1, patt:{__Rule} :> Association[patt], {#2}],
      data,
      Reverse[Range@Depth[data] - 1]],
    {assc_?AssociationQ} :> assc]];
Protect[ImportVistaSoftMATFile];

(* #ImportMatlabDataDirectory *********************************************************************)
ImportMatlabDataDirectory[dir_?DirectoryQ] := With[
  {contents = FileNames["*", {dir}]},
  Association@Table[
    Which[
      StringStartsQ[file, "."], Nothing,
      DirectoryQ[file], With[
        {sym = TemporarySymbol["directory"],
         f = file},
        sym := sym = ImportMatlabDataDirectory[f];
        FileBaseName[file] :> sym],
      StringEndsQ[ToLowerCase[file], ".mat"], With[
        {sym = TemporarySymbol["file"],
         f = file},
        sym := sym = ImportVistaSoftMATFile[f];
        FileBaseName[file] :> sym],
      True, Nothing],
    {file, contents}]];
Protect[ImportMatlabDataDirectory];

(* #ImportVistaSoftGrayData ***********************************************************************
 * Notes:
 *  * ImportVistaSoftMATFile yields an association with the following fields:
 *    {"allLeftNodes", "allLeftEdges", "leftPath", "allRightNodes", "allRightEdges", "rightPath",
 *     "keepLeft", "keepRight", "coords", "nodes", "edges", "leftClassFile", "rightClassFile"}
 *)
$VistaSoftGrayDataFields = {
  "allLeftNodes",  "allLeftEdges",  "leftPath",  "keepLeft",  "leftClassFile",
  "allRightNodes", "allRightEdges", "rightPath", "keepRight", "rightClassFile",
  "coords", "nodes", "edges"};
ImportVistaSoftGrayData[path_?DirectoryQ] := With[
  {coords = Replace[
     Check[ImportVistaSoftMATFile@FileNameJoin[{path, "coords.mat"}], $Failed],
     Except[_?AssociationQ] :> Message[
       VistaSoftImport::badarg,
       ImportVistaSoftGrayData, "coords.mat import failed"]]},
  Null];
Protect[$VistaSoftGrayDataFields, ImportVistaSoftGrayData];

(* #VistaSoftPRFModel *****************************************************************************)
Options[VistaSoftPRFModel] = {Reverse -> True};
VistaSoftPRFModel[fit_, coords_, opts:OptionsPattern[]] := Which[
  StringQ[fit] && FileExistsQ[fit], VistaSoftPRFModel[
    ImportVistaSoftMATFile[fit],
    coords,
    opts],
  StringQ[coords] && FileExistsQ[coords], VistaSoftPRFModel[
    fit,
    ImportVistaSoftMATFile[coords],
    opts],
  AssociationQ[coords] && KeyExistsQ[coords, "coords"], VistaSoftPRFModel[
    fit,
    coords["coords"],
    opts],
  !AssociationQ[fit] || !KeyExistsQ[fit, "model"], (
    Message[VistaSoftPRFModel::badarg, "fit must be an associaiton with the key \"model\""];
    $Failed),
  True, With[
    {yc = Replace[OptionValue[Reverse], {True -> -1, _ -> 1}],
     colnames = {
       "Voxel", "X0", "Y0", "PRFSize",
       "SigmaMinor", "SigmaMajor", "SigmaTheta",
       "Exponent", "VarianceExplained"},
     mdl = fit["model"]},
    Dataset@Map[
      Association@Thread[colnames -> #] &,
      Transpose@List[
        Round[{#[[3]], 257 - #[[2]], 257 - #[[1]]} & /@ Transpose[coords]],
        mdl["x0"][[1]],
        mdl["y0"][[1]] * yc,
        (#3 * #1 / ((1 - #3) + #2))&[
           0.5*(mdl["sigma"]["major"][[1]] + mdl["sigma"]["minor"][[1]]),
           Sqrt[mdl["exponent"][[1]]],
           Unitize@Chop[mdl["exponent"][[1]]]],
        mdl["sigma"]["minor"][[1]],
        mdl["sigma"]["major"][[1]],
        mdl["sigma"]["theta"][[1]],
        mdl["exponent"][[1]],
        Unitize[mdl["rawrss"][[1]]]*(
          1 - Clip[
            mdl["rss"][[1]]/(# + (1 - Unitize[#])) &@mdl["rawrss"][[1]],
            {0, 1},
            {0, 1}])]]]];


End[];
EndPackage[];
