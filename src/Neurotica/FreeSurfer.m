(* FreeSurfer.m
 *
 * Basic utility functions for dealing with FreeSurfer data in Mathematica.
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
  "Neurotica`FreeSurfer`",
  {"Neurotica`Global`",
   "Neurotica`Util`",
   "Neurotica`Mesh`",
   "Neurotica`MRImage`"}];
Unprotect["Neurotica`FreeSurfer`*", "Neurotica`FreeSurfer`Private`*"];
ClearAll[ "Neurotica`FreeSurfer`*", "Neurotica`FreeSurfer`Private`*"];

ImportMGH::usage = "ImportMGH[filename, options...] is identical to Import[filename, \"MGH\", options].";
ImportSurface::usage = "ImportSurface[filename, options...] is identical to Import[filename, \"FreeSurferSurface\", options].";
ImportCurv::usage = "ImportCurv[filename, options...] is identical to Import[filename, \"FreeSurferCurv\", options].";
ImportWeights::usage = "ImportWeights[filename, options...] is identical to Import[filename, \"FreeSurferWeights\", options].";
ImportAnnotation::usage = "ImportAnnotation[filename, options...] is identical to Import[filename, \"FreeSurferAnnotation\", options].";
ImportLabel::usage = "ImportLabel[filename, options...] is identical to Import[filename, \"FreeSurferLabel\", options].";
ImportMGH::badfmt = "Bad format when reading `1`: `2`";
ImportMGH::version = "ImportMGH expects a version 1 file; this file is version `1`";
ImportMGH::wrongtype = "MGH object appears to be a `1` type, but `2` type requested.";
ImportMGH::nofile = "File `1` given to ImportMGH could not be found";
ImportSurface::badfmt = "Bad format when importing surface `1`: `2`";
ImportCurv::badfmt = "Bad format when importing curv `1`: `2`";
ImportWeights::badfmt = "Bad format when importing weights `1`: `2`";
ImportAnnotation::badfmt = "Bad format when importing annotation `1`: `2`";
ImportAnnotation::warning = "Warning from ImportAnnotation: `1`";

ExportMGH::usage = "ExportMGH[filename, data, options...] is equivalent to Export[filename, data, \"MGH\", options].";
ExportMGH::badfmt = "Bad format for argument to ExportMGH: `1`";
ExportSurface::usage = "ExportSurface[filename, data, options...] is equivalent to Export[filename, data, \"FreeSurferSurface\", options].";
ExportWeights::usage = "ExportWeights[filename, data, options...] is equivalent to Export[filename, data, \"FreeSurferWeights\", options].";
ExportCurv::usage = "ExportCurv[filename, data, options...] is equivalent to Export[filename, data, \"FreeSurferCurv\", options].";
ExportCurv::badfmt = "Bad format when exporting curv `1`: `2`";
ExportWeights::badfmt = "Bad format when exporting weights `1`: `2`";
(* Not yet supported:
ExportAnnotation::usage = "ExportAnnotation[filename, data, options...] is equivalent to Export[filename, data, \"FreeSurferAnnotation\", options].";
ExportAnnotation::badfmt = "Bad format when exporting annotation `1`: `2`";
ExportAnnotation::warning = "Warning from ExportAnnotation: `1`";
ExportLabel::usage = "ExportLabel[filename, data, options...] is identical to Export[filename, data, \"FreeSurferLabel\", options].";
ExportLabel::badfmt = "ExportLabel given data not in the form of {vertexID -> {p, {R, A, S}} ...}";
*)

$MGHHeaderSize::usage = "$MGHHeaderSize is the number of bytes in the header of an MGH file.";
$MGHOptionalData::usage = "$MGHOptionalData is a list of the optional pieces of data that may be stored at the end of an MGH file.";

$FreeSurferSubjectsDirectories::usage = "$FreeSurferSubjectsDirectories is a list of potential FreeSurfer subjects directories; it is obtained at runtime from the $SUBJECTS_DIR environment variable as well as by scanning common FreeSurfer paths.";
$FreeSurferSubjects::usage = "$FreeSurferSubjects is an association containing potential FreeSurfer subjects. If you wish to add an arbitrary subject to FreeSurfer, use AddFreeSurferSubject; otherwise, this association will contain, as keys, the subject names (named by the last element of their directory, e.g. 'bert' for $SUBJECTS_DIR/bert) AND the ful directory paths of any subject found in the $FreeSurferSubjectsDirectories list (see AddFreeSurferSubjectsDirectory), and, as values, lazily-loaded FreeSurfer subject objects..";
$FreeSurferHomes::usage = "$FreeSurferHomes is a list of potential home directories for FreeSurfer. Homes may be added with AddFreeSurferHome.";
$FreeSurferColorLUT::usage = "$FreeSurferColorLUT is a dispatch table that will replace either an integer or string volume label and yield the remaining data (i.e., a list of the string or integer label, whichever was not given, and a color.";

RHX::usage = "RHX is a keyword that represents the inverted right hemisphere as used by FreeSurfer.";
Protect[RHX];

AddFreeSurferHome::usage = "AddFreeSurferHome[dir] adds the directory dir to the FreeSurfer home structures.";
RemoveFreeSurferHome::usage = "RemoveFreeSurferHome[dir] removes the directories matching the pattern dir from the FreeSurfer subject directories structures.";
AddFreeSurferSubjectsDirectory::usage = "AddSubjectsDir[dir] adds the directory dir to the FreeSurfer home structures.";
RemoveFreeSurferSubjectsDirectory::usage = "RemoveFreeSurferHome[dir] removes the directories matching the pattern dir from the FreeSurfer subject directories structures.";
AddFreeSurferSubject::usage = "AddFreeSurferSubject[dir] adds the directory dir to the FreeSurfer subject structures.";
RemoveFreeSurferSubject::usage = "RemoveFreeSurferSurbject[dir] removes the directories matching the pattern dir from the FreeSurfer subject structures.";

FreeSurferSubject::usage = "FreeSurferSubject[directory] represents a FreeSurfer subject whose data is stored in the directory named by directory.";
FreeSurferSubject::notfound = "FreeSurferSubject directory `1` does not seem to exist";
FreeSurferSubject::baddata = "FreeSurferSubject directory `1` does not seem to contains FreeSurfer data: `2`";
FreeSurferSubjectQ::usage = "FreeSurferSubjectQ[directory] yields true if and only if directory is a string referring to a directory that contains a FreeSurfer subject.";
FreeSurferSubjectData::usage = "FreeSurferSubjectData[...] is a form that represents the data of a single FreeSurfer subject, as obtained via the FreeSurferSubject function.";

(* The Transforms *)
FreeSurferSubjectLinearTransform::usage = "FreeSurferSubjectLinearTransform[sub, name] yields the linear transform written in the file with the given name in the directory sub/mri/transforms/name.xfm.";
FreeSurferSubjectLinearTransform::nofile = "FreeSurferSubjectLinearTransform could not find transform file: `1`";

(* Volume Functions *)
FreeSurferSubjectSegments::usage = "FreeSurferSubjectSegments[sub] yields an Image3D object for subject sub whose values correspond to segments of the brain volume.";
FreeSurferSubjectSegment::usage = "FreeSurferSubjectSegment[sub, label] yeids a list of indices at which the subject's anatomical volume is labeled with label, which may be an integer or a string, either of which must be found in the FreeSurferColorLUT.";
FreeSurferSubjectSegment::badlbl = "Label given to FreeSurferSubjectSegment not found: `1`";
FreeSurferSubjectBaseImage::usage = "FreeSurferSubjectBaseImage[sub, id] yields the volume for the subject sub's imported base image with the given id; these are the files in the mri/orig/ directory with names like 001.mgz. If no id is given, then 1 is assumed.";
FreeSurferSubjectRawavg::usage = "FreeSurferSubjectRawavg[sub] yields the volume for the subject sub's brain as represented in FreeSurfer's rawavg.mgz file, which contains the volume prior to conformation by mri_convert.";
FreeSurferSubjectOriginalBrain::usage = "FreeSurferSubjectOriginalBrain[sub] yields the volume for the subject sub's brain as represented in FreeSurfer's orig.mgz file, which contains the volume immediately after conformation by mri_convert.";
FreeSurferSubjectNormalizedOriginalBrain::usage = "FreeSurferSubjectNormalizedOriginalBrain[sub] yields the volume for the subject sub's brain as represented in FreeSurfer's orig_nu.mgz file, which contains the volume immediately after conformation by mri_convert and normalization.";
FreeSurferSubjectBrain::usage = "FreeSurferSubjectBrain[sub] yields the volume for the subject sub's normalized brain (after skull stripping).";
FreeSurferSubjectWhiteMatter::usage = "FreeSurferSubjectWhiteMattter[sub] yields the volume for the subject sub's white matter.";
FreeSurferSubjectFilledBrain::usage = "FreeSurferSubjectFilledBrain[sub] yields the volume for the subject sub's brain in which the right hemisphere white matter has values of 127 and the left hemisphere has values of 255.";
FreeSurferSubjectFilledMask::usage = "FreeSurferSubjectHemisphere[sub, LH|RH] yields the volume for the subject sub's right or left hemisphere only (with values of 1 for the white matter and 0 elsewhere).
FreeSurferSubjectFilledMask[sub] yields the volume of the subject's while brain with a 1 where there is white matter and a 0 elsewhere.";
FreeSurferSubjectRibbon::usage = "FreeSurferSubjectRibbon[sub, LH|RH] yields the volume for the subject sub's right or left hemisphere ribbon (ie, the non-white matter).";
FreeSurferSubjectRibbon::notfound = "FreeSurferSubject's ribbon file (?h.ribbon.mgz or ?h.ribbon.mgh) not found.";

(* Surface Functions *)
FreeSurferSubjectOriginalSurface::usage = "FreeSurferSubjectOriginalSurface[sub, hemi] yields the original cortical surface tesselation for subject sub's specified hemishere.";
FreeSurferSubjectPialSurface::usage = "FreeSurferSubjectPialSurface[sub, hemi] yields the pial surface object for subject sub's specified hemishere.";
FreeSurferSubjectWhiteSurface::usage = "FreeSurferSubjectWhiteSurface[sub, hemi] yields the white matter surface object for subject sub's specified hemishere.";
FreeSurferSubjectMiddleSurface::usage = "FreeSurferSubjectMiddleSurface[sub, hemi] yields the surface formed by the mean positions of the vertices in the white matter surface and the pial surface as a surface object for subject sub's specified hemishere.";
FreeSurferSubjectInflatedSurface::usage = "FreeSurferSubjectInflatedSurface[sub, hemi] yields the inflated surface object for subject sub's specified hemishere.";
FreeSurferSubjectSphereSurface::usage = "FreeSurferSubjectSphereSurface[sub, hemi] yields the spherical surface object for subject sub's specified hemishere; note that this is not the registration of the subject to the FSAverage hemisphere. For that, see SubjectRegisteredSurface.";
FreeSurferSubjectRegisteredSurface::usage = "FreeSurferSubjectRegisteredSurface[sub, hemi] yields the subject's cortical surface registered to the spherical fsaverage hemisphere for subject sub's specified hemishere.";
FreeSurferSubjectSymSurface::usage = "FreeSurferSubjectSymSurface[sub, hemi] yields the subject's cortical surface registered to the spherical fsaverage_sym hemisphere for subject sub's specified hemishere.";

(* Surface Overlay Functions *)
FreeSurferSubjectJacobian::usage = "FreeSurferSubjectJacobian[sub, hemi] yields the Jacobian field for the subejct sub's given hemisphere.";
FreeSurferSubjectCurvature::usage = "FreeSurferSubjectCurvature[sub, hemi] yields the curvature for subject sub's given hemisphere.";
FreeSurferSubjectSulci::usage = "FreeSurferSubjectSulci[sub, hemi] yields the sulcal depth for subject sub's given hemisphere.";
FreeSurferSubjectThickness::usage = "FreeSurferSubjectThickness[sub, hemi] yields the thickness of each vertex for subject sub's given hemisphere.";
FreeSurferSubjectVertexArea::usage = "FreeSurferSubjectVertexArea[sub, hemi] yields the area of each vertex for subject sub's given hemisphere (?h.mid.area).";
FreeSurferSubjectVertexAreaPial::usage = "FreeSurferSubjectVertexAreaPial[sub, hemi] yields the pial area of each vertex for subject sub's given hemisphere (?h.area.pial).";
FreeSurferSubjectVertexAreaWhite::usage = "FreeSurferSubjectVertexAreaWhite[sub, hemi] yields the white matter area of each vertex for subject sub's given hemisphere (?h.area.white).";
FreeSurferSubjectVertexVolume::usage = "FreeSurferSubjectVertexVolume[sub, hemi] yields the volume of each vertex for subject sub's given hemisphere";
FreeSurferSubjectParcellation::usage = "FreeSurferSubjectParcellation[sub, hemi] yields the 2009 cortical surface parcellation for subject sub's given hemisphere.";
FreeSurferSubjectParcellation2009::usage = "FreeSurferSubjectParcellation[sub, hemi] yields the 2009 cortical surface parcellation for subject sub's given hemisphere.";
FreeSurferSubjectParcellation2005::usage = "FreeSurferSubjectParcellation[sub, hemi] yields the 2005 cortical surface parcellation for subject sub's given hemisphere.";
FreeSurferSubjectBrodmannLabelList::usage = "FreeSurferSubjectBrodmannLabelList[sub, hemi] yields a list of the names of the available Brodmann area labels for the given subject and hemisphere.";
FreeSurferSubjectBrodmannThresholdedLabels::usage = "FreeSurferSubjectBrodmannLabels[sub, hemi] yields a list of the Brodmann labels for each vertex in the surface of the given FreeSurfer subject's given hemisphere.";
FreeSurferSubjectBrodmannsLabels::usage = "FreeSurferSubjectBrodmannThresholds[sub, hemi] yields an association of each Brodmann area number supported by freesurfer, paired to a list of the probability thresholds for each vertex in the surface of the given FreeSurfer subject's given hemisphere.";
FreeSurferSubjectV1ThresholdedLabel::usage = "FreeSurferSubjectV1ThresholdedLabel[sub, hemi] yields a binary (1/0) mask of the location of V1, according to FreeSurfer, for all vertices in the given FreeSurfer subject's given hemisphere.";
FreeSurferSubjectV1Label::usage = "FreeSurferSubjectV1Label[sub, hemi] yields an overlay of the probability of each vertex belonging to V1, according to FreeSurfer, for all vertices in the given FreeSurfer subject's given hemisphere.";
FreeSurferSubjectV2ThresholdedLabel::usage = "FreeSurferSubjectV1ThresholdedLabel[sub, hemi] yields a binary (1/0) mask of the location of V1, according to FreeSurfer, for all vertices in the given FreeSurfer subject's given hemisphere.";
FreeSurferSubjectV2Label::usage = "FreeSurferSubjectV1Label[sub, hemi] yields an overlay of the probability of each vertex belonging to V1, according to FreeSurfer, for all vertices in the given FreeSurfer subject's given hemisphere.";
FreeSurferSubjectMTThresholdedLabel::usage = "FreeSurferSubjectV1ThresholdedLabel[sub, hemi] yields a binary (1/0) mask of the location of V1, according to FreeSurfer, for all vertices in the given FreeSurfer subject's given hemisphere.";
FreeSurferSubjectMTLabel::usage = "FreeSurferSubjectV1Label[sub, hemi] yields an overlay of the probability of each vertex belonging to V1, according to FreeSurfer, for all vertices in the given FreeSurfer subject's given hemisphere.";

FreeSurferSubjectVertexVoxelMapping::usage = "FreeSurferSubjectVertexVoxelMapping[sub, hemi] yields a list with one element for each vertex in the given hemisphere's surfaces of the given FreeSurfer subject, sub, containing a list of the voxel indices (for the native orientation of the subject's ribbon) for the voxels that are aligned with the given vertex.";
FreeSurferSubjectVoxelVertexMapping::usage = "FreeSurferSubjectVertexVoxelMapping[sub, hemi] yields the inverse of the FreeSurferSubjectVertexVoxelMapping; the result is an Associative array mapping each voxel index in the ribbon to the vertex assigned to that voxel.";

FreeSurfer::nolabel = "The FreeSurferColorLUT.txt file could not be found. This file is normally in your FreeSurfer home directory. Try adding a FreeSurfer home directory (via AddFreeSurferHome[]) and re-evaluating.";
FreeSurfer::notfound = "Data not found: `1`";
Protect[FreeSurfer];

$FSAverage::usage = "$FSAverage is the subject directory for the fsaverage FreeSurfer subject. If you add your freesurfer directory (or it is auto-detected), this will be automatically discovered (see AddFreeSurferHome, AddFreeSurferSubjectsDirectory, and AddFreeSurferSubject).";
FSAverageSubject::usage = "FSAverageSubject[hemi] yields the subject data for the fsaverage subject's given hemisphere hemi.";

$FSAverageSym::usage = "$FSAverageSym is the subject directory for the fsaverage_sym FreeSurfer subject. If you add your freesurfer directory (or it is auto-detected), this will be automatically discovered (see AddFreeSurferHome, AddFreeSurferSubjectsDirectory, and AddFreeSurferSubject).";
FSAverageSymSubject::usage = "FSAverageSymSubject yields the subject data for the fsaverage_sym subject.";

FreeSurferSaveConfiguration::usage = "FreeSurferSaveConfiguration[] yields True if the FreeSurfer configuration (the FreeSurfer home directories and subject directories) are successfully saved as part of the Neurotica permanent data association.";
FreeSurferClearConfiguration::usage = "FreeSurferClearConfiguration[] yields True if the FreeSurfer configuration (the FreeSurfer home directories and subject directories) are successfully removed from the Neurotica permanent data association.";

RibbonToCortex::usage = "RibbonToCortex[sub, hemi, data] converts the 3D image data (must be an MRImage3D, an Image3D, or a 3D array the same size as the subject's ribbon) to a list of vertex data the same size as the subject's given hemisphere's VertexList. The option Method may be passed and must be a function such as Mean or Median that determines how to aggregate competing data.";
RibbonToCortex::badarg = "Bad argument given to RibbonToCortex: `1`";
CortexToRibbon::usage = "CortexToRibbon[sub, hemi, data] converts the given list of cortical surface data (which must be the same size as the subject's given hemisphere's VertexList) to an MRImage3D object by replacing the data in the subject's ribbon with that of the vertex closest to the given ribbon element.";

(**************************************************************************************************)
Begin["`Private`"];

$MGHHeaderSize = 284;
Protect[$MGHHeaderSize];

$MGHOptionalData = {{"TR", "Real32"},
                    {"FlipAngle", "Real32"},
                    {"TE", "Real32"},
                    {"TI", "Real32"},
                    {"FoV", "Real32"}};
Protect[$MGHOptionalData];

(* Translate between types and  *)
$MGHTypesToMMA = {
  0 -> "UnsignedInteger8",
  1 -> "Integer32",
  3 -> "Real32",
  4 -> "Integer16",
  t_ :> (
    Message[ImportMGH::badfmt, "header (type)", "Type ("<>ToString[t]<>") not recognized"];
    Throw[$Failed])};
$MMATypesToMGH = Append[
  Map[Rule[#[[2]],#[[1]]]&, Most[$MGHTypesToMMA]],
  t_ :> (
    Message[ExportMGH::badfmt, "header (type)", "Type ("<>ToString[t]<>") not recognized"];
    Throw[$Failed])];

$BytesPerType = {"Real32" -> 4, "UnsignedInteger8" -> 1, "Integer32" -> 4, "Integer16" -> 2};
Protect[$BytesPerType];

(* These are private; for use in registering FreeSurfer's import file formats. *)
Unprotect[ImportMGHHeader, ImportMGHFrames, ImportMGHFooter, ImportMGHData];
ClearAll[ImportMGHHeader, ImportMGHFrames, ImportMGHFooter, ImportMGHData];
ImportMGHHeader[stream_InputStream, opts___] := "Header" -> Catch[
  Block[
    {$ByteOrdering = 1},
    SetStreamPosition[stream, 0];
    With[
      {version = BinaryRead[stream, "Integer32"],
       sane = TrueQ["SanityChecks" /. Append[{opts}, "SanityChecks" -> True]]},
      If[version != 1, Message[ImportMGH::version, version]];
      With[
        {width = BinaryRead[stream, "Integer32"],
         height = BinaryRead[stream, "Integer32"],
         depth = BinaryRead[stream, "Integer32"],
         nframes = BinaryRead[stream, "Integer32"],
         type = Replace[
           BinaryRead[stream, "Integer32"],
           $MGHTypesToMMA],
         dof = BinaryRead[stream, "Integer32"],
         goodRASFlag = BinaryRead[stream, "Integer16"]},
        With[
          {spacings = If[goodRASFlag == 1, BinaryReadList[stream, "Real32", 3], {1.0, 1.0, 1.0}],
           Xras = If[goodRASFlag == 1, BinaryReadList[stream, "Real32", 3], {-1.0, 0.0, 0.0}],
           Yras = If[goodRASFlag == 1, BinaryReadList[stream, "Real32", 3], {0.0, 0.0, -1.0}],
           Zras = If[goodRASFlag == 1, BinaryReadList[stream, "Real32", 3], {0.0, 1.0, 0.0}],
           Cras = If[goodRASFlag == 1, BinaryReadList[stream, "Real32", 3], {0.0, 0.0, 0.0}]},
          {"MGHFileVersion" -> version,
           "Dimensions" -> {width, height, depth},
           "Frames" -> nframes,
           "ImageBufferType" -> type,
           "DegreesOfFreedom" -> dof,
           "Spacings" -> spacings,
           "VOXToRASMatrix" -> Transpose[
             {Append[Xras, 0.0],
              Append[Yras, 0.0],
              Append[Zras, 0.0],
              Append[Cras, 1.0]}]}]]]]];
ImportMGHFrames[stream_InputStream, opts___] := "Frames" -> Catch[
  Block[
    {$ByteOrdering = 1},
    With[
      {header = Replace[
         "Header",
         Append[
           {opts}, 
           "Header" :> Replace[
             ("Header" /. ImportMGHHeader[stream, opts]),
             $Failed :> Throw[$Failed]]]],
       lopts = {opts}},
      With[
        {dims = "Dimensions" /. header,
         type = "ImageBufferType" /. header,
         nframes = "Frames" /. header},
        With[
          {volsz = Apply[Times, dims[[1;;3]]]},
          If[TrueQ["SanityChecks" /. Append[lopts, "SanityChecks" -> True]],
            (* make sure the number of things we're about to read is reasonable *)
            If[volsz < 1 || volsz > 10^9,
              Message[
                ImportMGH::badfmt,
                "import stream",
                StringJoin[
                  "size of volume (", ToString[volsz], 
                  ") seems unreasonable; run with option \"SanityChecks\" -> False",
                  " to ignore this error"]];
              Throw[$Failed]];
            If[nframes < 1 || nframes > 10^5,
              Message[
                ImportMGH::badfmt,
                "import stream",
                StringJoin[
                  "number of frames (", ToString[nframes], 
                  ") seems unreasonable; run with option \"SanityChecks\" -> False",
                  " to ignore this error"]];
              Throw[$Failed]]];
          SetStreamPosition[stream, $MGHHeaderSize];
          If[Count[dims, Except[1]] == 1,
            Table[
              BinaryReadList[stream, type, volsz],
              {nframes}],
            Table[
              Map[
                Reverse,
                Partition[
                  Partition[
                    BinaryReadList[stream, type, volsz],
                    dims[[1]]],
                  dims[[2]]],
                {0,1}],
              {nframes}]]]]]]];
ImportMGHFooter[stream_InputStream, opts___] := "OptionalData" -> Catch[
  Block[
    {$ByteOrdering = 1,
     lopts = {opts}},
    With[
      {header = Replace[
         "Header",
         Append[
           {opts}, 
           "Header" :> Replace[
             ("Header" /. ImportMGHHeader[stream, opts]),
             $Failed :> Throw[$Failed]]]]},
      With[
        {dims = "Dimensions" /. header,
         type = "ImageBufferType" /. header,
         nframes = Replace["Frames", header]},
        With[
          {volsz = Apply[Times, dims[[1;;3]]]},
          If[TrueQ["SanityChecks" /. Append[lopts, "SanityChecks" -> True]],
            (* make sure the number of things we're about to read is reasonable *)
            If[volsz < 1 || volsz > 10^9,
              Message[
                ImportMGH::badfmt,
                "import stream",
                StringJoin[
                  "size of volume (", ToString[volsz], 
                  ") seems unreasonable; run with option \"SanityChecks\" -> False",
                  " to ignore this error"]];
              Throw[$Failed]];
            If[nframes < 1 || nframes > 10^5,
              Message[
                ImportMGH::badfmt,
                "import stream",
                StringJoin[
                  "number of frames (",
                  ToString[nframes],
                  ") seems unreasonable; run with option \"SanityChecks\" -> False",
                  " to ignore this error"]];
              Throw[$Failed]]];
            SetStreamPosition[
              stream,
              $MGHHeaderSize + (volsz * nframes * (type /. $BytesPerType)) - 1];
            If[BinaryRead[stream, "Integer8"] === EndOfFile, Throw[None]];
            With[
              {data = Fold[
                 Function[
                   With[
                     {name = #2[[1]],
                      type = #2[[2]]},
                     Replace[
                       BinaryRead[stream, type],
                       {EndOfFile :> #1, x_ :> Append[#1, name -> x]}]]],
                 {},
                 $MGHOptionalData]},
              If[Length[data] < Length[$MGHOptionalData], 
                data,
                Replace[
                  BinaryReadList[stream, "UnsignedInteger8"],
                  {EndOfFile :> data,
                   chars_ :> Append[
                     data,
                     "tags" -> Cases[
                       Map[
                         Apply[StringJoin, FromCharacterCode /@ Cases[#, Except[0]]]&,
                         Split[chars, (#1 != 0 && #2 != 0)&]],
                       s_String /; StringLength[s] > 0]]}]]]]]]]];
ImportMGHMetaInformation[stream_InputStream, opts___] := "MetaInformation" -> Catch[
  With[
    {header = Replace[
       "Header",
       Append[
         {opts}, 
         "Header" :> Replace[
           ("Header" /. ImportMGHHeader[stream, opts]),
           $Failed :> Throw[$Failed]]]]},
    With[
      {footer = Replace[
         "OptionalData",
         Append[
           {opts},
           "OptionalData" :> Replace[
             ("OptionalData" /. ImportMGHFooter[stream, opts]),
             $Failed :> {}]]]},
      Join[header, footer]]]];
ImportMGHData[stream_InputStream, opts___] := "Data" -> With[
  {header = "Header" /. Append[{opts}, "Header" :> ImportMGHHeader[stream, opts][[2]]]},
  If[header === $Failed,
    $Failed,
    Replace[
      {ImportMGHFrames[stream, "Header" -> header, opts],
       ImportMGHMetaInformation[stream, "Header" -> header, opts]},
      l_List /; Position[l, $Failed] != {} -> $Failed]]];
MGHInterpret[data_] := With[
  {frames = If[ArrayQ[#, 4], Transpose[#, {4,1,2,3}], #]&["Frames" /. data],
   header = "MetaInformation" /. data},
  With[
    {dims = Dimensions[frames],
     spacings = "Spacings" /. header,
     mtx = "VOXToRASMatrix" /. header},
    If[Count[dims, 1, {1}] == Length[dims] - 1,
      Flatten[frames],
      MRImage3D[
        frames,
        Center -> Most @ Dot[mtx, Append[-Dimensions[frames][[1;;3]]/2, 1]],
        RightDirectionVector -> mtx[[1, 1;;3]],
        AnteriorDirectionVector -> mtx[[2, 1;;3]],
        SuperiorDirectionVector -> mtx[[3, 1;;3]],
        VoxelDimensions -> spacings,
        MetaInformation -> ("MetaInformation" /. data)]]]];
ImportMGHObject[stream_InputStream, opts___] := MGHInterpret["Data" /.ImportMGHData[stream, opts]];
ImportMGH[filename_, opts___] := Import[filename, "MGH", opts];
Protect[ImportMGHHeader, ImportMGHFrames, ImportMGHFooter, ImportMGHMetaInformation, ImportMGHData, 
        MGHInterpret, ImportMGHObject, ImportMGH];
(* Register the importer *)
ImportExport`RegisterImport[
  "MGH",
  {"MetaInformation" :> ImportMGHMetaInformation, 
   "Header" :> ImportMGHHeader,
   "Frames" :> ImportMGHFrames,
   "OptionalData" :> ImportMGHFooter,
   "Data" :> ImportMGHData,
   ImportMGHObject},
  "FunctionChannels" -> {"Streams"},
  "BinaryFormat" -> True];

(* Exporting MGH files ****************************************************************************)
ExportMGH[filename_, data_, opts___] := Block[
  {$ByteOrdering = 1},
  (* Make sure we have all the data we need... *)
  With[
    {datExtr = Which[
         ImageQ[data] && Head[data] === Image3D, ImageData[data],
         MRImageQ[data], ImageData[data],
         ArrayQ[data, 3|4], Normal[data],
         VectorQ[data], {{Normal[data]}},
         MatrixQ[data], {{Normal[data]}},
         True, Message[
           ExportMGH::badfmt,
           "Export data for MGH must be a list, MRImage, or Image3D"]],
     meta = Join[
       Replace[
         Cases[{opts}, Rule[(MetaInformation|"MetaInformation"), info_] :> info, {1}],
         {(l_List /; Length[l] > 0) :> First[l],
          _ :> Which[
            ImageQ[data], MetaInformation /. Options[data, MetaInformation],
            MRImageQ[data], {
              "DegreesOfFreedom" -> ("DegreesOfFreedom" /. Options[data, MetaInformation]),
              "Spacings" -> VoxelDimensions[data],
              "VOXToRASMatrix" -> VoxelIndexToCoordinateMatrix[data]},
            True, {}]}],
       {"DegreesOfFreedom" -> 0,
        "Spacings" -> {1.0,1.0,1.0}, 
        "VOXToRASMatrix" -> {{-1., 0., 0., 0.}, {0., 0., -1., 0.}, {0., 1., 0., 0.}}}]},
    If[!ListQ[meta] || Length[meta] == 0,
      (Message[ExportMGH::badfmt, "Invalid MetaInformation"]; $Failed),
      With[
        {outtype = Replace[
           "OutputFormat", 
           Join[Flatten[{opts}], meta, {"OutputFormat" -> "Real32"}]],
         dat = If[ArrayQ[datExtr, 3], {datExtr}, Transpose[datExtr, {2,3,4,1}]],
         fl = OpenWrite[filename, BinaryFormat -> True]},
        With[
          {res = Catch[
             If[fl === $Failed,
               Message[ExportMGH::nofile, filename];
               Throw[$Failed]];
             (* Write header... *)
             BinaryWrite[fl, 1, "Integer32"];
             BinaryWrite[fl, Rest[Dimensions[dat]], "Integer32"];
             BinaryWrite[fl, Length[dat], "Integer32"];
             BinaryWrite[fl, outtype /. $MMATypesToMGH, "Integer32"];
             BinaryWrite[fl, "DegreesOfFreedom" /. meta, "Integer32"];
             BinaryWrite[fl, 1, "Integer16"];
             BinaryWrite[fl, "Spacings" /. meta, "Real32"];
             BinaryWrite[fl, Join@@Transpose@Part["VOXToRASMatrix" /. meta, 1;;3, All], "Real32"];
             BinaryWrite[fl, Table[0, {$MGHHeaderSize - StreamPosition[fl]}], "Integer8"];
             (* write frames... *)
             BinaryWrite[
               fl,
               Switch[Count[Dimensions[dat], Except[1]],
                 1|2, Flatten@dat,
                 _, Flatten[Map[Reverse, dat, {1,2}]]],
               outtype];
             (* Optional data is not currently supported; zeros are written *)
             Scan[BinaryWrite[fl, 0, #[[2]]]&, $MGHOptionalData];
             True]},
          Close[fl];
          If[res === $Failed, $Failed, filename]]]]]];
Protect[ExportMGH];
ImportExport`RegisterExport["MGH", ExportMGH];


(* Importing surface files ************************************************************************)
$SurfaceFileTrianglesID = -2;
ImportNLNLTerminatedString[stream_, opts___] := Catch[
  Apply[
    StringJoin,
    Map[
      FromCharacterCode,
      Reap[
        Module[
          {back2 = BinaryRead[stream, "Integer8"],
           back1 = BinaryRead[stream, "Integer8"]},
          While[back1 != 10 || back2 != 10, 
            Sow[back2];
            With[
              {c = BinaryRead[stream, "Integer8"]},
              If[c === EndOfFile,
                Message[ImportSurface::badfmt, "input stream", "string not terminated"];
                Throw[$Failed]];
              back2 = back1;
              back1 = c;]]];
       ][[2,1]]]]];
ImportSurfaceMetaInformation[stream_, opts___] := "MetaInformation" -> Block[
  {$ByteOrdering = 1},
  SetStreamPosition[stream, 0];
  Catch[
    With[
      {format = BinaryRead[stream, "Integer24"]},
      Which[
        format == $SurfaceFileTrianglesID,
        {"SurfaceFormat" -> "Triangle",
         "CreatorString" -> Replace[
           ImportNLNLTerminatedString[stream, opts],
           $Failed :> Throw[$Failed]],
         "VertexCount" -> BinaryRead[stream, "Integer32"],
         "FaceCount" -> BinaryRead[stream, "Integer32"]},
        True,
        (Message[
           ImportSurface::badfmt,
           "input stream",
           "Input type not recognized: "<>ToString[format]];
         Throw[$Failed])]]]];
ImportSurfaceVertexCoordinates[stream_, opts__] := "VertexCoordinates" -> Block[
  {$ByteOrdering = 1},
  Catch[
    With[
      {header = Replace[
         "MetaInformation",
         Append[
           {opts},
           "MetaInformation" :> ("MetaInformation" /. ImportSurfaceMetaInformation[stream,opts])]]},
      If[header === $Failed, Throw[$Failed]];
      Partition[
        Replace[
          BinaryReadList[stream, "Real32", 3 * ("VertexCount" /. header)],
          {EndOfFile :> (
             Message[
               ImportSurface::badfmt,
               "input stream",
               "EndOfFile reached before end of vertices"];
             Throw[$Failed]),
           x_ /; Not[NumberQ[x]] :> (
             Message[
               ImportSurface::badfmt,
               "input stream",
               "NaN read from file while importing vertices"];
             Throw[$Failed])},
          {1}],
        3]]]];
ImportSurfaceFaceList[stream_, opts__] := "FaceList" -> Block[
  {$ByteOrdering = 1},
  Catch[
    With[
      {header = Replace[
         "MetaInformation",
         Append[
           {opts},
           "MetaInformation" :> ("MetaInformation" /. ImportSurfaceMetaInformation[stream,opts])]]},
      If[header === $Failed, Throw[$Failed]];
      If[!ListQ["MetaInformation" /. {opts}],
        Skip[stream, "Real32", 3 * ("VertexCount" /. header)]];
      With[
        {n = ("VertexCount" /. header)},
        Switch[("SurfaceFormat" /. header),
          "Triangle", Partition[
            Replace[
              BinaryReadList[stream, "Integer32", 3 * ("FaceCount" /. header)],
              {EndOfFile :> (
                 Message[
                   ImportSurface::badfmt,
                   "input stream",
                   "EndOfFile reached before end of vertices"];
                 Throw[$Failed]),
               x_Integer /; 0 <= x < n :> (x + 1),
               x_ :> (
                 Message[
                   ImportSurface::badfmt,
                   "input stream",
                   "out of range value ("<>ToString@x<>") read from file while importing vertices"];
                 Throw[$Failed])},
              {1}],
            3],
          _, (
            Message[
              ImportSurface::badfmt,
              "input stream",
              "Unsupported format: " <> ("SurfaceFormat" /. header)];
            Throw[$Failed])]]]]];
ImportSurfaceData[stream_, opts___] := "Data" -> Catch[
  With[
    {header = Replace[
       "MetaInformation",
       Append[
         {opts}, 
         "MetaInformation" :> ("MetaInformation" /. ImportSurfaceMetaInformation[stream,opts])]]},
    {"MetaInformation" -> header,
     "VertexCoordinates" -> Replace[
       "VertexCoordinates" /. ImportSurfaceVertexCoordinates[
         stream,
         "MetaInformation" -> header,
         opts],
       $Failed :> Throw[$Failed]],
     "FaceList" -> Replace[
       "FaceList" /. ImportSurfaceFaceList[stream, "MetaInformation" -> header, opts],
       $Failed :> Throw[$Failed]]}]];
ImportSurfaceObject[stream_, opts___] := Catch[
  With[
    {dat = Replace[
        "Data" /. ImportSurfaceData[stream, opts],
        {$Failed :> Throw[$Failed]}]},
    CorticalMesh[
      "VertexCoordinates" /. dat,
      "FaceList" /. dat,
      MetaInformation -> ("MetaInformation" /. dat)]]];
ImportSurface[filename_String, opts___] := Import[filename, "FreeSurferSurface", opts];
Protect[ImportSurface, ImportSurfaceMetaInformation, ImportSurfaceVertexCoordinates,
        ImportSurfaceFaces, ImportSurfaceData, ImportSurfaceObject];
(* Register the importer *)
ImportExport`RegisterImport[
  "FreeSurferSurface",
  {"MetaInformation" :> ImportSurfaceMetaInformation,
   "VertexCoordinates" :> ImportSurfaceVertexCoordinates,
   "FaceList" :> ImportSurfaceFaceList,
   "Data" :> ImportSurfaceData,
   ImportSurfaceObject},
  "FunctionChannels" -> {"Streams"},
  "BinaryFormat" -> True];

(* Exporting Surface Format ***********************************************************************)
ExportSurface[filename_, data_, opts___] := Block[
  {$ByteOrdering = 1},
  (* Make sure we have all the data we need... *)
  With[
    {dat = If[CorticalMeshQ[data],
       {VertexCoordinates[data], FaceList[data]},
       (Message[
          ExportSurface::badfmt,
          "Export data must be a cortical mesh object or a list with " <>
           "\"Header\", \"Faces\", and \"VertexList\" rules"];
        $Failed)],
     meta = With[
       {metaOpt = Cases[{opts}, Rule[MetaInformation|"MetaInformation", x_] :> x, {1}]},
       Which[
         Length[metaOpts] > 0, First@metaOpts,
         MetaInformation[data] =!= None && MetaInformation[data] =!= {}, MetaInformation[data],
         True, {}]],
     dfltMeta = {
       "SurfaceFormat" -> "Triangle",
       "CreatorString" -> "Created by Neurotica for Mathematica",
       "VertexCount" -> VertexCount[data],
       "FaceCount" -> FaceCount[data]},
     outtype = Replace["OutputFormat", Append[Flatten[{opts}], "OutputFormat" -> "Real32"]],
     fl = OpenWrite[filename, BinaryFormat -> True]},
    If[fl === $Failed,
      $Failed,
      With[
        {header = meta,
         X = dat[[1]],
         F = dat[[2]],
         createdString = Fold[
           If[#1 == "CreatorString", Replace[#1, #2], #1]&,
           "CreatorString",
           {{opts}, meta, dfltMeta}]},
        With[
          {res = Check[
             (BinaryWrite[fl, $SurfaceFileTrianglesID, "Integer24"];
              BinaryWrite[
                fl,
                Join[ToCharacterCode /@ Characters[createdString], {10, 10}],
                "Integer8"];
              BinaryWrite[fl, Length[X], "Integer32"];
              BinaryWrite[fl, Length[F], "Integer32"];
              BinaryWrite[fl, Flatten[X], "Real32"];
              BinaryWrite[fl, Flatten[F], "Integer32"]),
             $Failed]},
        Close[fl];
        If[res === $Failed, $Failed, filename]]]]]];
Protect[ExportSurface];
ImportExport`RegisterExport["FreeSurferSurface", ExportSurface];


(* Importing Weights Format ***********************************************************************)
ImportWeightsMetaInformation[stream_, opts___] := "MetaInformation" -> Block[
  {$ByteOrdering = 1},
  SetStreamPosition[stream, 0];
  {"Latency" -> BinaryRead[stream, "Integer16"],
   "Count" -> BinaryRead[stream, "Integer24"]}];
ImportWeightsField[stream_, opts___] := "Field" -> Block[
  {$ByteOrdering = 1},
  With[
    {header = Replace[
       "MetaInformation",
        Append[
          {opts}, 
          "MetaInformation" :> ("MetaInformation" /. ImportWeightsMetaInformation[stream, opts])]]},
    If[header === $Failed,
      $Failed,
      Replace[
        Range["Count" /. header],
        Append[
          Map[
            Function[#[[1]] -> #[[2]]],
            BinaryReadList[stream, {"Integer24", "Real32"}, "Count" /. header]],
          _ -> None],
        {1}]]]];
ImportWeightsData[stream_, opts___] := "Data" -> Catch[
  Block[
    {$ByeOrdering = 1},
    With[
      {header = Replace[
         Replace[
           "MetaInformation", 
           {opts,
            "MetaInformation" :> Replace[
              "MetaInformation",
              ImportWeightsMetaInformation[stream, opts]]}],
         $Failed :> Throw[$Failed]]},
      {"MetaInformation" -> header,
       "Field" -> Replace[
         "Field" /. ImportWeightsField[stream, "MetaInformation" -> header, opts],
         $Failed :> Throw[$Failed]]}]]];
ImportWeightsObject[stream_, opts___] := With[
  {data = "Data" /. ImportWeightsData[stream, opts]},
  If[data === $Failed,
    $Failed,
      "Field" /. data]];
ImportWeights[filename_String, opts___] := Import[filename, "FreeSurferWeights", opts];
Protect[ImportWeights, ImportWeightsMetaInformation, ImportWeightsField, ImportWeightsData, ImportWeightsObject];
(* Register the importer *)
ImportExport`RegisterImport[
  "FreeSurferWeights",
  {"MetaInformation" :> ImportWeightsMetaInformation,
   "Field" :> ImportWeightsField,
   "Data" :> ImportWeightsData,
   ImportWeightsObject},
  "FunctionChannels" -> {"Streams"},
  "BinaryFormat" -> True];

(* Exporting weights ******************************************************************************)
ExportWeights[filename_, data_, opts___] := Block[
  {$ByteOrdering = 1},
  (* Make sure we have all the data we need... *)
  With[
    {impdat = If[MatchQ[data, {Rule[_Integer, _?NumericQ]..}],
       data,
       (Message[
          ExportWeights::badfmt,
          "Only 1D lists of numerical data may be exported as weights"];
        $Failed)],
     latency = With[
       {meta = Replace[
          "MetaInformation",
          {opts,
           "MetaInformation" :> Replace[MetaInformation, {opts, _ :> {}}]}]},
       Replace["Latency", Join[{opts}, meta]]]},
    If[impdat === $Failed,
      $Failed,
      With[
        {fl = OpenWrite[filename, BinaryFormat -> True],
         dat = Map[{#[[1]] - 1, #2}&, impdat]},
        With[
          {res = Catch[
             If[fl === $Failed, Message[ExportWeights::nofile, filename]; Throw[$Failed]];
             BinaryWrite[fl, latency, "Integer16"];
             BinaryWrite[fl, Length[impdat], "Integer16"];
             BinaryWrite[fl, dat, {"Integer24", "Real32"}];
             filename]},
          Close[fl];
          res]]]]];
Protect[ExportWeights];
ImportExport`RegisterExport["FreeSurferWeights", ExportWeights];

(* Importing curv files ***************************************************************************)
ImportCurvMetaInformation[stream_, opts___] := "MetaInformation" -> Catch[
  Block[
    {$ByteOrdering = 1},
    With[
      {firstInt = BinaryRead[stream, "Integer24"]},
      Which[
        firstInt == -1, {
          "VertexCount" -> BinaryRead[stream, "Integer32"],
          "FaceCount" -> BinaryRead[stream, "Integer32"],
          "ValuesPerVertex" -> BinaryRead[stream, "Integer32"],
          "CurvType" -> "New"},
        IntegerQ[firstInt] && firstInt > 0, {
          "VertexCount" -> firstInt,
          "FaceCount" -> BinaryRead[stream, "Integer24"],
          "ValuesPerVertex" -> 1,
          "CurvType" -> "Old"},
        True, (
          Message[ImportCurv::badfmt, "input stream", "Could not read header"];
          $Failed)]]]];
ImportCurvField[stream_, opts___] := "Field" -> Catch[
  Block[
    {$ByteOrdering = 1},
    With[
      {header = "MetaInformation" /. Append[
          {opts},
          "MetaInformation" :> Replace[
            "MetaInformation" /. ImportCurvMetaInformation[stream, opts],
            $Failed :> Throw[$Failed]]]},
      Switch["CurvType" /. header,
        "Old", BinaryReadList[stream, "Integer16", "VertexCount" /. header] / 100.0,
        "New", With[
          {k = "ValuesPerVertex" /. header},
          If[k == 1,
            BinaryReadList[stream, "Real32", "VertexCount" /. header],
            BinaryReadList[stream, Table["Real32", {k}], "VertexCount" /. header]]]]]]];
ImportCurvData[stream_, opts___] := "Data" -> Catch[
  With[
    {header = Replace[
       "MetaInformation" /. ImportCurvMetaInformation[stream, opts],
       $Failed :> Throw[$Failed]]},
    {"MetaInformation" -> header,
     "Field" -> Replace[
       "Field" /. ImportCurvField[stream, "MetaInformation" -> header, opts],
       $Failed :> Throw[$Failed]]}]];
ImportCurvObject[stream_, opts___] := Catch[
  With[
    {data = Replace[
       "Data" /. ImportCurvData[stream, opts],
       $Failed :> Throw[$Failed]]},
    "Field" /. data]];
ImportCurv[filename_String, opts___] := Import[filename, "FreeSurferCurv", opts];
Protect[ImportCurv, ImportCurvObject, ImportCurvData, ImportCurvField, ImportCurvMetaInformation];
(* Register the importer *)
ImportExport`RegisterImport[
  "FreeSurferCurv",
  {"MetaInformation" :> ImportCurvMetaInformation,
   "Field" :> ImportCurvField,
   "Data" :> ImportCurvData,
   ImportCurvObject},
  "FunctionChannels" -> {"Streams"},
  "BinaryFormat" -> True];

(* Exporting curv filess **************************************************************************)
ExportCurv[filename_, data_, opts___] := Block[
  {$ByteOrdering = 1},
  (* Make sure we have all the data we need... *)
  With[
    {dat = If[ArrayQ[data, 1|2, NumericQ],
       data,
       (Message[ExportCurv::badfmt, "output must be a field object or list or matrix"];
        $Failed)]},
    With[
      {header = "MetaInformation" /. {
         opts,
         "MetaInformation" :> (MetaInformation /. {
           opts,
           MetaInformation :> {}})}},
      With[
        {curvType = "CurvType" /. Join[{opts}, header, {"CurvType" -> "New"}],
         faceCount = "FaceCount" /. Join[{opts}, header, {"FaceCount" -> -1}]},
        Which[
          dat === $Failed, $Failed,
          curvType == "Old" && ("ValuesPerVertex" /. header) > 1, (
            Message[
              ExportCurv::badfmt, 
              "Cannot output old header type with more than 1 value per vertex"];
            $Failed),
          True, With[
            {fl = OpenWrite[filename, BinaryFormat -> True]},
            With[
              {res = Catch[
                 If[fl === $Failed, Message[ExportCurv::nofile, filename]; Throw[$Failed]];
                 If[curvType == "Old",
                   (BinaryWrite[fl, Length[dat], "Integer24"];
                    BinaryWrite[fl, faceCount, "Integer24"];
                    BinaryWrite[fl, Round[100.0 * Flatten[dat]], "Integer16"]),
                   (BinaryWrite[fl, -1, "Integer24"];
                    BinaryWrite[fl, Length[dat], "Integer32"];
                    BinaryWrite[fl, faceCount, "Integer32"];
                    BinaryWrite[fl, Dimensions[dat] /. {{_} :> 1, {_,c_} :> c}, "Integer32"];
                    BinaryWrite[fl, Flatten[dat], "Real32"])];
                 filename]},
              Close[fl];
              res]]]]]]];
Protect[ExportCurv];
ImportExport`RegisterExport["FreeSurferCurv", ExportCurv];

(* Importing Annotation files *********************************************************************)
ImportAnnotationData[stream_, opts___] := Block[
  {$ByteOrdering = 1},
  Catch[
    With[
      {n = Replace[
         BinaryRead[stream, "Integer32"],
         (x_Integer /; x < 1) :> (
           Message[ImportAnnotation::badfmt, "input stream", "number of vertices < 1"];
           Throw[$Failed])]},
      With[
        {vertices = Map[
           RGBColor[#[[1]]/255.0, #[[2]]/255.0, #[[3]]/255.0, 1.0 - #[[4]]/255.0]&,
           Reverse /@ SortBy[
             BinaryReadList[
               stream, 
               {"Integer32", 
                "UnsignedInteger8", 
                "UnsignedInteger8",
                "UnsignedInteger8",
                "UnsignedInteger8"}, 
               n],
             First
            ][[All, 2;;5]]]},
        With[
          {LUT = If[BinaryRead[stream, "Integer32"] == 0,
             ("LookupTable" /. Append[
                {opts},
                "LookupTable" :> (
                  Message[ImportAnnotation::warning, "Resorting to global label LookupTable"];
                  $FreeSurferColorLUT)]),
             With[
               {maxstruct = (
                  If[BinaryRead[stream, "Integer32"] >= 0,
                    Message[
                      ImportAnnotation::badfmt,
                      "input stream", "old annotation version not supported"];
                    Throw[$Failed]];
                  Replace[
                    BinaryRead[stream, "Integer32"],
                    (x_Integer /; x < 1) :> (
                      Message[ImportAnnotation::badfmt, "input stream", "max label < 1"];
                      Throw[$Failed])]),
                strlen = BinaryRead[stream, "Integer32"]},
               With[
                 {flname = StringJoin@@Most[BinaryReadList[stream, "Character8", strlen]],
                  entries = BinaryRead[stream, "Integer32"]},
                 Dispatch[
                   Table[
                     With[
                       {label = BinaryRead[stream, "Integer32"],
                        name = Apply[
                          StringJoin,
                          Most@BinaryReadList[
                            stream,
                            "Character8",
                            BinaryRead[stream, "Integer32"]]],
                        clr = RGBColor@@{
                          BinaryRead[stream, "UnsignedInteger32"] / 255.0,
                          BinaryRead[stream, "UnsignedInteger32"] / 255.0,
                          BinaryRead[stream, "UnsignedInteger32"] / 255.0,
                          1.0 - BinaryRead[stream, "UnsignedInteger32"] / 255.0}},
                       clr -> {name, label, clr}],
                     {entries}]]]]]},
          If[LUT === $Failed,
            (Message[
               ImportAnnotation::badfmt,
               "input stream", "could not find a color LookupTable"];
             $Failed),
            {"Labels" -> Replace[
               (vertices /. LUT),
               r_RGBColor :> (
                 Message[
                   ImportAnnotation::badfmt, 
                   "input stream", 
                   "some labels not found in annotation lookup table"];
                 $Failed),
               {1}],
             "LookupTable" -> LUT}]]]]]];
ImportAnnotationObject[stream_, opts___] := With[
  {data = ImportAnnotationData[stream, opts]},
  If[data === $Failed,
    $Failed,
    "Labels" /. data]];
ImportAnnotationLookupTable[stream_, opts___] := With[
  {data = ImportAnnotationData[stream, opts]},
  If[data === $Failed,
    $Failed,
    "LookupTable" /. data]];
ImportAnnotation[filename_String, opts___] := Import[filename, "FreeSurferAnnotation", opts];
Protect[ImportAnnotation, ImportAnnotationData, ImportAnnotationLookupTable, ImportAnnotationObject];
(* Register the importer *)
ImportExport`RegisterImport[
  "FreeSurferAnnotation",
  {"Data" :> ImportAnnotationData,
   "LookupTable" :> ImportAnnotationLookupTable,
   ImportAnnotationObject},
  "FunctionChannels" -> {"Streams"},
  "BinaryFormat" -> True];

(* Exporting Annotation filess ********************************************************************)
(*
ExportAnnotation[filename_, data_, opts___] := Block[
  {$ByteOrdering = 1},
  With[
    {dat = Which[
       !ListQ[data], (
         Message[ExportAnnotation::badfmt, "data must be a list of {structID, name, clr}'s"];
         $Failed),
       Count[data, {_Integer, _String, _RGBColor}, {1}] < Length[data], (
         Message[ExportAnnotation::badfmt, "data must be a list of {structID, name, clr}'s"];
         $Failed),
       True, data]},
    With[
      {res = Catch[
         With[
           {fl = OpenWrite[filename, BinaryFormat -> True]},
           BinaryWrite[fl, Length[dat], "Integer32"];
           Do[
             (BinaryWrite[fl, i - 1, "Integer32"];
              BinaryWrite[fl, {255.0,255.0,255.0,-255.0} * (List@@dat[[i,3]]), "UnsignedInteger8"]),
             {i,1,Length@dat}];
           BinaryWrite[fl, 0, "Integer32"];
           

            Close[fl];
            res]]]]]];
Protect[ExportAnnotation];;
ImportExport`RegisterExport["FreeSurferAnnotation", ExportAnnotation];
*)

(* Importing Label Files **************************************************************************)
ImportLabelData[filename_, options___] := "Data" -> Check[
  Flatten[
    Last@Reap[
      Replace[
        Import[filename, "Table"],
        {u_Integer /; u > 0, r_, a_, s_, p_} :> Sow[(u + 1) -> {{r,a,s},p}],
        {1}]]],
  $Failed];
ImportLabel[filename_, options___] := Check[
  With[
    {dat = "Data" /. ImportLabelData[filename, options],
     opts = {options},
     inval = True /. Append[{options}, True -> 1]},
    SparseArray[
      If[inval == "Unthresholded",
        Map[(#[[1]] -> #[[2,2]])&, dat],
        Map[(#[[1]] -> inval)&, dat]],
      {Max /. Append[opts, Max :> Max[dat[[All,1]]]]},
      False /. Append[opts, False -> 0]]],
  $Failed];
ImportExport`RegisterImport[
  "FreeSurferLabel",
  {"Data" :> ImportLabelData,
   ImportLabel}];

(* Exporting Label Files **************************************************************************)
(*
ExportLabel[filename_, data, options___] := Check[
  Export[
    filename,
    Prepend[
      Replace[
        If[ArrayQ[data], data, ImportedData[data]],
        {Rule[id_Integer, {p_, {r_,a_,s}}] :> {id, r, a, s, p},
         _ :> Message[ExportLabel::badfmt]},
        {1}],
      {"#", "Auto-generated", "label", "file", "exported", "from", "Mathematica"}],
    "Table"],
  $Failed];
*)

(* FreeSurfer Directory/filesystem data ***********************************************************)
$FreeSurferHomes = Union[
  Flatten[
    Last[
      Reap[
        Replace[
          {Environment["FREESURFER_HOME"],
           "/usr/local/freesurfer",
           "/opt/freesurfer",
           "/usr/freesurfer",
           "/Applications/freesurfer"},
          s_String /; DirectoryQ[s] && FileExistsQ[s <> "/FreeSurferEnv.sh"] :> Sow[s],
          {1}]]]]];
Protect[$FreeSurferHomes];

$FreeSurferSubjectsDirectories = Union[
  Flatten[
    Last[
      Reap[
        Replace[
          Prepend[
            Map[FileNameJoin[{#, "subjects"}]&, $FreeSurferHomes],
            Environment["SUBJECTS_DIR"]],
          s_String /; FreeSurferSubjectQ[s] :> Sow[s],
          {1}]]]]];
Protect[$FreeSurferSubjectsDirectories];

AutoFindFreeSurferSubjects[] := With[
  {subs = Union @ Flatten @ Map[
     Function @ If[DirectoryQ[#],
       Select[
         FileNames[FileNameJoin[{#, "*"}]],
         FreeSurferSubjectQ]],
     $FreeSurferSubjectsDirectories]},
  Association @ Flatten @ Map[
    Function @ With[
      {sym = TemporarySymbol["fssub"]},
      sym := (sym = FreeSurferSubject[#]);
      {# :> sym, Last@FileNameSplit[#] :> sym}],
    subs]];
Protect[AutoFindFreeSurferSubjects];

$FreeSurferSubjects = AutoFindFreeSurferSubjects[];
Protect[$FreeSurferSubjects];

UpdateSubjectsDirectories[] := (
  Unprotect[$FreeSurferSubjectsDirectories];
  $FreeSurferSubjectsDirectories = First/@Gather[
    Flatten[
      Last[
        Reap[
          Replace[
            Append[
              Join[
                Map[FileNameJoin[{#, "subjects"}]&, $FreeSurferHomes],
                $FreeSurferSubjectsDirectories],
              Environment["SUBJECTS_DIR"]],
            s_String /; FreeSurferSubjectQ[s] :> Sow[s],
            {1}]]]]];
  Protect[$FreeSurferSubjectsDirectories]);
UpdateSubjects[] := With[
  {cur = $FreeSurferSubjects,
   auto = AutoFindFreeSurferSubjects[]},
  Unprotect[$FreeSurferSubjects];
  $FreeSurferSubjects = Fold[
    If[KeyExistsQ[#1, #2[[1]]], #1, Association[#1, #2]]&,
    cur,
    Normal[auto]];
  Protect[$FreeSurferSubjects];
  $FreeSurferSubjects];
Protect[UpdateSubjects];

AddFreeSurferHome[s_String] := (
  Unprotect[$FreeSurferHomes];
  $FreeSurferHomes = Prepend[Complement[$FreeSurferHomes, {s}], s];
  UpdateSubjectsDirectories[];
  UpdateSubjects[];
  Protect[$FreeSurferHomes]);
RemoveFreeSurferHome[s_] := (
  Unprotect[$FreeSurferHomes];
  $FreeSurferHomes = Complement[$FreeSurferHomes, Cases[$FreeSurferHomes, s]];
  Protect[$FreeSurferHomes]);
AddFreeSurferSubjectsDirectory[s_String] := (
  Unprotect[$FreeSurferSubjectsDirectories];
  $FreeSurferSubjectsDirectories = Prepend[Complement[$FreeSurferSubjectsDirectories, {s}], s];
  UpdateSubjects[];
  Protect[$FreeSurferSubjectsDirectories]);
RemoveFreeSurferSubjectsDirectory[s_] := (
  Unprotect[$FreeSurferSubjectsDirectories];
  $FreeSurferSubjectsDirectories = Complement[
    $FreeSurferSubjectsDirectories,
    Select[$FreeSurferSubjectsDirectories, StringMatchQ[#, s]&]];
  Unprotect[$FreeSurferSubjects];
  $FreeSurferSubjects = Select[
    $FreeSurferSubjects,
    Function[
      With[
        {dir = First[StringCases[#, RegularExpression["^(.+)/[^/]"] -> "$1"]]},
        StringMatchQ[dir, s]]]];
  Protect[$Subjets];
  Protect[$FreeSurferSubjectsDir]);
AddFreeSurferSubjects[s_String] := (
  Unprotect[$FreeSurferSubjects];
  $FreeSurferSubjects = Prepend[Complement[$FreeSurferSubjects, {s}], s];
  Protect[$FreeSurferSubjects]);
RemoveFreeSurferSubject[s_] := (
  Unprotect[$FreeSurferSubjects];
  $FreeSurferSubjects = Complement[$FreeSurferSubjects, Cases[$FreeSurferSubjects, s]];
  Protect[$FreeSurferSubjects]);
FreeSurferSubjectQ[s_] := False;
FreeSurferSubjectQ[s_String] := And[
  DirectoryQ[s],
  DirectoryQ[FileNameJoin[{s, "surf"}]],
  DirectoryQ[FileNameJoin[{s, "mri"}]]];
FreeSurferSubjectDirectory[s_?SurfaceQ] := If[SurfaceName[s] =!= s,
  SubjectDirectory[SurfaceName[s]],
  None];

FreeSurferSaveConfiguration[] := And[
  NeuroticaPermanentDatum[
    "FreeSurferHomes", 
    StringJoin[Riffle[$FreeSurferHomes, ":"]]],
  NeuroticaPermanentDatum[
    "FreeSurferSubjectsDirectories", 
    StringJoin[Riffle[$FreeSurferSubjectsDirectories, ":"]]]];
FreeSurferClearConfiguration[] := And[
  NeuroticaPermanentDatum["FreeSurferHomes", None];
  NeuroticaPermanentDatum["FreeSurferSubjectsDirectories", None]];

(* Auto-load the configuration at startup *)
With[
  {homes = NeuroticaPermanentDatum["FreeSurferHomes"],
   subdirs = NeuroticaPermanentDatum["FreeSurferSubjectsDirectories"]},
  If[homes =!= None,
    Scan[
      AddFreeSurferHome,
      Reverse @ StringSplit[homes, ":"]]];
  If[subdirs =!= None,
    Scan[
      AddFreeSurferSubjectsDirectory,
      Reverse @ StringSplit[subdirs, ":"]]]];

Protect[AddFreeSurferHome, RemoveFreeSurferHome, 
        AddFreeSurferSubjectsDirectory, RemoveFreeSurferSubjectsDirectory, 
        AddFreeSurferSubject, RemoveFreeSurferSubject,
        FreeSurferSubjectQ, FreeSurferSaveConfiguration, FreeSurferClearConfiguration];


(* Volume labels from the LUT (for use with aseg.mgz) *********************************************)
$FreeSurferColorLUT := With[
  {dat = Catch[
     Scan[
       Function[
         If[FileExistsQ[# <> "/FreeSurferColorLUT.txt"], 
           Throw[
             Select[
               Import[# <> "/FreeSurferColorLUT.txt", "Table"],
               Length[#] == 6 && NumberQ[#[[1]]] && StringQ[#[[2]]] &]]]],
       $FreeSurferHomes];
     $Failed]},
  If[ListQ[dat],
    Set[
      $FreeSurferColorLUT,
      Dispatch[
        Append[
          Flatten[
            Map[
              Function[
                {#[[1]] -> {#[[2]], RGBColor@@Append[#[[3;;5]] / 255.0, 1 - #[[6]]/255.0]}, 
                 #[[2]] -> {#[[1]], RGBColor@@Append[#[[3;;5]] / 255.0, 1 - #[[6]]/255.0]}}],
              dat]],
          _ -> Indeterminate]]],
    (Message[FreeSurfer::nolabel]; $Failed)]];

(* FreeSurfer's Volume data ***********************************************************************)
FreeSurferSubjectSegments[sub_String] := With[
  {dat = Check[
    If[FileExistsQ[sub <> "/mri/aseg.mgh"],
      Import[FileNameJoin[{sub, "mri", "aseg.mgh"}],  "MGH"],
      Import[FileNameJoin[{sub, "mri", "aseg.mgz"}], {"GZIP", "MGH"}]],
    $Failed]},
   If[dat === $Failed, 
     $Failed, 
     Set[
       FreeSurferSubjectSegments[sub],
       dat]]];
FreeSurferSubjectSegment[sub_String, label:(_String | _Integer)] := Check[
  With[
    {aseg = FreeSurferSubjectSegments[sub],
     lbl = label /. $FreeSurferColorLUT},
    Which[
      lbl === Indeterminate, Message[FreeSurferSubjectSegment::badlbl, label],
      IntegerQ[label], Position[aseg, label, {3}],
      True, Position[aseg, lbl[[1]], {3}]]],
  $Failed];
FreeSurferSubjectWhiteMatter[sub_String] := With[
  {dat = Check[
     If[FileExistsQ[sub <> "/mri/wm.mgh"],
      Import[FileNameJoin[{sub, "mri", "wm.mgh"}],  "MGH"],
      Import[FileNameJoin[{sub, "mri", "wm.mgz"}], {"GZIP", "MGH"}]],
    $Failed]},
   If[dat === $Failed,
     $Failed, 
     Set[
       FreeSurferSubjectWhiteMatter[sub],
       dat]]];
FreeSurferSubjectBaseImage[sub_String, id_Integer] := With[
  {idname = IntegerString[id, 10, 3]},
  With[
    {dat = Check[
       If[FileExistsQ[FileNameJoin[{sub, "mri", "orig", idname <> ".mgh"}]],
         Import[FileNameJoin[{sub, "mri", "orig", idname <> ".mgh"}],  "MGH"],
         Import[FileNameJoin[{sub, "mri", "orig", idname <> ".mgz"}], {"GZIP", "MGH"}]],
       $Failed]},
    If[dat === $Failed, 
      $Failed, 
      Set[
        FreeSurferSubjectBaseImage[sub, id],
        dat]]]];
FreeSurferSubjectBaseImage[sub_String] := FreeSurferSubjectBaseImage[sub, 1];
FreeSurferSubjectRawavg[sub_String] := With[
  {dat = Check[
     If[FileExistsQ[FileNameJoin[{sub, "mri", "rawavg.mgh"}]],
      Import[FileNameJoin[{sub, "mri", "rawavg.mgh"}],  "MGH"],
      Import[FileNameJoin[{sub, "mri", "rawavg.mgz"}], {"GZIP", "MGH"}]],
    $Failed]},
   If[dat === $Failed, 
     $Failed, 
     Set[
       FreeSurferSubjectRawavg[sub], 
       dat]]];
FreeSurferSubjectOriginalBrain[sub_String] := With[
  {dat = Check[
     If[FileExistsQ[FileNameJoin[{sub, "mri", "orig_nu.mgh"}]],
      Import[FileNameJoin[{sub, "mri", "orig_nu.mgh"}],  "MGH"],
      Import[FileNameJoin[{sub, "mri", "orig_nu.mgz"}], {"GZIP", "MGH"}]],
    $Failed]},
   If[dat === $Failed, 
     $Failed, 
     Set[
       FreeSurferSubjectOriginalBrain[sub], 
       dat]]];
FreeSurferSubjectNormalizedOriginalBrain[sub_String] := With[
  {dat = Check[
     If[FileExistsQ[FileNameJoin[{sub, "mri", "T1.mgh"}]],
      Import[FileNameJoin[{sub, "mri", "T1.mgh"}],  "MGH"],
      Import[FileNameJoin[{sub, "mri", "T1.mgz"}], {"GZIP", "MGH"}]],
    $Failed]},
   If[dat === $Failed, 
     $Failed, 
     Set[
       FreeSurferSubjectNormalizedOriginalBrain[sub], 
       dat]]];
FreeSurferSubjectBrain[sub_String] := With[
  {dat = Check[
     If[FileExistsQ[FileNameJoin[{sub, "mri", "brain.mgh"}]],
      Import[FileNameJoin[{sub, "mri", "brain.mgh"}],  "MGH"],
      Import[FileNameJoin[{sub, "mri", "brain.mgz"}], {"GZIP", "MGH"}]],
    $Failed]},
   If[dat === $Failed, 
     $Failed, 
     Set[
       FreeSurferSubjectBrain[sub, opts],
       dat]]];
FreeSurferSubjectFilledBrain[sub_String] := With[
  {dat = Check[
     If[FileExistsQ[FileNameJoin[{sub, "mri", "filled.mgh"}]],
      Import[FileNameJoin[{sub, "mri", "filled.mgh"}],  "MGH"],
      Import[FileNameJoin[{sub, "mri", "filled.mgz"}], {"GZIP", "MGH"}]],
    $Failed]},
   If[dat === $Failed, 
     $Failed, 
     Set[
       FreeSurferSubjectFilledBrain[sub, opt],
       dat]]]; (* 127 -> RH, 255 -> LH *)
FreeSurferSubjectFilledMask[sub_String, hem:LH|RH] := Check[
  With[
    {dat = FreeSurferSubjectFilledBrain[sub],
     v = If[hemi === LH, 255, 127]},
    If[dat =!= $Failed,
      Set[
        FreeSurferSubjectHemisphere[sub, hem],
        ImageApply[If[Flatten[{#}][[1]] == v, 1, 0]&, dat]],
      $Failed]],
  $Failed];
FreeSurferSubjectFilledMask[sub_String, hem:LH|RH] := Check[
  With[
    {dat = FreeSurferSubjectFilledBrain[sub]},
    If[dat =!= $Failed,
      Set[
        FreeSurferSubjectHemisphere[sub, hem],
        ImageApply[If[Flatten[{#}][[1]] > 0, 1, 0]&, dat]],
      $Failed]],
  $Failed];
FreeSurferSubjectRibbon[sub_String, hem:LH|RH] := With[
  {mgh = Check[
     Which[
       FileExistsQ[FileNameJoin[{sub, "mri", ToLowerCase@ToString[hem] <> ".ribbon.mgz"}]], Import[
         FileNameJoin[{sub, "mri", ToLowerCase[ToString[hem]] <> ".ribbon.mgz"}],
         {"GZip", "MGH"}],
       FileExistsQ[FileNameJoin[{sub, "mri", ToLowerCase@ToString[hem] <> ".ribbon.mgh"}]], Import[
         FileNameJoin[{sub, "mri", ToLowerCase@ToString[hem] <> ".ribbon.mgh"}],
         "MGH"],
       True, Message[FreeSurferSubjectRibbon::notfound]],
     $Failed]},
  If[mgh === $Failed, 
    $Failed,
    Set[
      FreeSurferSubjectRibbon[sub, hem], 
      mgh]]];
FreeSurferSubjectRibbon[sub_String] := With[
  {mgh = Check[
     Which[
       FileExistsQ[FileNameJoin[{sub, "mri", "ribbon.mgz"}]], Import[
         FileNameJoin[{sub, "mri", "ribbon.mgz"}],
         {"GZip", "MGH"}],
       FileExistsQ[FileNameJoin[{sub, "mri", "ribbon.mgh"}]], Import[
         FileNameJoin[{sub, "mri", "ribbon.mgh"}],
         "MGH"],
       True, Message[FreeSurferSubjectRibbon::notfound]],
     $Failed]},
  If[mgh === $Failed, 
    $Failed,
    Set[
      FreeSurferSubjectRibbon[sub], 
      mgh]]];

(* FreeSurferSubject Surface Data *****************************************************************)
FreeSurferSubjectSimpleSurface[sub_?DirectoryQ, hemi:(LH|RH|RHX), surf_String] := With[
  {dat = Check[
    With[
      {hemistr = Replace[hemi, {(LH|RHX) -> "lh", RH -> "rh"}],
       dirstr = If[hemi === RHX, 
         FileNameJoin[{sub, "xhemi", "surf"}],
         FileNameJoin[{sub, "surf"}]]},
      Import[FileNameJoin[{dirstr, hemistr <> "." <> surf}], "FreeSurferSurface"]],
    $Failed]},
  If[dat === $Failed, 
    $Failed,
    Set[FreeSurferSubjectSimpleSurface[sub, hemi, surf], dat]]];
FreeSurferSubjectSimpleCurv[sub_String /; DirectoryQ[sub], hemi:(LH|RH|RHX), surf_String] := With[
  {dat = Check[
    With[
      {hemistr = Replace[hemi, {(LH|RHX) -> "lh", RH -> "rh"}],
       dirstr = If[hemi === RHX, 
         FileNameJoin[{sub, "xhemi", "surf"}],
         FileNameJoin[{sub, "surf"}]]},
      Import[FileNameJoin[{dirstr, hemistr <> "." <> surf}], "FreeSurferCurv"]],
    $Failed]},
  If[dat === $Failed, 
    $Failed,
    Set[FreeSurferSubjectSimpleCurv[sub, hemi, surf], dat]]];
FreeSurferSubjectSimpleWeights[sub_String /; DirectoryQ[sub], hemi:(LH|RH|RHX), surf_String] := With[
  {dat = Check[
    With[
      {hemistr = Replace[hemi, {(LH|RHX) -> "lh", RH -> "rh"}],
       dirstr = If[hemi === RHX, 
         FileNameJoin[{sub, "xhemi", "surf"}],
         FileNameJoin[{sub, "surf"}]]},
      Import[FileNameJoin[{dirstr, hemistr <> "." <> surf}], "FreeSurferWeights"]],
    $Failed]},
  If[dat === $Failed, 
    $Failed,
    Set[FreeSurferSubjectSimpleWeights[sub, hemi, surf], dat]]];
FreeSurferSubjectSimpleAnnot[sub_String /; DirectoryQ[sub], hemi:(LH|RH|RHX), label_String] := With[
  {dat = Check[
    With[
      {hemistr = Replace[hemi, {(LH|RHX) -> "lh", RH -> "rh"}],
       dirstr = If[hemi === RHX, 
         FileNameJoin[{sub, "xhemi", "label"}],
         FileNameJoin[{sub, "label"}]]},
      Import[FileNameJoin[{dirstr, hemistr <> "." <> label}], "FreeSurferAnnotation"]],
    $Failed]},
  If[dat === $Failed, 
    $Failed,
    Set[FreeSurferSubjectSimpleAnnot[sub, hemi, surf], dat]]];
FreeSurferSubjectSimpleLabel[sub_String /; DirectoryQ[sub], hemi:(LH|RH|RHX), label_String] := With[
  {dat = Check[
     With[
       {hemistr = Replace[hemi, {(LH|RHX) -> "lh", RH -> "rh"}],
        dirstr = If[hemi === RHX, 
          FileNameJoin[{sub, "xhemi", "label"}],
          FileNameJoin[{sub, "label"}]]},
       Import[FileNameJoin[{dirstr, hemistr <> "." <> label <> ".label"}], 
              "FreeSurferLabel",
              True -> "Unthresholded"]],
     $Failed]},
  If[dat === $Failed, 
    $Failed,
    Set[FreeSurferSubjectSimpleLabel[sub, hemi, surf], dat]]];
FreeSurferSubjectSimpleThresholdedLabel[sub_String /; DirectoryQ[sub], 
                                        hemi:(LH|RH|RHX),
                                        label_String] := With[
  {dat = Check[
     With[
       {hemistr = Replace[hemi, {(LH|RHX) -> "lh", RH -> "rh"}],
        dirstr = If[hemi === RHX, 
          FileNameJoin[{sub, "xhemi", "label"}],
          FileNameJoin[{sub, "label"}]]},
       Import[FileNameJoin[{dirstr, hemistr <> "." <> label <> ".label"}], 
              "FreeSurferLabel",
              True -> 1]],
     $Failed]},
  If[dat === $Failed, 
    $Failed,
    Set[FreeSurferSubjectSimpleThresholdedLabel[sub, hemi, surf], dat]]];

(* FreeSurferSubjectLinearTransform ***************************************************************)
FreeSurferSubjectLinearTransform[sub_String, name_String] := Check[
  With[
    {lines = Import[
       FileNameJoin[{sub, "mri", "transforms", "talairach.xfm"}],
       "String"]},
    ImportString[
      First@StringCases[
        lines,
        "Linear_Transform = \n" ~~ s__ ~~ ";\n" :> s],
      "Table"]],
  $Failed];
Protect[FreeSurferSubjectLinearTransform];

(* FreeSurferSubject specific surfaces ************************************************************)
FreeSurferSubjectOriginalSurface[sub_String, hemi:(LH|RH|RHX)] := Check[
  FreeSurferSubjectSimpleSurface[sub, hemi, "orig"],
  $Failed];
FreeSurferSubjectPialSurface[sub_String, hemi:(LH|RH|RHX)] := Check[
  FreeSurferSubjectSimpleSurface[sub, hemi, "pial"],
  $Failed];
FreeSurferSubjectWhiteSurface[sub_String, hemi:(LH|RH|RHX)] := Check[
  FreeSurferSubjectSimpleSurface[sub, hemi, "white"],
  $Failed];
FreeSurferSubjectInflatedSurface[sub_String, hemi:(LH|RH|RHX)] := Check[
  FreeSurferSubjectSimpleSurface[sub, hemi, "inflated"],
  $Failed];
FreeSurferSubjectSphereSurface[sub_String, hemi:(LH|RH|RHX)] := Check[
  FreeSurferSubjectSimpleSurface[sub, hemi, "sphere"],
  $Failed];
FreeSurferSubjectRegisteredSurface[sub_String, hemi:(LH|RH|RHX)] := Check[
  FreeSurferSubjectSimpleSurface[sub, hemi, "sphere.reg"],
  $Failed];
FreeSurferSubjectSymSurface[sub_String, hemi:(LH|RH)] := Check[
  FreeSurferSubjectSimpleSurface[sub, If[hemi === RH, RHX, hemi], "fsaverage_sym.sphere.reg"],
  $Failed];
Protect[FreeSurferSubjectOriginalSurface,
        FreeSurferSubjectPialSurface,
        FreeSurferSubjectInflatedSurface,
        FreeSurferSubjectSphereSurface,
        FreeSurferSubjectRegisteredSurface,
        FreeSurferSubjectSymSurface, FreeSurferSubjectWhiteSurface];

FreeSurferSubjectMiddleSurface[sub_String, hemi:(LH|RH|RHX)] := With[
  {res = Check[
     Clone[
       FreeSurferSubjectPialSurface[sub, hemi],
       VertexCoordinatesTr -> 0.5 * Plus[
         VertexCoordinatesTr[FreeSurferSubjectPialSurface[sub,hemi]],
         VertexCoordinatesTr[FreeSurferSubjectWhiteSurface[sub,hemi]]]],
     $Failed]},
  If[res === $Failed, 
    res,
    (Unprotect[FreeSurferSubjectMiddleSurface];
     FreeSurferSubjectMiddleSurface[sub,hemi] = res;
     Protect[FreeSurferSubjectMiddleSurface];
     res)]];
Protect[FreeSurferSubjectMiddleSurface];


(* Data that can be merged with surfaces **********************************************************)
FreeSurferSubjectJacobian[sub_String, hemi:(LH|RH|RHX)] := Check[
  FreeSurferSubjectSimpleCurv[sub, hemi, "jacobian_white"],
  $Failed];
FreeSurferSubjectCurvature[sub_String, hemi:(LH|RH|RHX)] := Check[
  FreeSurferSubjectSimpleCurv[sub, hemi, "curv"],
  $Failed];
FreeSurferSubjectSulci[sub_String, hemi:(LH|RH|RHX)] := Check[
  FreeSurferSubjectSimpleCurv[sub, hemi, "sulc"],
  $Failed];
FreeSurferSubjectThickness[sub_String, hemi:(LH|RH|RHX)] := Check[
  FreeSurferSubjectSimpleCurv[sub, hemi, "thickness"],
  $Failed];
FreeSurferSubjectRegisteredCurvature[sub_String, hemi:(LH|RH|RHX)] := Check[
  FreeSurferSubjectSimpleCurv[sub, hemi, "avg_curv"],
  $Failed];
FreeSurferSubjectVertexArea[sub_String, hemi:(LH|RH|RHX)] := Check[
  FreeSurferSubjectSimpleCurv[sub, hemi, "mid.area"],
  $Failed];
FreeSurferSubjectVertexAreaPial[sub_String, hemi:(LH|RH|RHX)] := Check[
  FreeSurferSubjectSimpleCurv[sub, hemi, "area.pial"],
  $Failed];
FreeSurferSubjectVertexAreaWhite[sub_String, hemi:(LH|RH|RHX)] := Check[
  FreeSurferSubjectSimpleCurv[sub, hemi, "area.white"],
  $Failed];
FreeSurferSubjectVertexVolume[sub_String, hemi:(LH|RH|RHX)] := Check[
  FreeSurferSubjectSimpleCurv[sub, hemi, "volume"],
  $Failed];
FreeSurferSubjectParcellation[sub_String, hemi:(LH|RH|RHX)] := Check[
  FreeSurferSubjectSimpleAnnot[sub, hemi, "aparc.a2009s.annot"],
  $Failed];
FreeSurferSubjectParcellation2009[sub_String, hemi:(LH|RH|RHX)] := FreeSurferSubjectParcellation[
  sub, hemi];
FreeSurferSubjectParcellation2005[sub_String, hemi:(LH|RH|RHX)] := Check[
  FreeSurferSubjectSimpleAnnot[sub, hemi, "aparc.annot"],
  $Failed];
FreeSurferSubjectV1Label[sub_String, hemi:(LH|RH|RHX)] := Check[
  FreeSurferSubjectSimpleLabel[sub, hemi, "v1.prob"],
  $Failed];
FreeSurferSubjectV1ThresholdedLabel[sub_String, hemi:(LH|RH|RHX)] := Check[
  With[
    {label = FreeSurferSubjectV1Label[sub, hemi]},
    SparseArray[
      Map[
        #[[1]] -> 1&,
        Select[
          Most@ArrayRules@label,
          #[[2]] > 0.5&]],
      Dimensions[label],
      0]],
  $Failed];
FreeSurferSubjectV2ThresholdedLabel[sub_String, hemi:(LH|RH|RHX)] := Check[
  FreeSurferSubjectSimpleThresholdedLabel[sub, hemi, "V2.thresh"],
  $Failed];
FreeSurferSubjectV2Label[sub_String, hemi:(LH|RH|RHX)] := Check[
  FreeSurferSubjectSimpleLabel[sub, hemi, "V2"],
  $Failed];
FreeSurferSubjectMTThresholdedLabel[sub_String, hemi:(LH|RH|RHX)] := Check[
  FreeSurferSubjectSimpleThresholdedLabel[sub, hemi, "MT.thresh"],
  $Failed];
FreeSurferSubjectMTLabel[sub_String, hemi:(LH|RH|RHX)] := Check[
  FreeSurferSubjectSimpleLabel[sub, hemi, "MT"],
  $Failed];
FreeSurferSubjectBrodmannLabelList[sub_String, hemi:(LH|RH|RHX)] := Check[
  Map[
    StringSplit[Last[FileNameSplit[#]], "."][[2]] &,
    FileNames @ FileNameJoin[
      Join[
        If[hemi === RHX, {sub, "xhemi", "label"}, {sub, "label"}],
        {If[hemi === RH, "rh", "lh"] <> ".BA*.thresh.label"}]]],
  $Failed];
FreeSurferSubjectBrodmannThresholdedLabels[sub_String, hemi:(LH|RH|RHX)] := With[
  {res = Check[
     With[
       {list = FreeSurferSubjectBrodmannLabelList[sub, hemi]},
       If[list === $Failed,
         $Failed,
         Association @ Map[
           Function[# -> FreeSurferSubjectSimpleThresholdedLabel[sub, hemi, # <> ".thresh"]],
           list]]],
     $Failed]},
  If[res === $Failed, res, (FreeSurferSubjectBrodmannLabels[sub, hemi] = res)]];
FreeSurferSubjectBrodmannLabels[sub_String, hemi:(LH|RH|RHX)] := With[
  {res = Check[
     With[
       {list = FreeSurferSubjectBrodmannLabelList[sub, hemi]},
       If[list === $Failed,
         $Failed,
         Association @ Map[
           Function[# -> FreeSurferSubjectSimpleLabel[sub, hemi, #]],
           list]]],
     $Failed]},
  If[res === $Failed, res, (FreeSurferSubjectBrodmannThresholds[sub, hemi] = res)]];

FreeSurferSubjectOP[sub_, hemi:(LH|RH)] := With[
  {V = Check[VertexCoordinates[FreeSurferSubjectInflatedSurface[sub, hemi]], $Failed]},
  If[V === $Failed,
    $Failed,
    Set[
      FreeSurferSubjectOP[sub, hemi],
      First[Ordering[V[[All, 2]], 1]]]]];
FreeSurferSubjectPialOP[sub_, hemi:(LH|RH)] := With[
  {idx = FreeSurferSubjectOP[sub, hemi]},
  If[idx === $Failed,
    $Failed, 
    VertexCoordinates[FreeSurferSubjectPialSurface[sub, hemi]][[idx]]]];
FreeSurferSubjectWhiteOP[sub_, hemi:(LH|RH)] := With[
  {idx = FreeSurferSubjectOP[sub, hemi]},
  If[idx === $Failed,
    $Failed, 
    VertexCoordinates[FreeSurferSubjectWhiteSurface[sub, hemi]][[idx]]]];
FreeSurferSubjectMiddleOP[sub_, hemi:(LH|RH)] := With[
  {idx = FreeSurferSubjectOP[sub, hemi]},
  If[idx === $Failed,
    $Failed, 
    VertexCoordinates[FreeSurferSubjectMiddleSurface[sub, hemi]][[idx]]]];
FreeSurferSubjectInflatedOP[sub_, hemi:(LH|RH)] := With[
  {idx = FreeSurferSubjectOP[sub, hemi]},
  If[idx === $Failed,
    $Failed, 
    VertexCoordinates[FreeSurferSubjectInflatedSurface[sub, hemi]][[idx]]]];
FreeSurferSubjectSphereOP[sub_, hemi:(LH|RH)] := With[
  {idx = FreeSurferSubjectOP[sub, hemi]},
  If[idx === $Failed,
    $Failed,
    VertexCoordinates[FreeSurferSubjectSphereSurface[sub, hemi]][[idx]]]];
FreeSurferSubjectRegisteredOP[sub_, hemi:(LH|RH)] := With[
  {idx = FreeSurferSubjectOP[sub, hemi]},
  If[idx === $Failed,
    $Failed,
    VertexCoordinates[FreeSurferSubjectRegisteredSurface[sub, hemi]][[idx]]]];
FreeSurferSubjectSymOP[sub_, hemi:(LH|RH)] := With[
  {idx = FreeSurferSubjectOP[sub, hemi]},
  If[idx === $Failed,
    $Failed,
    VertexCoordinates[FreeSurferSubjectSymSurface[sub, hemi]][[idx]]]];

Protect[FreeSurferSubjectJacobian, FreeSurferSubjectCurvature,
        FreeSurferSubjectVertexArea, FreeSurferSubjectVertexAreaPial,
        FreeSurferSubjectVertexAreaWhite, FreeSurferSubjectVolume,
        FreeSurferSubjectParcellation,
        FreeSurferSubjectParcellation2009,
        FreeSurferSubjectParcellation2005,
        FreeSurferSubjectV1Label, FreeSurferSubjectV1ThresholdedLabel,
        FreeSurferSubjectV2Label, FreeSurferSubjectV2ThresholdedLabel,
        FreeSurferSubjectMTLabel, FreeSurferSubjectMTThresholdedLabel,
        FreeSurferSubjectBrodmannLabelList,
        FreeSurferSubjectThickness,
        FreeSurferSubjectPialOP, FreeSurferSubjectWhiteOP, 
        FreeSurferSubjectMiddleOP, FreeSurferSubjectInflatedOP, FreeSurferSubjectSphereOP,
        FreeSurferSubjectRegisteredOP, FreeSurferSubjectSymOP];


(* #FreeSurferSubject immutable *******************************************************************)

(* These are meta-data about the various surfaces *)
$FreeSurferSurfaceData = Association[
  {"WhiteSurface" -> <|
     "Pattern" -> "white"|"whitemesh"|"inner"|"innersurface"|"innermesh",
     "MapSurface" -> "SphereSurface"|>,
   "MiddleSurface" -> <|
     "Pattern" -> "mid"|"middle"|"middlemesh"|"midgray"|"midgraysurface"|"midgraymesh",
     "MapSurface" -> "SphereSurface"|>,
   "PialSurface" -> <|
     "Pattern" -> "pial"|"pialmesh"|"outer"|"outersurface"|"outermesh",
     "MapSurface" -> "SphereSurface"|>,
   "InflatedSurface" -> <|
     "Pattern" -> "inflated"|"inflatedmesh",
     "MapSurface" -> "SphereSurface"|>,
   "SphereSurface" -> <|
     "Pattern" -> "sphere"|"spheremesh"|"sphericalsurface"|"sphericalmesh",
     "MapSurface" -> "Sphere"|>,
   "RegisteredSurface" -> <|
     "Pattern" -> ("reg"|"fsaverage"|"fs"|"fsaveragemesh"|"fsaveragesurface"|"registered"
                   |"registeredmesh"),
     "MapSurface" -> "RegisteredSurface"|>,
   "SymRegisteredSurface" -> <|
     "Pattern" -> ("sym"|"symregistered"|"fsaveragesym"|"fsaverage_sym"|"fsaveragesymsurface"
                   |"fsaveragesymmesh"|"symsurface"|"symmesh"|"symregisteredmesh"|"symmetric"
                   |"symmetricsurface"|"symmetricmesh"|"symmetricregistered"
                   |"symmetricregisteredsurface"|"symmetricregisteredmesh"),
     "MapSurface" -> "SymRegisteredSurface"|>,
   "OriginalSurface" -> <|
     "Pattern" -> "orig"|"original"|"originalmesh",
     "MapSurface" -> "SphereSurface"|>}];
Protect[$FreeSurferSurfaceData];

(* And these are meta-data about the various images *)
$FreeSurferImageData = Association[
  {"SegmentImage" -> <|
     "Pattern" -> "seg"|"segs"|"segments"|"segmentsimage"|"segment"|>,
   "BaseImage" -> <|
     "Pattern" -> "base"|"base"|>,
   "OriginalBrainImage" -> <|
     "Pattern" -> Alternatives @@ Flatten[
       Outer[StringJoin, {"orig", "original"}, {"brain",""}, {"image",""}]]|>,
   "NormalizedOriginalBrainImage" -> <|
     "Pattern" -> Alternatives @@ Flatten[
       Outer[StringJoin, {"norm", "normalized"}, {"orig", "original"}, {"brain",""}, {"image",""}]]|>,
   "RawAverageImage" -> <|
     "Pattern" -> Alternatives @@ Flatten[
       Outer[StringJoin, {"raw", ""}, {"average", "avg"}, {"brain", ""}, {"image", ""}]]|>,
   "BrainImage" -> <|
     "Pattern" -> "brain"|>,
   "WhiteMatterImage" -> <|
     "Pattern" -> "wm"|"white"|"whiteimage"|"whitematter"|>,
   "FilledMaskImage" -> <|
     "Pattern" -> "mask"|"maskimage"|"filledmask"|>,
   "FilledBrainImage" -> <|
     "Pattern" -> "filledbrain"|"fill"|>,
   "RibbonImage" -> <|
     "Pattern" -> "ribbon"|"graymatter"|"graymatterimage"|"graymask"|"ribbonmask"|>}];
Protect[$FreeSurferImageData];


DefineImmutable[
  FreeSurferSubject[path_?FreeSurferSubjectQ] :> sub,
  {(* Path gives the subject's directory *)
   Path[sub] -> path,
   (* MetaInformation is just a placeholder *)
   MetaInformation[sub] = {},

   (* The data itself is stored in the subject's association *)
   Association[sub] -> Apply[
     Function[{dat},
       If[$VersionNumber >= 10.0, 
         Association[
           Map[
             If[Head[#] === Rule && ListQ[#[[2]]], #[[1]] -> Association[#[[2]]], #]&,
             dat]],
         (* make a fake association *)
         MimicAssociation[
           Map[
             If[Head[#] === Rule && ListQ[#[[2]]], #[[1]] -> Association[#[[2]]], #]&,
             dat]]]],
     {{"SegmentImage"                 :> FreeSurferSubjectSegments[path],
       "BaseImage"                    :> FreeSurferSubjectBaseImage[path],
       "RawAverageImage"              :> FreeSurferSubjectRawavg[path],
       "OriginalBrainImage"           :> FreeSurferSubjectOriginalBrain[path],
       "NormalizedOriginalBrainImage" :> FreeSurferSubjectNormalizedOriginalBrain[path],
       "BrainImage"                   :> FreeSurferSubjectBrain[path],
       "WhiteMatterImage"             :> FreeSurferWhiteMatter[path],
       "FilledBrainImage"             :> FreeSurferSubjectFilledBrain[path],
       "FilledMaskImage"              -> {LH :> FreeSurferSubjectFilledMask[path, LH],
                                          RH :> FreeSurferSubjectFilledMask[path, RH],
                                          RHX :> FreeSurferSubjectFilledMask[path, RHX]},
       "RibbonImage"                  -> {LH :> FreeSurferSubjectRibbon[path, LH],
                                          RH :> FreeSurferSubjectRibbon[path, RH],
                                          RHX :> FreeSurferSubjectRibbon[path, RHX],
                                          Full :> FreeSurferSubjectRibbon[path]},
       "OriginalSurface"              -> {LH :> FreeSurferSubjectOriginalSurface[path, LH],
                                          RH :> FreeSurferSubjectOriginalSurface[path, RH],
                                          RHX :> FreeSurferSubjectOriginalSurface[path, RHX]},
       "PialSurface"                  -> {LH :> FreeSurferSubjectPialSurface[path, LH],
                                          RH :> FreeSurferSubjectPialSurface[path, RH],
                                          RHX :> FreeSurferSubjectPialSurface[path, RHX]},
       "WhiteSurface"                 -> {LH  :> FreeSurferSubjectWhiteSurface[path, LH],
                                          RH  :> FreeSurferSubjectWhiteSurface[path, RH],
                                          RHX :> FreeSurferSubjectWhiteSurface[path, RHX]},
       "MiddleSurface"                -> {LH  :> FreeSurferSubjectMiddleSurface[path, LH],
                                          RH  :> FreeSurferSubjectMiddleSurface[path, RH],
                                          RHX :> FreeSurferSubjectMiddleSurface[path, RHX]},
       "InflatedSurface"              -> {LH  :> FreeSurferSubjectInflatedSurface[path, LH],
                                          RH  :> FreeSurferSubjectInflatedSurface[path, RH],
                                          RHX :> FreeSurferSubjectInflatedSurface[path, RHX]},
       "SphereSurface"                -> {LH  :> FreeSurferSubjectSphereSurface[path, LH],
                                          RH  :> FreeSurferSubjectSphereSurface[path, RH],
                                          RHX :> FreeSurferSubjectSphereSurface[path, RHX]},
       "RegisteredSurface"            -> {LH  :> FreeSurferSubjectRegisteredSurface[path, LH],
                                          RH  :> FreeSurferSubjectRegisteredSurface[path, RH],
                                          RHX :> FreeSurferSubjectRegisteredSurface[path, RHX]},
       "SymRegisteredSurface"         -> {LH  :> FreeSurferSubjectSymSurface[path, LH],
                                          RH  :> FreeSurferSubjectSymSurface[path, RH],
                                          RHX :> FreeSurferSubjectSymSurface[path, RHX]},
       "Jacobian"                     -> {LH  :> FreeSurferSubjectJacobian[path, LH],
                                          RH  :> FreeSurferSubjectJacobian[path, RH],
                                          RHX :> FreeSurferSubjectJacobian[path, RHX]},
       "Curvature"                    -> {LH  :> FreeSurferSubjectCurvature[path, LH],
                                          RH  :> FreeSurferSubjectCurvature[path, RH],
                                          RHX :> FreeSurferSubjectCurvature[path, RHX]},
       "SulcalDepth"                  -> {LH  :> FreeSurferSubjectSulci[path, LH],
                                          RH  :> FreeSurferSubjectSulci[path, RH],
                                          RHX :> FreeSurferSubjectSulci[path, RHX]},
       "Thickness"                    -> {LH  :> FreeSurferSubjectThickness[path, LH],
                                          RH  :> FreeSurferSubjectThickness[path, RH],
                                          RHX :> FreeSurferSubjectThickness[path, RHX]},
       "VertexArea"                   -> {LH  :> FreeSurferSubjectVertexArea[path, LH],
                                          RH  :> FreeSurferSubjectVertexArea[path, RH],
                                          RHX :> FreeSurferSubjectVertexArea[path, RHX]},
       "VertexAreaPial"               -> {LH  :> FreeSurferSubjectVertexAreaPial[path, LH],
                                          RH  :> FreeSurferSubjectVertexAreaPial[path, RH],
                                          RHX :> FreeSurferSubjectVertexAreaPial[path, RHX]},
       "VertexAreaWhite"              -> {LH  :> FreeSurferSubjectVertexAreaWhite[path, LH],
                                          RH  :> FreeSurferSubjectVertexAreaWhite[path, RH],
                                          RHX :> FreeSurferSubjectVertexAreaWhite[path, RHX]},
       "Parcellation"                 -> {LH  :> FreeSurferSubjectParcellation[path, LH],
                                          RH  :> FreeSurferSubjectParcellation[path, RH],
                                          RHX :> FreeSurferSubjectParcellation[path, RHX]},
       "Parcellation2009"             -> {LH  :> FreeSurferSubjectParcellation2009[path, LH],
                                          RH  :> FreeSurferSubjectParcellation2009[path, RH],
                                          RHX :> FreeSurferSubjectParcellation2009[path, RHX]},
       "Parcellation2005"             -> {LH  :> FreeSurferSubjectParcellation2005[path, LH],
                                          RH  :> FreeSurferSubjectParcellation2005[path, RH],
                                          RHX :> FreeSurferSubjectParcellation2005[path, RHX]},
       "V1Probability"                -> {LH  :> FreeSurferSubjectV1Label[path, LH],
                                          RH  :> FreeSurferSubjectV1Label[path, RH],
                                          RHX :> FreeSurferSubjectV1Label[path, RHX]},
       "V2Probability"                -> {LH  :> FreeSurferSubjectV2Label[path, LH],
                                          RH  :> FreeSurferSubjectV2Label[path, RH],
                                          RHX :> FreeSurferSubjectV2Label[path, RHX]},
       "MTProbability"                -> {LH  :> FreeSurferSubjectMTLabel[path, LH],
                                          RH  :> FreeSurferSubjectMTLabel[path, RH],
                                          RHX :> FreeSurferSubjectMTLabel[path, RHX]},
       "V1Label"                      -> {LH  :> FreeSurferSubjectV1ThresholdedLabel[path, LH],
                                          RH  :> FreeSurferSubjectV1ThresholdedLabel[path, RH],
                                          RHX :> FreeSurferSubjectV1ThresholdedLabel[path, RHX]},
       "V2Label"                      -> {LH  :> FreeSurferSubjectV2ThresholdedLabel[path, LH],
                                          RH  :> FreeSurferSubjectV2ThresholdedLabel[path, RH],
                                          RHX :> FreeSurferSubjectV2ThresholdedLabel[path, RHX]},
       "MTLabel"                      -> {LH  :> FreeSurferSubjectMTThresholdedLabel[path, LH],
                                          RH  :> FreeSurferSubjectMTThresholdedLabel[path, RH],
                                          RHX :> FreeSurferSubjectMTThresholdedLabel[path, RHX]},
       "BrodmannLabels"               -> {LH  :> FreeSurferSubjectBrodmannThresholdedLabels[path, LH],
                                          RH  :> FreeSurferSubjectBrodmannThresholdedLabels[path, RH],
                                          RHX :> FreeSurferSubjectBrodmannThresholdedLabels[path, RHX]},
       "BrodmannProbabilities"        -> {LH  :> FreeSurferSubjectBrodmannLabels[path, LH],
                                          RH  :> FreeSurferSubjectBrodmannLabels[path, RH],
                                          RHX :> FreeSurferSubjectBrodmannLabels[path, RHX]},
       "OccipitalPoleIndex"           -> {LH  :> FreeSurferSubjectOP[path, LH],
                                          RH  :> FreeSurferSubjectOP[path, RH],
                                          RHX :> FreeSurferSubjectOP[path, RHX]},
       "TalairachTransform"           :> FreeSurferSubjectLinearTransform[path, "talairach"]}}],
                                 
   (* Now we make some accessors for this subject *)
   Cortex[sub, name_, hemi:(LH|RH|RHX)] := With[
     {assoc = Association[sub],
      id = If[name === Automatic, "Sphere", ToLowerCase[name]] // Function @ FirstCase[
        Normal @ $FreeSurferSurfaceData,
        (r_Rule /; Or[MatchQ[ToLowerCase[#], ToLowerCase[r[[1]]]],
                      MatchQ[ToLowerCase[#], r[[2]]["Pattern"]]]) :> r[[1]]]},
     With[
       {mapMeshName = $FreeSurferSurfaceData[id]["MapSurface"]},
       With[
         {mesh = assoc[id][hemi] // CorticalMesh[
            #,
            MetaInformation -> Join[
              Options[#, MetaInformation],
              {"CorticalMap" -> {
                 Method -> "Mollweide",
                 Center :> OccipitalPole[sub, mapMeshName, hemi],
                 Radius -> Full},
               "SphericalMesh" :> Cortex[sub, mapMeshName, hemi],
               "Subject" -> sub,
               "SurfaceName" -> id,
               "Hemisphere" -> hemi}]] &},
         SetVertexProperties[
           mesh,
           {"Curvature" :> Quiet@Check[assoc["Curvature"][hemi], $Failed],
            "SulcalDepth" :> Quiet@Check[assoc["SulcalDepth"][hemi], $Failed],
            "Thickness" :> Quiet@Check[assoc["Thickness"][hemi], $Failed],
            "VertexArea" :> Quiet@Check[assoc["VertexArea"][hemi], $Failed],
            "Parcellation" :> Quiet@Check[assoc["Parcellation2009"][hemi], $Failed],
            "Parcellation2005" :> Quiet@Check[assoc["Parcellation"][hemi], $Failed],
            "RibbonIndices" :> Quiet@Check[Normal[VertexToVoxelMap[sub, hemi]][[All,2]], $Failed],
            "V1Label" :> Quiet@Check[assoc["V1Label"][hemi], $Failed],
            "V2Label" :> Quiet@Check[assoc["V2Label"][hemi], $Failed],
            "MTLabel" :> Quiet@Check[assoc["MTLabel"][hemi], $Failed],
            "BrodmannLabels" :> Quiet@Check[assoc["BrodmannLabels"][hemi], $Failed],
            "V1Probability" :> Quiet@Check[assoc["V1Probability"][hemi], $Failed],
            "V2Probability" :> Quiet@Check[assoc["V2Probability"][hemi], $Failed],
            "MTProbability" :> Quiet@Check[assoc["MTProbability"][hemi], $Failed],
            "BrodmannThresholds" :> Quiet@Check[assoc["BrodmannProbabilities"][hemi], $Failed]}]]]],
   (* We also want an accessor for MRImages *)
   MRImage[sub, name_, hemi:(None|LH|RH|RHX)] := With[
     {assoc = Association[sub],
      id = ToLowerCase[name] // Function @ FirstCase[
        Normal @ $FreeSurferImageData,
        (r_Rule /; Or[MatchQ[#, ToLowerCase[r[[1]]]],
                      MatchQ[#, r[[2]]["Pattern"]]] :> r[[1]])]},
     With[
       {vol0 = assoc[id]},
       With[
         {vol = Which[
            AssociationQ[vol0], If[hemi === None,
              If[KeyExistsQ[vol0, Full], 
                vol0[Full],
                ImageAdd[vol0[LH], vol0[RH]]],
              vol0[hemi]],
            ListQ[vol0], If[hemi === None,
              (Full /. vol0) /. Full :> ImageAdd[LH /. vol0, RH /. vol0],
              hemi /. vol0],
            MRImageQ[vol0], If[hemi === None,
              vol0,
              ImageMultiply[vol0, Image3D@assoc["FilledMaskImage"][hemi]]],
            True, Message[FreeSurferSubject::baddata, Path[sub], "Image not found for given FreeSurfer subject"]]},
         MRImage3D[
           vol,
           MetaInformation -> Join[
             Options[vol, MetaInformation],
             {"Subject" -> sub,
              "Hemisphere" -> hemi,
              "ImageName" -> id}]]]]],
   MRImage[sub, name_] := MRImage[sub, name, None],
               
   (* We can also get the occipital pole in a similar way... *)
   OccipitalPoleIndex[sub, hemi:(LH|RH|RHX)] := Check[
     FreeSurferSubjectOP[Path[sub], hemi],
     $Failed],
   (* We can get certain labels this way also *)
   LabelVertexList[sub, hemi:(LH|RH|RHX), name_] := Check[
     If[ListQ[name],
       (* custom label *)
       If[Length@Union[name, {0.0, 0, 1, 1.0, True, False}] == 6,
         Pick[
           Range[Length@name],
           Replace[name, {0.0|False -> 0, Except[0|0.0|False] -> 1}, {1}],
           1],
         name],
       (* builtin label *)
       With[
         {propAndPatt = Switch[
            name,
            "V1"|"V1Label", {"V1Label", 1},
            "V2"|"V2Label", {"V2Label", 1},
            "MT"|"MTLabel", {"MTLabel", 1},
            _, $Failed]},
         If[!ListQ[propAndPatt],
           propAndPatt,
           Flatten @ Position[
             Normal[Association[sub][propAndPatt[[1]]][hemi]],
             propAndPatt[[2]],
             {1},
             Heads -> False]]]],
     $Failed],
        
   SubjectLabels[sub] := Join[
     Intersection[
       FreeSurferSubjectBrodmannLabelList[Path[sub], LH],
       FreeSurferSubjectBrodmannLabelList[Path[sub], RH]],
     Select[
       Map[
         Function @ With[
           {files = Table[
              FileNameJoin[{Path[sub], "label", hem <> "." <> # <> ".thresh.label"}],
              {hem, {"lh","rh"}}]},
           If[FileExistsQ[files[[1]]] && FileExistsQ[files[[2]]], #, None]],
         {"V1","V2","MT"}],
       StringQ]],

   (* We want to be able to grab vertex/voxel mappings as well... *)
   FreeSurferSubjectVoxelVertexMapping[sub, hemi:(LH|RH)] := Check[
     With[
       {white = Cortex[sub, "White", hemi],
        pial = Cortex[sub, "Pial", hemi],
        ribbon = MRImage[sub, "Ribbon", hemi]},
       With[
         {idcs = Position[ImageData[ribbon], 1|1.0, {3,4}][[All, 1;;3]]},
         With[
           {xyz = VoxelIndexToCoordinate[ribbon, idcs]},
           Association @ MapThread[
             (#1 -> #2[[1]]) &,
             {idcs,
              Nearest[
                MapThread[
                  Rule,
                  {Join[VertexCoordinates[white], VertexCoordinates[pial]],
                   Join[#, #]& @ VertexList[white]}],
                xyz]}]]]],
     $Failed],
   FreeSurferSubjectVertexVoxelMapping[sub, hemi:(LH|RH)] := Check[
     With[
       {surf = Cortex[sub, "Middle", hemi],
        ribbon = MRImage[sub, "Ribbon", hemi]},
       With[
         {idcs = Position[ImageData[ribbon], 1|1.0, {3,4}][[All, 1;;3]]},
         Association @ MapThread[
           (#1 -> #2[[1]]) &,
           {VertexList[surf],
            Nearest[
              MapThread[
                Rule,
                {VoxelIndexToCoordinate[ribbon, idcs], idcs}],
              VertexCoordinates[surf]]}]]],
     $Failed],
   VertexToVoxelMaps[sub] :> Association[
     {LH -> FreeSurferSubjectVertexVoxelMapping[sub, LH],
      RH -> FreeSurferSubjectVertexVoxelMapping[sub, RH]}],
   VertexToVoxelMap[sub, hemi:(LH|RH)] := VertexToVoxelMaps[sub][hemi],
   VoxelToVertexMaps[sub] :> Association[
     {LH -> FreeSurferSubjectVoxelVertexMapping[sub, LH],
      RH -> FreeSurferSubjectVoxelVertexMapping[sub, RH]}],
   VoxelToVertexMap[sub, hemi:(LH|RH)] := VoxelToVertexMaps[sub][hemi]},
  SetSafe -> True,
  Symbol -> FreeSurferSubjectData];

(* We want to have a nice box-form for the subjects *)
MakeBoxes[s_FreeSurferSubjectData, form_] := MakeBoxes[#]&[
  With[
    {style = {
       FontSize -> 11,
       FontColor -> Gray,
       FontFamily -> "Arial",
       FontWeight -> "Thin"}},
    Row[
      {"FreeSurferSubject"[
         Panel @ Grid[
           {{Style[Path[s], Sequence@@style]}}, 
           Alignment -> {{Center}}]]},
      BaseStyle -> Darker[Gray]]]];

Protect[FreeSurferSubjectVoxelVertexMapping, FreeSurferSubjectVertexVoxelMapping];


(* FSAverage and FSAverageSym *********************************************************************)
$FSAverage := With[
  {possibles = Select[
     Keys[$FreeSurferSubjects],
     (Last[StringSplit[#, $PathnameSeparator]] == "fsaverage")&]},
  If[possibles == {},
    (Message[
       FreeSurfer::notfound,
       "no fsaverage subject found; you may beed to add freesurfer homes or subjects, or set" <> 
       " $FSAverage to the fsaverage subject directory manually."];
     $Failed),
    Set[$FSAverage, First[possibles]]]];
FSAverageSubject := With[
  {res = Check[
     FreeSurferSubject[$FSAverage], 
     $Failed]},
  If[res === $Failed, 
    res,
    (Unprotect[FSAverageSubject];
     FSAverageSubject = res;
     Protect[FSAverageSubject];
     res)]];

$FSAverageSym := With[
  {possibles = Select[
     Keys[$FreeSurferSubjects],
     (Last[StringSplit[#, $PathnameSeparator]] == "fsaverage_sym")&]},
  If[possibles == {},
    (Message[
       FreeSurfer::notfound,
       "no fsaverage_sym subject found; you may beed to add freesurfer homes or subjects, or" <> 
        " set $FSAverageSym to the fsaverage_sym subject directory manually."];
     $Failed),
    Set[$FSAverageSym, First[possibles]]]];
FSAverageSymSubject := With[
  {res = Check[FreeSurferSubject[$FSAverageSym], $Failed]},
  If[res === $Failed, 
    res,
    (Unprotect[FSAverageSymSubject];
     FSAverageSymSubject = res;
     Protect[FSAverageSymSubject];
     res)]];

Protect[FSAverageSubject, FSAverageSymSubject];

FSAverageOP[hemi:(LH|RH)] := With[
  {V = Check[VertexCoordinates[FSAverageSubject["InflatedSurface"][hemi]], $Failed]},
  If[V === $Failed,
    $Failed,
    Set[
      FSAverageOP[hemi],
      First[Ordering[V[[All, 2]], 1]]]]];
FSAverageSymOP := With[
  {V = Check[VertexVertexCoordinates[FSAverageSymSubject["InflatedSurface"][LH]], $Failed]},
  If[V === $Failed,
    $Failed,
    Set[
      FSAverageSymOP,
      First[Ordering[V[[All, 2]], 1]]]]];

(* #CortexToRibbon ********************************************************************************)
Options[CortexToRibbon] = {Filling -> 0};
CortexToRibbon[sub_, hemi:(LH|RH), dat_List, OptionsPattern[]] := With[
  {vtx2vox = VertexToVoxelMap[sub, hemi],
   ribbon = MRImage[sub, "Ribbon", hemi],
   pial = Cortex[sub, "Pial", hemi],
   fill = OptionValue[Filling]},
  MRImage3D[
    ribbon,
    ImageData -> ReplacePart[
      ImageData[ribbon],
      With[
        {rules = MapThread[
           (vtx2vox[#1] -> If[NumericQ[#2], #2, fill])&,
           {VertexList[pial], dat}]},
        If[MatchQ[Dimensions@ImageData[ribbon], {_,_,_,1}],
          Map[(Append[#[[1]],1] -> #[[2]])&, rules],
          rules]]]]];
Protect[CortexToRibbon];

(* #RibbonToCortex ********************************************************************************)
Options[RibbonToCortex] = {Method -> Mean};
RibbonToCortex[sub_, hemi:(LH|RH), img_] := With[
  {vox2vtx = VoxelToVertexMap[sub, hemi],
   pial = Cortex[sub, "Pial", hemi],
   aggf = OptionValue[Method],
   ribbon = MRImage[sub, "Ribbon", hemi]},
  With[
    {dat = Which[
       MRImageQ[img] && ImageDimensions[img] == ImageDimensions[ribbon], ImageData[img],
       ImageQ[img] && ImageDimensions[img] == ImageDimensions[ribbon], ImageData[img],
       ArrayQ[img, 3] && Dimensions[img] == ImageDimensions[ribbon], img,
       True, Message[
         RibbonToCortex::badarg,
         "image must be an MRImage, an Image3D, or a 3D array the same size as sub's ribbon"]]},
  Map[
    aggf[Extract[dat, vox2vtx[#]]]&,
    VertexList[pial]]]];
Protect[RibbonToCortex];

End[];
EndPackage[];

