(* Global.m
 *
 * The Neurotica`Global namespace includes global symbols and values that must be included by many
 * Neurotical sub-packages.
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
BeginPackage["Neurotica`Global`"];
Unprotect["Neurotica`Global`*", "Neurotica`Global`Private`*"];
ClearAll[ "Neurotica`Global`*", "Neurotica`Global`Private`*"];

Curvature::usage = "Curvature is a keyword that can be used to refer to curvature values; it is automatically defined by the CorticalSurface package to include a CorticalColor function as well.";

RH::usage = "RH is a keyword that represents the right hemisphere.";
LH::usage = "LH is a keyword that represents the left hemisphere.";

Cortex::usage = "Cortex is a keyword that is used by various Neurotica functions to represent the cortical surface of a subject.";
SubjectLabels::usage = "SubjectLabels[sub] yields a list of the labels supported by the given subject subject sub.";

OccipitalPole::usage = "OccipitalPole[subject, mesh, hemisphere] is usually defined by subject modalities (e.g., FreeSurferSubject[]) such that the function yields the vertex coordinate for the occipital pole in the particular mesh and hemisphere requested.
Note that if you define a new subject modality, then defining OccipitalPoleIndex[] for the subject should be sufficient.";
OccipitalPoleIndex::usage = "OccipitalPoleIndex[subject, hemisphere] is usually defined by subject modalities (e.g., FreeSurferSubject[]) such that the function yields the index for the occipital pole in the particular subject and hemisphere requested.";

LabelVertexList::usage = "LabelVertexList[sub, hemi, name] is defined by subject modalities (e.g., FreeSurferSubject[]) such that it yields a list of the vertices that are in the label with the given name for the given subject and hemisphere.
Note that the labels supported will vary by subject modality; use SubjectLabels to see the supported labels.";
LabelVertexList::nolab = "No such label for given subject: `1`";
LabelVertexCoordinates::usage = "LabelVertexCoordinates[sub, mesh, hemi, name] yields the vertex coordinates for the given subject, mesh, and hemisphere of the vertices that lie in the label with the given name.
Note that if you are defining a new subject modality, you should define LabelVertexList[] instead of LabelVertexCoordinates.";
LabelVertexCoordinatesTr::usage = "LabelVertexCoordinatesTr[sub, mesh, hemi, name] is equivalent to Transpose[LabelVertexCoordinates[sub, mesh, hemi, name]].
Note that if you are defining a new subject modality, you should define LabelVertexList[] instead of LabelVertexCoordinatesTr.";

LabelEdgePairsTr::usage = "LabelEdgePairsTr[sub, hemi, name] yields the equivalent of Transpose @ LabelEdgePairs[sub, hemi, name].";
LabelEdgePairs::usage = "LabelEdgePairs[sub, hemi, name] yields a list of the edge pairs that compose the label with the given name.";
LabelEdgeList::usage = "LabelEdgeList[sub, hemi, name] is equivalent to LabelEdgePairs[sub, hemi, name] except that it yields UndirectedEdge's instead of lists of vertex pairs.";

LabelFaceListTr::usage = "LabelFaceListTr[sub, hemi, name] yields the equivalent of Transpose @ LabelFaceList[sub, hemi, name].";
LabelFaceList::usage = "LabelFaceList[sub, hemi, name] yields the list of faces that are part of the label with the given name in the given subject and hemisphere.";

LabelBoundaryVertexList::usage = "LabelBoundaryVertexList[sub, hemi, name] yields a list of the vertices, in counter-clockwise order, that are part of the boundary between the given label and the outside of the label.";
LabelBoundaryVertexCoordinates::usage = "LabelBoundaryVertexCoordinates[sub, mesh, hemi, name] yields a list of the coordinates of the vertices, in counter-clockwise order, that appear in on the boundary of the given label in the given subject, mesh, and hemisphere.";
LabelBoundaryVertexCoordinatesTr::usage = "LabelBoundaryVertexCoordinatesTr[sub, mesh, hemi, name] is equivalent to Transpose @ LabelBoundaryVertexCoordinates[sub, mesh, hemi, name].";

LabelBoundaryEdgeList::usage = "LabelBoundaryEdgeList[sub, hemi, name] yields a list of the edges that form an approximate boundary to the region with the given label name for the given subject and hemisphere.
Note that if you are defining a new subject modality, then you only need to define the LabelVertexList function.";
LabelBoundaryEdgePairsTr::usage = "LabelBoundaryEdgePairsTr[sub, hemi, name] yields a list equivalent to Transpose[LabelBoundaryEdgePairs[sub, hemi]].";
LabelBoundaryEdgePairs::usage = "LabelBoundaryEdgePairs[sub, hemi, name] yields a list of the edge pairs rather than the edges themselves that are returned by LabelEdgeList.";

If[!ValueQ[Radius],
  (Radius::usage = "Radius is a keyword that is used to specify the distance from the center of a cortical map projection that should be included in the map projection.";
   Protect[Radius])];

Begin["`Private`"];

Protect[Curvature, RH, LH, OccipitalPoleIndex, LabelVertexList, SubjectLabels];

(* #Cortex ****************************************************************************************)
Cortex[sub_, hemi_] := Cortex[sub, Automatic, hemi];
Protect[Cortex];

(* #OccipitalPole *********************************************************************************)
OccipitalPole[sub_, mesh_, hemi_] := Check[
  VertexCoordinatesTr[Cortex[sub, mesh, hemi]][[All, OccipitalPoleIndex[sub, hemi]]],
  $Failed];
Protect[OccipitalPole];

(* #LabelVertexCoordinatesTr **********************************************************************)
LabelVertexCoordinatesTr[sub_, mesh_, hemi_, name_] := Check[
  With[
    {cortex = Cortex[sub, hemi, hemi]},
    Part[
      VertexCoordinatesTr[cortex],
      All,
      VertexIndex[cortex, LabelVertexList[sub, hemi, name]]]],
  $Failed];
Protect[LabelVertexCoordinatesTr];

(* #LabelVertexCoordinates ************************************************************************)
LabelVertexCoordinates[sub_, mesh_, hemi_, name_] := Check[
  Transpose[LabelVertexCoordinatesTr[sub, mesh, hemi, name]],
  $Failed];
Protect[LabelVertexCoordinates];

(* #LabelEdgePairsTr ******************************************************************************)
LabelEdgePairsTr[sub_, hemi_, name_] := Check[
  With[
    {U = LabelVertexList[sub, hemi, name],
     cortex = Cortex[sub, Automatic, hemi]},
    With[
      {UE = VertexEdgeList[cortex][[VertexIndex[cortex, U]]],
       Et = EdgePairsTr[cortex]},
      Et[[All, Select[Tally[Join@@UE], Last[#] == 2&][[All, 1]]]]]],
  $Failed];
Protect[LabelEdgePairsTr];

(* #LabelEdgePairs ********************************************************************************)
LabelEdgePairs[sub_, hemi_, name_] := Check[
  Transpose @ LabelEdgePairsTr[sub, hemi, name],
  $Failed];
Protect[LabelEdgePairs];

(* #LabelEdgeList *********************************************************************************)
LabelEdgeList[sub_, hemi_, name_] := Check[
  Apply[UndirectedEdge, #]& /@ Transpose[LabelEdgePairsTr[sub, hemi, name]],
  $Failed];
Protect[LabelEdgeList];

(* #LabelFaceListTr ******************************************************************************)
LabelFaceListTr[sub_, hemi_, name_] := Check[
  With[
    {U = LabelVertexList[sub, hemi, name],
     cortex = Cortex[sub, Automatic, hemi]},
    With[
      {UF = VertexFaceList[cortex][[VertexIndex[cortex, U]]],
       Ft = FaceListTr[cortex]},
      Ft[[All, Select[Tally[Join@@UF], Last[#] == 3&][[All, 1]]]]]],
  $Failed];
Protect[LabelFaceListTr];

(* #LabelFaceList *********************************************************************************)
LabelFaceList[sub_, hemi_, name_] := Check[
  Transpose @ LabelFaceListTr[sub, hemi, name],
  $Failed];
Protect[LabelFaceList];

(* #LabelBoundaryEdgePairsTr **********************************************************************)
LabelBoundaryEdgePairsTr[sub_, hemi_, name_] := Check[
  With[
    {Ft = LabelFaceListTr[sub, hemi, name]},
    With[
      {allEt = {
         Join@@Ft,
         Join[Ft[[2]], Ft[[3]], Ft[[1]]]}},
      With[
        {idcs = Select[
           Tally[
             Range[Length@allEt],
             (Length[Union[Fs[[All, #1]], Fs[[All, #2]]]] == 2)&],
           Last[#] == 1&]},
        allEt[[All, idcs]]]]],
  $Failed];
Protect[LabelBoundaryEdgePairsTr];

(* #LabelBoundaryEdgePairs ************************************************************************)
LabelBoundaryEdgePairs[sub_, hemi_, name_] := Check[
  Transpose @ LabelBoundaryEdgePairsTr[sub, hemi, name],
  $Failed];
Protect[LabelBoundaryEdgePairs];

(* #LabelBoundaryEdgeList *************************************************************************)
LabelBoundaryEdgeList[sub_, hemi_, name_] := Check[
  Apply[UndirectedEdge, #]& /@ Transpose[LabelBoundaryEdgePairsTr[sub, hemi, name]],
  $Failed];
Protect[LabelBoundaryEdgeList];

(* #LabelBoundaryVertexList ***********************************************************************)
LabelBoundaryVertexList[sub_, hemi_, name_] := Check[
  Append[#[[1]], #[[2, -1]]]& @ LabelBoundaryEdgePairsTr[sub, hemi, name],
Protect[LabelBoundaryVertexList];

(* #LabelBoundaryVertexCoordinatesTr **************************************************************)
LabelBoundaryVertexCoordinatesTr[sub_, mesh_, hemi_, name_] := Check[
  With[
    {cortex = Cortex[sub, mesh, hemi]},
    Part[
      VertexCoordinatesTr[cortex],
      All, 
      VertexIndex[cortex, LabelBoundaryVertexList[sub, hemi, name]]]],
  $Failed];
Protect[LabelBoundaryVertexCoordinatesTr];

(* #LabelBoundaryVertexCoordinates ****************************************************************)
LabelBoundaryVertexCoordinates[sub_, mesh_, hemi_, name_] := Check[
  Transpose @ LabelBoundaryVertexCoordinatesTr[sub, mesh_, hemi, name],
  $Failed];
Protect[LabelBoundaryVertexCoordinates];



End[];
EndPackage[];
