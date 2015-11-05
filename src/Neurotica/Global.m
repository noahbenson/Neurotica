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
LR::usage = "LR is a keyword that represents a pseudo-hemisphere, such as when both hemispheres are included in an MRImage (e.g., MRImage[sub, LR, \"Ribbon\"]).";
RHX::usage = "RHX is a keyword that represents the inverted right hemisphere as used by programs like FreeSurfer.";
LHX::usage = "LHX is a keyword that represents the inverted left hemisphere as used by programs like FreeSurfer.";
HemisphereQ::usage = "HemisphereQ[x] yields true if x is a valid hemisphere, i.e., RH, LH, RHX, or LHX.";

Anterior::usage = "Anterior is a keyword that represents the forward part of the brain; it is generally a synonym for Front and is an antonym with Posterior and Back.";
Posterior::usage = "Posterior is a keyword that represents the rear part of the brain; it is generally a synonym for Back and is an antonym with Anterior and Front.";
Superior::usage = "Superior is a keyword that represents the upper part of the brain; it is generally a synonym for Top and is an antonym with Inferior and Bottom.";
Inferior::usage = "Superior is a keyword that represents the lower part of the brain; it is generally a synonym for Bottom and is an antonym with Superior and Top.";

Eccentricity::usage = "Eccentricity is a key used by the visual cortex package to represent eccentricity, as measured in degrees of visual angle from the foveal confluence (center of the visual field).";
PolarAngle::usage = "PolarAngle is a key used by the visual cortex package to represent the polar angle, as measured in degrees of rotation about the foveal confluence (center of the visual field) from the upper to the lower vertical meridia.";
VisualArea::usage = "VisualArea is a key used by the retinotopy package to represent the visual area ID of a particular patch of cortex. See also VisualAreasData.";

Cortex::usage = "Cortex[sub, hemi, name] yields a cortical mesh object associated with the given hemisphere and name for the given subject. Cortex is a keyword that is used by various Neurotica functions to represent the cortical surface of a subject, so behavior is not guaranteed for non-standard subject modalities.";
SubjectLabels::usage = "SubjectLabels[sub] yields a list of the labels supported by the given subject subject sub.";

OccipitalPoleIndex::usage = "OccipitalPoleIndex[subject, hemisphere] is usually defined by subject modalities (e.g., FreeSurferSubject[]) such that the function yields the index for the occipital pole in the particular subject and hemisphere requested.";

LabelVertexList::usage = "LabelVertexList[sub, hemi, name] is defined by subject modalities (e.g., FreeSurferSubject[]) such that it yields a list of the vertices that are in the label with the given name for the given subject and hemisphere.
Note that the labels supported will vary by subject modality; use SubjectLabels to see the supported labels.";
LabelVertexList::nolab = "No such label for given subject: `1`";

VertexToVoxelMap::usage = "VertexToVoxelMap[sub, hemi] yields a mapping of vertices to voxels for the given subject. This is defined by each subject modality so may not be identical depending on how the subject is loaded.";
VoxelToVertexMap::usage = "VoxelToVertexMap[sub, hemi] yields a mapping of voxels to vertices for the given subject. This is defined by each subject modality so may not be identical depending on how the subject is loaded.";

If[!ValueQ[Radius],
  (Radius::usage = "Radius is a keyword that is used to specify the distance from the center of a cortical map projection that should be included in the map projection.";
   Protect[Radius])];

Begin["`Private`"];

HemisphereQ[x_] := False;
HemisphereQ[LH] = True;
HemisphereQ[RH] = True;
HemisphereQ[LR] = True;
HemisphereQ[LHX] = True;
HemisphereQ[RHX] = True;

Protect[RH, LH, LR, RHX, LHX, HemisphereQ, Anterior, Posterior, Inferior, Superior, PolarAngle,
        Eccentricity, VisualArea, Curvature, OccipitalPoleIndex, LabelVertexList, SubjectLabels,
        VertexToVoxelMap, VoxelToVertexMap];

(* #Cortex ****************************************************************************************)
Cortex[sub_, hemi_] := Cortex[sub, Automatic, hemi];
Protect[Cortex];

End[];
EndPackage[];
