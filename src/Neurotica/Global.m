(* Global.m
 *
 * The Neurotica`Global namespace includes global symbols and values that must be included by many
 * Neurotical sub-packages. To a large extent, this file fleshes out interfaces that are filled in
 * by other sub-packages.
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

Hemi::usage = "Hemi is a keyword used by Neurotica to indicate whether a piece of data relates to the left (LH) or right (RH) hemisphere. When the Hemi keyword is an option to a function, it indicates that one of these (LH or RH) is the appropriate value, as opposed to a hemisphere cortical mesh object.";
HemiQ::usage = "HemiQ[x] yields True if x is a valid hemisphere ID, i.e., RH, LH, RHX, or LHX. Note that a \"Hemi\" is a hemisphere tag while a \"Hemisphere\" is an object that encapsulates the data relevant to a subject's hemisphere. If you wish to test if an object is a hemisphere object, use HemisphereQ.";

RH::usage = "RH is a keyword that represents the right hemisphere.";
LH::usage = "LH is a keyword that represents the left hemisphere.";
LR::usage = "LR is a keyword that represents a pseudo-hemisphere, such as when both hemispheres are included in an MRImage (e.g., MRImage[sub, LR, \"Ribbon\"]).";
RHX::usage = "RHX is a keyword that represents the inverted right hemisphere as used by programs like FreeSurfer.";
LHX::usage = "LHX is a keyword that represents the inverted left hemisphere as used by programs like FreeSurfer.";

(* The Hemisphere interface: *)
Chirality::usage = "Chirality[hemi] yields the appropriate chirality of the given hemisphere object (which may be a hemisphere ID such as RH).";

(* The Subject interface: *)
SubjectQ::usage = "SubjectQ[s] yields True if s is a Neurotica subject object, otherwise False.";
Cortex::usage = "Cortex[sub, hemi] yields the cortical surface object associated with the given hemisphere.";
CortexList::usage = "CortexList[sub] yields a list of the cortical surface names assocaited with the given subject sub.";
HemiList::usage = "HemiList[sub] yields a list of hemi objects that are valid for the given subject sub.";

MRImage::usage = "MRImage[sub, hemi, name] yields the MRImage3D data for the given subject, hemisphere, and name (e.g., \"Ribbon\" or \"T1\"). The name argument varies by subject modality; they can be found via MRImageList[sub] or MRImageList[sub, hemi]. For FreeSurfer, available names include \"Brain\", \"Ribbon\", \"GrayMask\" and \"WhiteMatter\". Most subject modalities should define at least a \"GrayMask\".
As a general rule, all subject modalities should accept the names \"T1\", \"Brain\", and \"GrayMask\".
The option hemi may be LR or All to indicate the whole brain.
The options hemi, name, and orientation all have default parameters:
MRImage[sub, hemi] is equivalent to MRImage[sub, hemi, Automatic].
MRImage[sub, name] is equivalent to MRImage[sub, LR, Automatic].
MRImage[sub] is equivalent to MRImage[sub, LR, Automatic].

MRImage[data] yields an Image-like form that can be used with Image functions but which stores additional relevant data regarding MR images. Note this form of MRImage is intended for MRI slices rather than for 3D images; for a 3D image representation, use MRImage3D.";
MRImage::badarg = "Bad argument given to MRImage: `1`";
MRImageList::usage = "MRImageList[sub] yields a list of MRImage names for the given subject sub.";

LabelList::usage = "LabelList[sub] yields a list of the labels supported by the given subject sub.";
LabelVertexList::usage = "LabelVertexList[sub, hemi, name] is defined by subject modalities (e.g., FreeSurferSubject[]) such that it yields a list of the vertices that are in the label with the given name for the given hemisphere.

Note that the labels supported will vary by subject modality; use SubjectLabels to see the supported labels.";
LabelVertexList::nolab = "No such label found: `1`";
LabelVoxelList::usage = "LabelVoxelList[sub, hemi, name] yields a list of voxel indices for the label of the given subject, hemisphere, and name. If hemi is LR or All, then this is the union of LabelVoxelList[sub, LH, name] and LabelVoxelList[sub, RH, name].";

VertexPropertyAssociation::usage = "VertexPropertyAssociation[subject, hemi] yields the vertex property association for the given subject and hemisphere; the keys of this association are the property names for the subject.
VertexPropertyAssociation[mesh] yields an association whose keys are the property names of the vertices of the given mesh and whose values are the lists of properties for each vertex.";

OccipitalPoleIndex::usage = "OccipitalPoleIndex[hemisphere] is usually defined by subject modalities (e.g., FreeSurferSubject[]) such that the function yields the index for the occipital pole in the particular hemisphere requested.
OccipitalPoleIndex[sub, hemi] is equivalent to OccipitalPoleIndex@Hemisphere[sub, hemi].";

GrayMask::usage = "GrayVoxels[sub, hemi] yields a SparseArray representation of the gray-matter voxels in the given subject sub. In the SparseArray, gray voxels are given the value 1 and other voxels the value 0. The orientation may be any orientation accepted by MRImage. If a subject modality does not define this explicitly, it will default to MRImage[sub, hemi, \"GrayMask\"].";

CortexToRASMatrix::usage = "CortexToRASMatrix[sub] yields the transformation matrix that converts a vertex on the given subject's cortical surface into a RAS-oriented position.
CortexToRASMatrix[mesh] yields the equivalent matrix for a mesh, assuming that the element \"CortexToRASMatrix\" has been set in its meta-information.";

VertexToVoxelMap::usage = "VertexToVoxelMap[sub, hemi, name] yields a mapping of vertices of the cortical surface with the given name to ribbon voxels for the given subject and hemisphere. This is defined by each subject modality so may not be identical depending on how the subject is loaded.";
VoxelToVertexMap::usage = "VoxelToVertexMap[sub, hemi] yields a mapping of voxels to vertices for the given subject. This is defined by each subject modality so may not be identical depending on how the subject is loaded.";

Anterior::usage = "Anterior is a keyword that represents the forward part of the brain; it is generally a synonym for Front and is an antonym with Posterior and Back.";
Posterior::usage = "Posterior is a keyword that represents the rear part of the brain; it is generally a synonym for Back and is an antonym with Anterior and Front.";
Superior::usage = "Superior is a keyword that represents the upper part of the brain; it is generally a synonym for Top and is an antonym with Inferior and Bottom.";
Inferior::usage = "Superior is a keyword that represents the lower part of the brain; it is generally a synonym for Bottom and is an antonym with Superior and Top.";

Eccentricity::usage = "Eccentricity is a key used by the visual cortex package to represent eccentricity, as measured in degrees of visual angle from the foveal confluence (center of the visual field).";
PolarAngle::usage = "PolarAngle is a key used by the visual cortex package to represent the polar angle, as measured in degrees of rotation about the foveal confluence (center of the visual field) from the upper to the lower vertical meridia.";
VisualArea::usage = "VisualArea is a key used by the retinotopy package to represent the visual area ID of a particular patch of cortex. See also VisualAreasData.";

If[!ValueQ[Radius],
  (Radius::usage = "Radius is a keyword that is used to specify the distance from the center of a cortical map projection that should be included in the map projection.";
   Protect[Radius])];

Begin["`Private`"];

Protect[RH, LH, LR, RHX, LHX, Hemi, Anterior, Posterior, Inferior, Superior, PolarAngle,
        Eccentricity, VisualArea, Curvature, VertexToVoxelMap, VoxelToVertexMap];

(* #HemiQ *****************************************************************************************)
HemiQ[x_] := False;
HemiQ[LH] = True;
HemiQ[RH] = True;
HemiQ[LR] = True;
HemiQ[LHX] = True;
HemiQ[RHX] = True;
HemiQ[Left] = True;
HemiQ[Right] = True;
Protect[HemiQ];

(* #Chirality *************************************************************************************)
Chirality[LH] = LH;
Chirality[RH] = RH;
Chirality[LR] = None;
Chirality[LHX] = RH;
Chirality[RHX] = LH;
Chirality[Left] = LH;
Chirality[Right] = RH;
Protect[Chirality];

(* #SubjectQ **************************************************************************************)
SubjectQ[x_] := False;
Protect[SubjectQ];

(* #Cortex ****************************************************************************************)
Cortex[sub_, hemi_?HemiQ] := Cortex[sub, hemi, Automatic];
Cortex[hemi_?HemisphereQ] := Cortex[hemi, Automatic];
Protect[Cortex];

(* #GrayMask **************************************************************************************)
GrayMask[sub_, hemi_?HemiQ] := Check[MRImage[sub, hemi, "GrayMask"], $Failed];
Protect[GrayMask];

End[];
EndPackage[];
