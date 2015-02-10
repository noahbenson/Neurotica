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

GraphicsOptions::usage = "GraphicsOptions is a keyword that is used by immutable objects such as CorticalMesh to store default options for plotting the object.";

RH::usage = "RH is a keyword that represents the right hemisphere.";
LH::usage = "LH is a keyword that represents the left hemisphere.";

ProjectionArea::usage = "ProjectionArea is an option to CorticalProjection and related functions that specifies that the area of the projection should be limited to the given value; if this is not Full or All, then the projection algorithm excludes those faces, edges, and vertices that are farthest from the center of the projection, trimming them until the area restriction is met.";
ProjectionRadius::usage = "ProjectionRadius is an option to CorticalProjection and related fucntions that specifies that the radius of the projection should be limited to the given value; if this is not Full or All, then the projection algorithm excludes those faces, edges, and vertices that are farther (along the cortical surface) from the center of the projection than the given value.";

Cortex::usage = "Cortex is a keyword that is used by various Neurotica functions to represent the cortical surface of a subject.";

Begin["`Private`"];

Protect[Curvature, GraphicsOptions, RH, LH, ProjectionArea, ProjectionRadius, Cortex];

End[];
EndPackage[];
