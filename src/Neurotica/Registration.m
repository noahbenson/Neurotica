(* Registration.m
 *
 * The Neurotica`Registration` namespace contains functions related to the registration of surface
 * meshes to models defined on the cortical surface.
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
  "Neurotica`Registration`", 
  {"JLink`", "Neurotica`Global`","Neurotica`Util`","Neurotica`Mesh`"}];
Unprotect["Neurotica`Registration`*", "Neurotica`Registration`Private`*"];
ClearAll[ "Neurotica`Registration`*", "Neurotica`Registration`Private`*"];

CorticalPotentialFunction::usage = "CorticalPotentialFunction[{F, G, H}] yields a cortical potential function object with potential function F, gradient function G, and Hessian function H (which may be None). These functions are only ever called with a single argument, a 2 or 3 by n matrix, when the potential/gradient/hessian is evaluated. The following options may be given:
  * Print (default: Subscript[\"F\", \"Potential\"]) specifies what symbol should be used to display the potential function.
  * MetaInformation (default: {}) specifies any optional meta-information to attach to the potential function.
  * CorticalMesh (default: None) specifies the (optional) cortical mesh for which this potential function was defined.

Note: This function, along with the RegistrationTrajectory package, has been depricated and may be removed in the future. As an alternative, see the MeshRegister function and the Potential function.";
PotentialFunction::usage = "PotentialFunction[f] yields a pure functional form of the cortical potential function instance, f.

Note: This function, along with the RegistrationTrajectory package, has been depricated and may be removed in the future. As an alternative, see the MeshRegister function and the Potential function.";
GradientFunction::usage = "GradientFunction[f] yields a pure functional form of the gradient of the cortical potential function instance, f.

Note: This function, along with the RegistrationTrajectory package, has been depricated and may be removed in the future. As an alternative, see the MeshRegister function and the Potential function.";
HessianFunction::usage = "HessianFunction[f] yields a pure functional form of the Hessian of the cortical potential function instance, f.

Note: This function, along with the RegistrationTrajectory package, has been depricated and may be removed in the future. As an alternative, see the MeshRegister function and the Potential function.";
CorticalPotentialFunctionInstance::usage = "CorticalPotentialFunctionInstance is the head given to CorticalPotentialFunction objects.

Note: This function, along with the RegistrationTrajectory package, has been depricated and may be removed in the future. As an alternative, see the MeshRegister function and the Potential function.";

CalculateAngleIntermediateData::usage = "CalculateAngleIntermediateData[a,b,c] yields the numerical array equivalent to Flatten[{u1,u2,n1,n2,{d1,d2,cos}}, 1] where u1, u2, n1, and n2 are the vectors b - a, c - a, Normalize[b - a], and Normalize[c - a], respectively, and where d1 and d2 are the lengths of the b-a and c-a, respectively, and where cos is the cosine of the angle centered at a. Note that a, b, and c are expected to be 2 or 3 by n arrays where each column of the array corresponds to a single point.";
CalculateAngle::usage = "CalculateAngle[a,b,c] yields the angle between the vectors b-a and c-a; the three arguments must be 2D arrays with length 2 or 3 (for 2d or 3d points) and the return value is a vector of angles, one for each column of a,b, and c."
CalculateAngleGradient::usage = "CalculateAngleGradient[a,b,c] yields the gradient of the angle between the vectors b-a and c-a in terms of the vertices of a, b, and c; the three arguments must be 2D arrays with length 2 or 3 (for 2d or 3d points) and the return value is a 3x3xn or 3x2xn (for 3D or 2D points, respectively) matrix of gradients, one for each column of a,b, and c."
CalculateAngleHessian::usage = "CalculateAngleHessian[a,b,c] yields the Hessian of the angle between the vectors b-a and c-a in terms of the vertices of a, b, and c; the three arguments must be 2D arrays with length 2 or 3 (for 2d or 3d points) and the return value is a 3x3x3x3xn or 3x2x3x2xn (for 3D or 2D points, respectively) matrix of Hessians, one for each column of a,b, and c."

HarmonicEdgePotential::usage = "HarmonicEdgePotential[mesh] yields a function symbol f such that f[X] is the harmonic edge potential where X is a possible vertex coordinate list for the given cortical mesh. The potential is calculated as the total of (d - d0)^2 where d is the distance between two vertices in X and d0 is the distance in the coordinates of the given mesh. Note that Grad[f, X] yields the numerical gradient of the potential at the vertex configuration given in X.
The potential function of a HarmonicEdgePotential is U[d] = 0.5 / n * (d - d0)^2 where d is the distance between a pair of vertices connected by an edge, n is the number of edges in the system, and d0 is the distance, in the initial mesh, between the two vertices.

Note: This function, along with the RegistrationTrajectory package, has been depricated and may be removed in the future. As an alternative, see the MeshRegister function and the Potential function.";

HarmonicAnglePotential::usage = "HarmonicAnglePotential[mesh] yields a function symbol f such that f[X] is the harmonic angle potential where X is a possible vertex coordinate list for the given cortical mesh. The potential is calculated as the total of (a - a0)^2 where a is the angle of a face in X and a0 is the angle of the same face corner in the original mesh. Note that Grad[f, X] yields the numerical gradient of the potential at the vertex configuration given in X.
The potential function of a HarmonicAnglePotential is U[a] = 0.5 / n * (a - a0)^2 where a is the angle of a corner of a face, n is the number of faces in the system, and a0 is the angle of the same face in the initial mesh.

Note: This function, along with the RegistrationTrajectory package, has been depricated and may be removed in the future. As an alternative, see the MeshRegister function and the Potential function.";

VDWEdgePotential::usage = "VDWEdgePotential[mesh] is identical to HarmonicEdgePotential[mesh], but uses a much different potential function, based on the physical model of the van der Waals force. This potential function is experimental.

Note: This function, along with the RegistrationTrajectory package, has been depricated and may be removed in the future. As an alternative, see the MeshRegister function and the Potential function.";
VDWEdgePotential::badarg = "Bad argument given to VDWEdgePotential: `1`";

VDWAnglePotential::usage = "VDWAnglePotential[mesh] is identical to HarmonicAnglePotential[mesh], but uses a much different potential function, based on the physical model of the van der Waals force. This potential function is experimental.

Note: This function, along with the RegistrationTrajectory package, has been depricated and may be removed in the future. As an alternative, see the MeshRegister function and the Potential function.";


GaussianPotentialWell::usage = "GaussianPotentialWell[mesh, u -> {x0, std}] yields a function symbol f such that f[X] is the potential and Grad[f, X] is the gradient of an inverted Gaussian potential well that draws vertex u toward position x0 with the standard deviation std. In addition to the center and standard deviation, the following rules may be appended to the list on the right hand side of the rule:
  * \"FWHM\" (default: False) when True indicates that std should be interpreted as a full-width-half-max specification instead of the standard deviation.
  * \"Normalize\" (default: True) when True indicates that the Gaussian should be multiplied by 1 / (Sqrt[2 Pi] std).
  * \"Shape\" (default: 2) specifies that the Gaussian function should be of the form Exp[Abs[t / std]^q/q] where q is the shape of the generalized Gaussian.
  * \"Weight\" (default: 1) specifies that the Gaussian function should be weighted by the given number.
GaussianPotentialWell[mesh, {u1, u2, ...} -> {x0, std}] yields the Gaussian potential well function in which all of the vertices in the list on the left hand side of the rule are attracted to the specified Gaussian potential well on the right side and the potential and gradient yielded are divided by the number of vertices. Note that with this argument, the option \"Weight\" may be a list of numbers, one for each of the individual vertices.
GaussianPotentialWell[mesh, {vertices1 -> gaussian1, vertices2 -> gaussian2, ...}] yields the Gaussian potential well function that operates over all of the given vertex and Gaussian specifications.
Given the following values:
  * \[Beta] is the shape parameter,
  * \[Sigma] is the standard deviation (\[Sigma] = FWHM \[Beta]^(-1/\[Beta]) Log[2]^(-1/\[Beta])),
  * w is the weight, which has been multiplied by the normalizing term if the \"Normalize\" parameter is not set to False,
  * d[x] is equal to Norm[x - x0], the distance of the vertex from the center of the well,
the Gaussian potential function f and its gradient are defined as:
  * f[x] = -w Exp[-(d[x]/\[Sigma])^\[Beta] / \[Beta]]
  * \[Gradient]f[x] = w (d[x]/\[Sigma])^(\[Beta] - 1) / \[Sigma] Exp[-(d[x]/\[Sigma])^\[Beta] / \[Beta]]

Note: This function, along with the RegistrationTrajectory package, has been depricated and may be removed in the future. As an alternative, see the MeshRegister function and the Potential function.";
GaussianPotentialWell::badarg = "Bad argument given to GaussianPotentialWell: `1`";

HarmonicPotentialWell::usage = "HarmonicPotentialWell[mesh, u -> x0] yields a function symbol f such that f[X] is the potential and Grad[f, X] is the gradient of a harmonic potential well that draws vertex u toward position x0 with the potential function Norm[x - x0]^2.
HarmonicPotentialWell[mesh, u -> {x0}] is identical to HarmonicPotentialWell[mesh, u -> x0], but the following options may be appended to the list on the right hand side of the rule following x0:
  * \"Width\" (default 1) specifies the point at which the derivative of the potential well, in terms of the distance from the center of the well, equals 1.
  * \"Shape\" (default: 2) specifies that the harmonic .
  * \"Weight\" (default: 1) specifies that the Gaussian function should be weighted by the given number.
HarmonicPotentialWell[mesh, {u1, u2, ...} -> {x0, std}] yields the harmonic potential well function in which all of the vertices in the list on the left hand side of the rule are attracted to the specified harmonic potential well on the right side and the potential and gradient yielded are divided by the number of vertices. Note that with this argument, the option \"Weight\" may be a list of numbers, one for each of the individual vertices.
HarmonicPotentialWell[mesh, {vertices1 -> harmonic1, vertices2 -> harmonic2, ...}] yields the harmonic potential well function that operates over all of the given vertex and harmonic specifications.
Given the following values:
  * k is the width parameter,
  * \[Beta\] is the shape parameter,
  * w is the weight,
  * d[x] is equal to Norm[x - x0], the distance of the vertex to the center of the well,
the harmonic potential function f and its gradient are defined as:
  * f[x] = w/\[Beta] (d[x]/k)^\[Beta]
  * \[Gradient]f[x] = w/k (d[x]/k)^(\[Beta] - 1)

Note: This function, along with the RegistrationTrajectory package, has been depricated and may be removed in the future. As an alternative, see the MeshRegister function and the Potential function.";
HarmonicPotentialWell::badarg = "Bad argument given to HarmonicPotentialWell: `1`";

RegionDistancePotential::usage = "RegionDistancePotential[mesh, reg, {F, G}] yields a cortical potential function with potential function F and gradient G. Both F and G must be functions such that F[dists] and G[dists], for a vector dists of the distances of the relevant vertices to the region, yield a total potential and a vector of gradient magnitudes, respectively. Note that the direction of the gradient is calculated automatically. The following additional options may be given:
  * Print (default: Style[\"\[GothicCapitalG]\", FontWeight -> Bold]) specifies the default display name for the function.
  * MetaInformation (default: {}) specifies extra meta information attached to the function.
  * VertexWeight (default: Automatic) specifies the relative strength of each vertex in the field; the region will have a repulsive, neutral, or attractive effect on any vertex with a weight less than, equal to, or greater than 0, respectively. If a property is named by this argument, then its values are used. The default value, Automatic, applies the field to all vertices.
See also SignedRegionDistancePotential.

Note: This function, along with the RegistrationTrajectory package, has been depricated and may be removed in the future. As an alternative, see the MeshRegister function and the Potential function.";
RegionDistancePotential::badarg = "Bad argument given to RegionDistancePotential: `1`";

SignedRegionDistancePotential::usage = "SignedRegionDistancePotential[mesh, reg, {F, G}] is identical to RegionDistancePotential[mesh, reg, {F, G}] except that the functions F and G are given signed distances to the BoundaryMeshRegion reg for the relevant vertices instead of the absolute distances.

Note: This function, along with the RegistrationTrajectory package, has been depricated and may be removed in the future. As an alternative, see the MeshRegister function and the Potential function.";
SignedRegionDistancePotential::badarg = "Bad argument given to SignedRegionDistancePotential: `1`";

HarmonicPerimeterPotential::usage = "HarmonicPerimeterPotential[map] yields a cortical potential function that operates on the vertices on the perimeter of the given map to hold them in place using a harmonic potential well tied to their initial positions.

Note: This function, along with the RegistrationTrajectory package, has been depricated and may be removed in the future. As an alternative, see the MeshRegister function and the Potential function.";

PotentialField::usage = "PotentialField[mesh, direction, shape] yields a cortical potential function for use with the given mesh based on the requiested direction and shape of the potential field. Direction specifies the value over which the potential is computed while shape specifies the kind of shape that the potential takes.

Possible values of direction include:
  * \"Edges\": A function of the changes in the lengths of the edges is minimized.
  * \"Angles\": a function of the changes in the angles is minimized.
  * \"Anchors\" -> {vertices, points}: a function of the distance of the given set of vertices to the given set of fixed points is minimized.
  * \"Perimeter\": like \"Anchors\" but assumes that the vertices are the vertices along the perimeter of the (2D) mesh and that the points are the initial starting points of those vertices.
  * \"Mesh\": a combination of three potential fields: a harmonic edge potential (default scale: 500), an infinite-well angle potential (default scale: 1), and a harmonic perimeter potential with scale 1. This combination is generally useful for mesh deformation along with an anchor potential. If mesh is a 3D mesh, then the perimeter potential is not included.

Possible shapes include the following; note that x is used as the value of the variable and x0 is the reference value:
  * \"Harmonic\": the formula for a harmonic potential is (1/q Abs[x - x0]^q) where q is the order of the harmonic. The following options are available:
    * Order (default: 2) the order (q) of the harmonic.
  * \"Gaussian\": the formula for a Gaussian potential is (1 - Exp[-Abs[(x - x0)/s]^q]) where q is the order of the Gaussian and s is the standard deviation of the Gaussian. Options include:
    * Order (default: 2) the order (q) of the Gaussian.
    * StandardDeviation (default: 1) the standard deviation of the Gaussian.
  * \"LennardJones\": the formula for a Lennard-Jones potential is (1 + (x0/x)^q - 2 (x0/x)^(q/2)) where q is the order of the Lennard-Jones function. Options include:
    * Order (deatuls: 2) the order (q) of the Lennard-Jones function.
  * \"InfiniteWell\": the formula for an infinite-well potential is (((x0 - min)/(x - min))^(1/q))^2 + (((max - x0)/(min - x))^(1/q))^2 where q is the order of the infinite-well function, min is the minimum of the well, and max is the maximum of the well. Options include:
    * Min: (default: 0) the minimum of the well.
    * Max: (default: Pi) the maximum of the well; the default is Pi because infinite-well potentials are primarily used to prevent angles from increasing past Pi or below 0.
    * Order: (default: 0.5) the order of the infinite-well function.
  * \"Mesh\": two parameters are supported:
    * EdgePotentialScale (default: 500) specifies the scale of the edge potential field.
    * AnglePotentialScale (default: 1) specifies the scale of the angle potential field.

Note that if the potential direction is given as \"Mesh\", then no other options are required and, in fact, all other options are ignored with the exception of ReferenceCoordinates.

Additionally, all potential functions accept the following arguments:
  * Scale (default: 1) the scale of the potential field (this value is multiplied by the overall potential formula
  * ReferenceCoordinates (default: VertexCoordinates[mesh]) the reference coordinates to use when determining reference values.

Note that in all cases, a list of values may be given for the parameters, in which case, a separate parameterization is setup for each unit of the potential field. The unit for edges, angles, anchors, and perimeter is the edge, face, vertex, or perimeter-vertex, thus the length of the list must match the length of the corresponding unit list.";
PotentialField::badarg = "Bad argument given to Potential: `1`";
ReferenceCoordinates::usage = "ReferenceCoordinates is an option for the Potential function that specifies the reference coordinates for the mesh.";
PotentialFieldQ::usage = "PotentialFieldQ[p] yields True if p is a potential field and false otherwise.";
EdgePotentialScale::usage = "EdgePotentialScale is an option to PotentialField that specifies the scale of the edge potential.";
AnglePotentialScale::usage = "AnglePotentialScale is an option to PotentialField that specifies the scale of the angle potential.";

MeshRegister::usage = "MeshRegister[field] yields the mesh that results from minimizing the given potential field by warping the vertices of its associated mesh.
MeshRegister[mesh, {fieldDescriptions...}] registers the given mesh using the given set of field descriptions; each description should be a list consisting of arguments to PotentialField[] (excepting the mesh argument).
The following options may be given:
  * MaxSteps (default: 10,000) the maximum number of minimization steps to take
  * MaxStepSize (default: 0.1) the maximum distance any vertex should move in a single step
  * MaxPotentialChange (default: 1) the fraction of change in the potential that should be allowed to occur";
MeshRegister::badarg = "Bad argument given to MeshRegistar: `1`";
MaxPotentialChange::usage = "";

CortexGradientPlot::usage = "CortexGradientPlot[mesh, functions] yields a plot of the edges in the given mesh with the arrows representing the gradient of the vertices in the mesh, according to the list of potential functions given in the list functions. In addition to all options that are valid for CortexPlot, the following options may be given:
  * Arrowheads (default: Small) indicates that the arrowheads should be the given size in the plot.
  * PlotStyle (default: Automatic) should be a list of style instructions for the arrows of the gradients; these are cycled across the potential functions as in ListPlot.
  * Scaled (default: Automatic) indicates the absolute plotting length of the largest single gradient vector for any of the vertices; effectively, all gradients are scaled such that the largest gradient is equal to this value. If Automatic is given, then the value used is 75% of the mean edge length.";

CurvaturePotential::usage = "CurvaturePotential[mesh, template, sigma] yields a cortical potential function that enforces an aligment of the curvature in the cortical mesh, mesh, to the curvature in the cortical mesh, template, using a Gaussian smoothing function with shape parameter sigma over the cortical surface.

Note: This function, along with the RegistrationTrajectory package, has been depricated and may be removed in the future. As an alternative, see the MeshRegister function and the Potential function.";

MapTangledQ::usage = "MapTangledQ[map] yields True if and only if the given cortical map is tangled (has faces that are inverted); otherwise yields False.
MapTangledQ[map, X] is identical to MapTangledQ[map] except that it uses the coordinates given in X.";
MapTangles::usage = "MapTangles[map] yields a list of the vertex indices in the given cortical map that are tangles (their neighbors are not in the correct counter-clockwise ordering). Note that this yields the indices in map as opposed to the vertex id's.
MapTangles[map, X] uses the coordinates given in X for the map.";
MapUntangle::usage = "MapUntangle[map] attempts to untangle the map given by repeatedly moving tangled vertices to their centroids (with respect to their neighbors); this may not succeed, but will return the coordinates regardless after 50 such attempts.
MapUntangle[map, X] uses the coordinates X as the map coordinates.
MapUntangle[map, X, max] attempts at most max times to untangle the map.";

RegistrationFrame::usage = "RegistrationFrame[PF,X0] yields a registration frame with the given potential function PF and starting coordinates X0.

Note: The RegistrationFrame, RegistrationTrajectory, and related functions are now depricated in favor of MeshRegister. They may be removed in the future.";
RegistrationFrame::badarg = "Bad argument given to RegistrationFrame: `1`";
RegistrationFrameData::usage = "RegistrationFrameData[...] is used to represent a registration frame object.";
RegistrationFrameQ::usage = "RegistrationFrameQ[frame] yields True if and only if frame is a valid registration frame; otherwise yields False.

Note: The RegistrationFrame, RegistrationTrajectory, and related functions are now depricated in favor of MeshRegister. They may be removed in the future.";
StepNumber::usage = "StepNumber[frame] yields the step number of the given registration frame.";
StepSize::usage = "StepSize[frame] yields the recommended step-size for the given frame.";
MinStepSize::usage = "MinStepSize[frame] yields the minimum step-size the registration will descend to; once the step-size is below this length, the registration is considered converged.";
MaxVertexChange::usage = "MaxVertexChange[frame] yields the maximum distance a vertex is allowed to travel in a single step of the registration.";

RegistrationTrajectory::usage = "RegistrationTrajectory[mesh, F] yields a registration trajectory object for the given cortical mesh and given potential field F. The following options may be given:
  * InitialVertexCoordinates (default: Automatic) specifies the starting coordinates in the registration; the default, Automatic, specifies that it should start with VertexCoordinates[mesh].
  * MaxVertexChange (default: 0.1) specifies the maximum distance a single vertex is allowed to move during a single step of the minimization.
  * MinStepSize (default: 10^-5) specifies that the minimum gradient length that must be reached for the minimization to be considered complete.
  * CacheFrequency (default: 50) specifies that the registration will save every <n>'th frame during minimization.

Note: The RegistrationFrame, RegistrationTrajectory, and related functions are now depricated in favor of MeshRegister. They may be removed in the future.";
RegistrationTrajectory::badarg = "Bad argument given to RegistrationTrajectory: `1`";
RegistrationTrajectory::nocnv = "Failure to converge: `1`";
RegistrationTrajectoryData::usage = "RegistrationTrajectoryData[...] is used to represent a registration trajectory object.";
RegistrationTrajectoryQ::usage = "RegistrationTrajectoryQ[traj] yields True if and only if traj is a RegistrationTrajectory object.";
CacheFrequency::usage = "CacheFrequency is an argument given to RegistrationTrajectory that specifies how often a frame in the registration should be cached; by default the value is 100.";
InitialFrame::usage = "InitialFrame[traj] yields the initial frame of the given registration trajectory, traj.";
FinalFrame::usage = "FinalFrame[traj] yields the final frame of the given registration trajectory, traj; note that this may take a long time to calculate.";
InitialVertexCoordinates::usage = "InitialVertexCoordinates is an option to RegistrationTrajectory and MapRegister that indicates the coordinates that should be used in the initial frame.";

Begin["`Private`"];

(* #Java ******************************************************************************************)
(* Here, we load the Java libraries that we require... start by putting nben.jar on the classpath *)
With[
  {nben = Select[JavaClassPath[], Last@FileNameSplit[#] == "nben.jar" &]},
  If[!ListQ[nben] || nben == {},
    (* we need to add it... we can do this via the $InputFileName variable *)
    AddToClassPath@FileNameJoin@Join[
      Most@FileNameSplit[$InputFileName], 
      {"lib", "nben", "target"}]]];
(* Okay, now we need to make sure to load some classes; if these raise exceptions, then probably
   git submodules were not initialized *)
If[!Check[
     Quiet[
       LoadJavaClass["nben.mesh.registration.Minimizer"];
       LoadJavaClass["nben.mesh.registration.Fields"];
       LoadJavaClass["nben.mesh.registration.Util"];
       True],
     False],
  Message[
    Neurotica::initwarn,
    StringJoin[
      "Neurotica's Java registration subsystem failed to initialize, which breaks the ",
      "MeshRegister and PotentialField functions. This is usually because git submodules ",
      "were not initialized and updated. See the Neurotica README for more information."]]];

(* #CorticalPotentialFunction *********************************************************************)
Options[CorticalPotentialFunction] = {
  Print -> Subscript["F", "Potential"],
  MetaInformation -> {},
  CorticalMesh -> None};
DefineImmutable[
  CorticalPotentialFunction[{F_, G_, H_}, OptionsPattern[]] :> P,
  {(* Let's start with options! *)
   Options[P] = Map[# -> OptionValue[#]&, Options[CorticalPotentialFunction][[All,1]]],
   Options[P, opt_] := Replace[
     opt,
     Append[Options[P], x_ :> (Message[Options::optionf, x, CorticalPotentialFunction]; {})]],
   (* We need to define the potential functions and gradients... *)
   PotentialFunction[P] -> F,
   GradientFunction[P]  -> G,
   HessianFunction[P]   -> H,
   (* And a call form for the gradient. *)
   Grad[P, M_ /; ArrayQ[M, 2, NumericQ]] := With[
     {f = GradientFunction[P]},
     Join @@ If[Length[M] <= Length@First[M], f[M], Transpose@f@Transpose[M]]],
   (* And a call form for the Hessian. *)
   Hessian[P, M_ /; ArrayQ[M, 2, NumericQ]] := With[
     {f = HessianFunction[P]},
     If[Length[M] <= Length[M[[1]]], f[M], Transpose[f[Transpose @ M]]]]},
  SetSafe -> True,
  Symbol -> CorticalPotentialFunctionInstance];
Protect[PotentialFunction, GradientFunction];

(* We have a few more edits for the potential function though: *)
Unprotect[CorticalPotentialFunctionInstance];

(* The call form for the potential is here, where it can be defined: *)
(CPF:CorticalPotentialFunctionInstance[__])[X_ /; ArrayQ[X, 2, NumericQ]] := With[
   {f = PotentialFunction[CPF]},
   If[Length[X] <= Length[X[[1]]], f[X], f[Transpose @ X]]];
(CPF:CorticalPotentialFunctionInstance[__])[X_Symbol] := With[
   {tmp = TemporarySymbol["CPF"],
    f = PotentialFunction[CPF]},
   tmp[M_ /; ArrayQ[M, 2, NumericQ]] := f[M];
   tmp[X]];

(* Critically, we want a nice display form for these potential functions: *)
MakeBoxes[P_CorticalPotentialFunctionInstance, form_] := MakeBoxes[#]& @ Options[P, Print];

(* Make sure that these functions combine properly: *)
CorticalPotentialFunctionInstance /: Plus[P0_CorticalPotentialFunctionInstance,
                                          P1_CorticalPotentialFunctionInstance,
                                          rest___] := Plus[
  CorticalPotentialFunction[
    {Function[PotentialFunction[P0][#] + PotentialFunction[P1][#]],
     Function[GradientFunction[P0][#] + GradientFunction[P1][#]],
     If[HessianFunction[P0] === None || HessianFunction[P1] === None,
       None,
       Function[HessianFunction[P0][#] + HessianFunction[P1][#]]]},
    Print -> Row[{Print /. Options[P0], "+", Print /. Options[P1]}],
    CorticalMesh -> With[
      {msh = CorticalMesh /. Options[P0]},
      If[msh == (CorticalMesh /. Options[P1]), msh, None]],
    MetaInformation -> ((MetaInformation /. #)& /@ {Options[P0], Options[P1]})],
  rest];
CorticalPotentialFunctionInstance /: Plus[x_?NumericQ,
                                          P_CorticalPotentialFunctionInstance,
                                          r___] := Plus[
  CorticalPotentialFunction[
    {Function[x + PotentialFunction[P][#]],
     GradientFunction[P],
     HessianFunction[P]},
    Sequence@@Options[P]],
  r];
CorticalPotentialFunctionInstance /: Times[x_?NumericQ,
                                           P_CorticalPotentialFunctionInstance,
                                           r___] := Times[
  CorticalPotentialFunction[
    {Function[x * PotentialFunction[P][#]],
     Function[x * GradientFunction[P][#]],
     If[HessianFunction[P] === None, None, Function[x * HessianFunction[P][#]]]},
    Sequence@@Options[P]],
  r];
Protect[CorticalPotentialFunction, CorticalPotentialFunctionInstance];


(* ============================================================================================== *)
(* ====================================== Private Functions ===================================== *)
(* ============================================================================================== *)

(* Compute distance intermediate data *)
CalculateDistance = Compile[
  {{a, _Real, 2}, {b, _Real, 2}},
  Sqrt @ Total[(b - a)^2],
  RuntimeOptions -> {"Speed", "EvaluateSymbolically" -> False}, 
  Parallelization -> True];
CalculateDistanceGradient = Compile[
  {{a, _Real, 2}, {b, _Real, 2}},
  With[
    {ab = b - a},
    With[
      {r = Sqrt @ Total[ab^2]},
      With[
        {da = -ab / ConstantArray[r, Length[a]]},
        {da, -da}]]],
  RuntimeOptions -> {"Speed", "EvaluateSymbolically" -> False}, 
  Parallelization -> True];
CalculateDistanceHessian = Compile[
  {{a, _Real, 2}, {b, _Real, 2}},
  With[
    {ba = a - b,
     dims = Length[a],
     n = Length@First[a]},
    With[
      {r = Sqrt @ Total[ba^2]},
      With[
        {denom = ConstantArray[1/r^3, {2*dims, 2*dims}],
         rdiag = With[
           {r2 = r^2, zeros = ConstantArray[0, n]},
           Table[If[i == j, r2, zeros], {i, 1, dims}, {j, 1, dims}]]},
        With[
          {raa = rdiag - Outer[Times, ba, ba, 1, 1],
           rab = Outer[Times, ba, ba, 1, 1] - rdiag},
          denom * Join[
            MapThread[Join, {raa, rab}],
            MapThread[Join, {Transpose@rab, raa}]]]]]],
  RuntimeOptions -> {"Speed", "EvaluateSymbolically" -> False},
  Parallelization -> True];
SetAttributes[CalculateDistance, NumericFunction];
SetAttributes[CalculateDistanceGradient, NumericFunction];
SetAttributes[CalculateDistanceHessian, NumericFunction];
Protect[CalculateDistance, CalculateDistanceGradient, CalculateDistanceHessian];

(* Compute the gradient of the angles in terms of the vertices *)
CalculateAngleIntermediateData = Compile[
  {{a, _Real, 2}, {b, _Real, 2}, {c, _Real, 2}},
  With[
    {u1 = b - a, u2 = c - a, u3 = c - b,
     dim = Length[a]},
    With[
      {d1 = (# + (1 - Unitize[#])) &@Chop@Re@Sqrt[Total[u1^2]],
       d2 = (# + (1 - Unitize[#])) &@Chop@Re@Sqrt[Total[u2^2]],
       d3 = (# + (1 - Unitize[#])) &@Chop@Re@Sqrt[Total[u3^2]]},
      With[
        {cos = Total[u1*u2]/(d1*d2),
         n1 = u1/ConstantArray[d1, dim],
         n2 = u2/ConstantArray[d2, dim],
         n3 = u3/ConstantArray[d3, dim]},
        Join[u1, u2, u3, n1, n2, n3, {d1, d2, d3, cos, Re@ArcCos[cos]}]]]],
  RuntimeOptions -> {"Speed", "EvaluateSymbolically" -> False}, 
  Parallelization -> True];
SetAttributes[CalculateAngleIntermediateData, NumericFunction];
Protect[CalculateAngleIntermediateData];
CalculateAngle = Compile[
  {{a, _Real, 2}, {b, _Real, 2}, {c, _Real, 2}},
  With[
    {data = CalculateAngleIntermediateData[a,b,c]},
    Last[data]],
  {{CalculateAngleIntermediateData[_,_,_], _Real, 2}},
  RuntimeOptions -> {"Speed", "EvaluateSymbolically" -> False}, 
  Parallelization -> True];
SetAttributes[CalculateAngle, NumericFunction];
Protect[CalculateAngle];
CalculateAngleGradientWithData = Compile[
  {{a, _Real, 2}, {b, _Real, 2}, {c, _Real, 2}, {data, _Real, 2}},
  With[
    {dims = Length[a]},
    With[
      {uAB = data[[If[dims == 2, {1,2},   {1,2,3}]]],
       uAC = data[[If[dims == 2, {3,4},   {4,5,6}]]],
       uBC = data[[If[dims == 2, {5,6},   {7,8,9}]]],
       nAB = data[[If[dims == 2, {7,8},   {10,11,12}]]],
       nAC = data[[If[dims == 2, {9,10},  {13,14,15}]]],
       nBC = data[[If[dims == 2, {11,12}, {16,17,18}]]],
       dAB = data[[6*dims + 1]],
       dAC = data[[6*dims + 2]],
       dBC = data[[6*dims + 3]],
       cos = data[[6*dims + 4]],
       th  = data[[6*dims + 5]]},
      With[
        {sin = Re@Sqrt[1 - cos^2],
         unit = ConstantArray[Unitize[Chop[1 - Abs[cos]]], dims]},
        With[
          {sinAB = Chop[dAB*sin],
           sinAC = Chop[dAC*sin]},
          With[
            {f1 = Times[
               unit,
               ConstantArray[cos, dims]*nAB - nAC,
               ConstantArray[1.0/(sinAB + (1 - Unitize[sinAB])), dims]],
             f2 = Times[
               unit,
               ConstantArray[cos, dims]*nAC - nAB,
               ConstantArray[1.0/(sinAC + (1 - Unitize[sinAC])), dims]]},
            {-(f1 + f2), f1, f2}]]]]],
  RuntimeOptions -> {"Speed", "EvaluateSymbolically" -> False}, 
  Parallelization -> True];
CalculateAngleGradient = Compile[
  {{a, _Real, 2}, {b, _Real, 2}, {c, _Real, 2}},
  CalculateAngleGradientWithData[a,b,c, CalculateAngleIntermediateData[a,b,c]],
  {{CalculateAngleIntermediateData[_,_,_], _Real, 2},
   {CalculateAngleGradientWithData[_,_,_,_], _Real, 3}},
  RuntimeOptions -> {"Speed", "EvaluateSymbolically" -> False}, 
  Parallelization -> True];
SetAttributes[CalculateAngleGradient, NumericFunction];
SetAttributes[CalculateAngleGradientWithData, NumericFunction];
Protect[CalculateAngleGradient, CalculateAngleGradientWithData];
(* Used in calculating the Hessian... *)
NormalizedVectorGradient = Compile[
  {{u, _Real, 2}, {v, _Real, 2}},
  (* Gives D[Normalize[v-u], {u}] == -D[Normalize[v-u], {v}] *)
  With[
    {diff = u - v},
    With[
      {norm = Sqrt@Total[diff^2]},
      Divide[
        Table[
          If[i == j,
            diff[[i]] * diff[[j]] - norm^2,
            diff[[i]] * diff[[j]]],
          {i,1,Length@u},
          {j,1,Length@v}],
        ConstantArray[norm^3, {Length@u, Length@v}]]]],
  RuntimeOptions -> {"Speed", "EvaluateSymbolically" -> False}, 
  Parallelization -> True];
SetAttributes[NormalizedVectorGradient, NumericFunction];
Protect[NormalizedVectorGradient];
VectorNormalizationFactorGradient = Compile[
  {{u, _Real, 2}, {v, _Real, 2}},
  (* Gives D[1 / Norm[v-u], {u}] == -D[1 / Norm[v-u], {v}]*)
  With[
    {diff = v - u},
    With[
      {norm = Sqrt@Total[diff^2]},
      diff / ConstantArray[norm^3, Length@diff]]],
  RuntimeOptions -> {"Speed", "EvaluateSymbolically" -> False}, 
  Parallelization -> True];
SetAttributes[VectorNormalizationFactorGradient, NumericFunction];
Protect[VectorNormalizationFactorGradient];
CalculateAngleHessianWithData = Function[{a,b,c,data},(*Compile[
  {{a, _Real, 2}, {b, _Real, 2}, {c, _Real, 2}, {data, _Real, 2}},*)
  With[
    {grad = CalculateAngleGradientWithData[a,b,c,data],
     dims = Length[a]},
    With[
      {uAB = data[[If[dims == 2, {1,2},   {1,2,3}]]],
       uAC = data[[If[dims == 2, {3,4},   {4,5,6}]]],
       uBC = data[[If[dims == 2, {5,6},   {7,8,9}]]],
       nAB = data[[If[dims == 2, {7,8},   {10,11,12}]]],
       nAC = data[[If[dims == 2, {9,10},  {13,14,15}]]],
       nBC = data[[If[dims == 2, {11,12}, {16,17,18}]]],
       dAB = data[[6*dims + 1]],
       dAC = data[[6*dims + 2]],
       dBC = data[[6*dims + 3]],
       cos = data[[6*dims + 4]],
       th  = data[[6*dims + 5]]},
      With[
        {csc = ConstantArray[Csc[th], dims],
         cot = ConstantArray[Cot[th], dims],
         grAnAB = NormalizedVectorGradient[a, b],
         grAnAC = NormalizedVectorGradient[a, c],
         grBhBC = VectorNormalizationFactorGradient[b, c]},
        With[
          {commonBC = nAC * csc - nAB * cot,
           commonCB = nAB * csc - nAC * cot,
           dBCfull = ConstantArray[dBC, dims]},
          With[
            {thBA = ConstantArray[1/dBCfull, dims] * Plus[
               -ConstantArray[csc, dims] * grAnAC,
                ConstantArray[csc*cot, dims] * Outer[Times, grad[[1]], nAC, 1, 1],
                ConstantArray[cot, dims] * grAnAB,
               -ConstantArray[csc*csc, dims] * Outer[Times, grad[[1]], nAB, 1, 1]],
             thBB = Plus[
                ConstantArray[cot, dims] * Outer[Times, grBhBC, nAB, 1, 1],
                ConstantArray[cot / dBCfull, dims] * (-grAnAB),
               -ConstantArray[csc*csc / dBCfull, dims] * Outer[Times, grad[[2]], nAB, 1, 1],
               -ConstantArray[csc, dims] * Outer[Times, grBhBC, nAC, 1, 1],
                ConstantArray[csc*cot / dBCfull, dims] * Outer[Times, grad[[2]], nAC, 1, 1]],
             thBC = Plus[
               Outer[Times, -nBC / ConstantArray[dBC^2, dims], commonBC, 1, 1],
               ConstantArray[dBC * csc[[1]], {dims, dims}] * Plus[
                 -grAnAC,
                 -Outer[Times, grad[[3]], cot * nAC, 1, 1],
                 Outer[Times, grad[[3]], nAB * csc, 1, 1]]],
             thCA = ConstantArray[1.0 / dBC, {dims, dims}] * Plus[
               grAnAB * ConstantArray[csc, dims],
               -ConstantArray[cot*csc, dims] * Outer[Times, nAB, grad[[1]], 1, 1],
               -grAnAC*ConstantArray[cot, dims],
               ConstantArray[csc^2, dims] * Outer[Times, nAC, grad[[1]], 1, 1]],
             thCC = Plus[
               Outer[Times, nBC / ConstantArray[dBC^2, dims], commonCB, 1, 1],
               ConstantArray[1.0/dBC, {dims, dims}] * Plus[
                 -Outer[Times, grad[[3]] * cot * csc, nAB, 1, 1],
                 grAnAC * ConstantArray[cot, dims],
                 Outer[Times, grad[[3]], nAC * csc^2, 1, 1]]]},
            {{-(thBA + thCA), Transpose@thBA, Transpose@thCA},
             {thBA, thBB, thBC},
             {thCA, Transpose@thBC, thCC}}]]]]](*,
  {{CalculateAngleGradientWithData[_,_,_,_], _Real, 3},
   {NormalizedVectorGradient[_,_], _Real, 3},
   {VectorNormalizationFactorGradient[_,_], _Real, 2}},
  RuntimeOptions -> {"Speed", "EvaluateSymbolically" -> False}, 
  Parallelization -> True*)];
CalculateAngleHessian = Compile[
  {{a, _Real, 2}, {b, _Real, 2}, {c, _Real, 2}},
  CalculateAngleHessianWithData[a,b,c, CalculateAngleIntermediateData[a,b,c]],
  {{CalculateAngleIntermediateData[_,_,_], _Real, 2},
   {CalculateAngleHessianWithData[_,_,_,_], _Real, 5}},
  RuntimeOptions -> {"Speed", "EvaluateSymbolically" -> False}, 
  Parallelization -> True];
SetAttributes[CalculateAngleHessian, NumericFunction];
SetAttributes[CalculateAngleHessianWithData, NumericFunction];
Protect[CalculateAngleHessian, CalculateAngleHessianWithData];

(* Here we compile functions for calculating harmonic angle potentials *)
CalculateHarmonicAnglePotential = Compile[
  {{a, _Real, 2}, {b, _Real, 2}, {c, _Real, 2}, {th0, _Real, 1}},
  With[
    {data = CalculateAngleIntermediateData[a, b, c]},
    Total[0.5 * (th0 - Last[data])^2]],
  {{CalculateAngleIntermediateData[_,_,_], _Real, 2}},
  RuntimeOptions -> {"Speed", "EvaluateSymbolically" -> False}, 
  Parallelization -> True];
CalculateHarmonicAngleGradient = Compile[
  {{a, _Real, 2}, {b, _Real, 2}, {c, _Real, 2}, {th0, _Real, 1}},
  With[
    {data = CalculateAngleIntermediateData[a, b, c],
     dims = Length[a]},
    With[
      {th = Last[data],
       dth = CalculateAngleGradientWithData[a, b, c, data]},
      dth * ConstantArray[th - th0, {3, dims}]]],
  {{CalculateAngleIntermediateData[_,_,_], _Real, 2},
   {CalculateAngleGradientWithData[_,_,_,_], _Real, 3}},
  RuntimeOptions -> {"Speed", "EvaluateSymbolically" -> False}, 
  Parallelization -> True];
CalculateHarmonicAngleHessian = Compile[
  {{a, _Real, 2}, {b, _Real, 2}, {c, _Real, 2}, {th0, _Real, 1}},
  -CalculateAngleHessian[a,b,c],
  {{CalculateAngleHessian[_,_,_], _Real, 4}},
  RuntimeOptions -> {"Speed", "EvaluateSymbolically" -> False}, 
  Parallelization -> True];
SetAttributes[CalculateHarmonicAnglePotential, NumericFunction];
SetAttributes[CalculateHarmonicAngleGradient, NumericFunction];
SetAttributes[CalculateHarmonicAngleHessian, NumericFunction];
Protect[CalculateHarmonicAnglePotential, CalculateHarmonicAngleGradient,
        CalculateHarmonicAngleHessian];

CalculateVDWAnglePotential = Compile[
  {{a, _Real, 2}, {b, _Real, 2}, {c, _Real, 2}, {th0, _Real, 1}},
  With[
    {data = CalculateAngleIntermediateData[a,b,c]},
    With[
      {ratio = (# + (1 - Unitize[#]))&[2.0 * Last[data] / th0]},
      Total[0.25 + 1.0/ratio^2 - 1.0/ratio]]],
  {{CalculateAngleIntermediateData[_,_,_], _Real, 2}},
  RuntimeOptions -> {"Speed", "EvaluateSymbolically" -> False},
  Parallelization -> True];
CalculateVDWAngleGradient = Compile[
  {{a, _Real, 2}, {b, _Real, 2}, {c, _Real, 2}, {th0, _Real, 1}},
  With[
    {data = CalculateAngleIntermediateData[a,b,c],
     dims = Length[a]},
    With[
      {th = Last[data],
       dth = CalculateAngleGradientWithData[a,b,c,data]},
      dth * ConstantArray[0.5 * (th0 / th^2 - th0^2 / th^3), {3, dims}]]],
  {{CalculateAngleIntermediateData[_,_,_], _Real, 2},
   {CalculateAngleGradientWithData[_,_,_,_], _Real, 3}},
  RuntimeOptions -> {"Speed", "EvaluateSymbolically" -> False},
  Parallelization -> True];
CalculateVDWAngleHessian = Compile[
  {{a, _Real, 2}, {b, _Real, 2}, {c, _Real, 2}, {th0, _Real, 1}},
  With[
    {data = CalculateAngleIntermediateData[a,b,c],
     dims = Length[a]},
    With[
      {th = Last[data]
       Hth = CalculateAngleHessianWithData[a,b,c,data]},
      Hth * ConstantArray[0.5 * th0 * (3*th0 - 2*th) / th^4, {3, dims, 3, dims}]]],
  {{CalculateAngleIntermediateData[_,_,_], _Real, 2},
   {CalculateAngleGradientWithData[_,_,_,_], _Real, 3}},
  RuntimeOptions -> {"Speed", "EvaluateSymbolically" -> False},
  Parallelization -> True];
SetAttributes[CalculateVDWAnglePotential, NumericFunction];
SetAttributes[CalculateVDWAngleGradient, NumericFunction];
SetAttributes[CalculateVDWAngleHessian, NumericFunction];
Protect[CalculateVDWAnglePotential, CalculateVDWAngleGradient, CalculateVDWAngleHessian];


(* These functions are helpers for the Gaussian potential well functions *)
Options[GaussianPotentialWellParseGaussian] = {
  "FWHM" -> False,
  "Normalize" -> True,
  "Shape" -> 2,
  "Weight" -> 1};
GaussianPotentialWellParseGaussian[x0_, sigma_, OptionsPattern[]] := With[
  {fwhm = OptionValue["FWHM"],
   normalize = OptionValue["Normalize"],
   shape = OptionValue["Shape"],
   weight = OptionValue["Weight"],
   dims = Length[x0]},
  With[
    {sig = If[TrueQ[fwhm], 2.0*Sqrt[2.0*Log[2.0]] * sigma, sigma]},
    Block[
      {X, U, distance},
      {ReplacePart[ (* Potential function *)
         Hold[
           {{X, _Real, 2}},
           Evaluate @ Hold[
             {distance = Sqrt @ Total[(X - x0)^2]},
             Evaluate @ Times[
               -weight * Exp[-(distance / sig)^shape / shape],
               If[TrueQ[normalize],
                 1 / (2 * sig * Gamma[(shape + 1) / shape] * shape ^ (1 / shape)),
                 1]]],
           RuntimeOptions -> {"Speed", "EvaluateSymbolically" -> False},
           Parallelization -> True],
        {{2,0} -> With,
         {0}   -> Compile}],
       ReplacePart[ (* Gradient function *)
         Hold[
           {{X, _Real, 2}},
           Evaluate @ Hold[
             {U = (x0 - X)},
             Evaluate @ Hold[
               {distance = Sqrt @ Total[U^2]},
               Evaluate @ Hold @ Evaluate[
                 Times[
                   weight * Exp[-(distance / sig)^shape / shape],
                   (* the extra distance in the denominator normalizes U *)
                   (distance / sig)^(shape - 1) / (sig * distance),
                   If[TrueQ[normalize],
                     1 / (2 * sig * Gamma[(shape + 1) / shape] * shape ^ (1 / shape)),
                     1]]]]],
           RuntimeOptions -> {"Speed", "EvaluateSymbolically" -> False},
           Parallelization -> True],
         {{2,2,2,0} -> Function[U * ConstantArray[#, dims]],
          {2,2,0}   -> With,
          {2,0}     -> With,
          {0}       -> Compile}],
       ReplacePart[ (* Hessian function *)
         Hold[
           {{X, _Real, 2}},
           Evaluate @ Hold[ (* part {2,0}: With *)
             {distance = Sqrt @ Total[(x0 - X)^2]},
             Evaluate @ Hold[ (* part {2,2,0}: With *)
               {F = Times[
                  -weight * Exp[-(distance / sig)^shape / shape],
                  If[TrueQ[normalize],
                    1 / (2 * sig * Gamma[(shape + 1) / shape] * shape ^ (1 / shape)),
                    1]]},
               Evaluate @ Hold[ (* part {2,2,2,0}: Times *)
                 (distance / sig)^shape  / distance^4,
                 ConstantArray[F, {dims, dims}],
                 Subtract[
                   Times[
                     ConstantArray[(2 - shape + (distance/sig)^2), {dims, dims}],
                     Outer[Times, U, U, 1, 1]],
                   Table[ConstantArray[If[i == j, distance^2, 0], dims]]]]]],
           RuntimeOptions -> {"Speed", "EvaluateSymbolically" -> False},
           Parallelization -> True],
         {{2,2,2,0} -> Times,
          {2,2,0}   -> With,
          {2,0}     -> With,
          {0}         -> Compile}]}]]];
Protect[GaussianPotentialWellParseGaussian];

GaussianPotentialWellParseSpec[mesh_, Rule[idcsArg_, gaussArg_]] := Check[
  {With[
     {CheckIndices = Replace @ {
        id_Integer :> {id},
        ids:{_Integer..} :> ids,
        _ :> Message[GaussianPotentialWell::badarg, "index "<>ToString[idcsArg]<>" not found"]}},
     Which[
       IntegerQ[idcsArg], CheckIndices @ VertexIndex[mesh, idcsArg],
       ArrayQ[idcsArg, 1, IntegerQ], CheckIndices @ VertexIndex[mesh, idcsArg],
       All, Range[VertexCount[mesh]],
       True, Message[
         GaussianPotentialWell::badarg,
         "Gaussian well spec's indices must be an integer or a list of integers"]]],
   If[ListQ[gaussArg],
     Apply[GaussianPotentialWellParseGaussian, gaussArg],
     Message[
       GaussianPotentialWell::badarg,
       "Gaussian well's spec must be a list of {x0, sigma} (possibly with options)"]]},
  $Failed];
Protect[GaussianPotentialWellParseSpec];

ParseGaussianPotentialWells[mesh_?CorticalObjectQ, spec_] := Check[
  Which[
    Head[spec] === Rule, {GaussianPotentialWellParseSpec[mesh, spec]},
    ArrayQ[spec, 1, Head[#] === Rule&], Map[GaussianPotentialWellParseSpec[mesh, #]&, spec],
    True, Message[GaussianPotentialWell::badarg, "spec must be a rule or list of rules"]],
  $Failed];
Protect[ParseGaussianPotentialWells];

(* #CalculateGaussianPotential ********************************************************************)
CalculateGaussianPotential[{indices_, fns_}, Xt_] := Total @ MapThread[
  Function[{idcs, f}, Total[(f[[1]])[Xt[[All, idcs]]]]],
  {indices, fns}];
Protect[CalculateGaussianPotential];

(* #CalculateGaussianGradient *********************************************************************)
CalculateGaussianGradient[{indices_, fns_}, Xt_] := Module[
  {dXt = ConstantArray[0.0, Dimensions[Xt]]},
  MapThread[
    Function[{idcs, f}, dXt[[All, idcs]] += (f[[2]])[Xt[[All, idcs]]]],
    {indices, fns}];
  dXt];
Protect[CalculateGaussianGradient];

(* #CalculateGaussianHessian **********************************************************************)
CalculateGaussianHessian[{indices_, fns_}, Xt_] := With[
  {dims = Length[Xt],
   n = Length[Xt[[1]]]},
  SparseArray[
    Join @@ #[[1]] -> Join @@ #[[2]]& @ Transpose @ MapThread[
       Function @ With[
         {idcs = #1, hessians = #2[[3]][Xt[[All, #1]]]},
         {Transpose[
            Flatten /@ {
              Table[idcs + i*n, {i,0,dims-1}, {dims-1}],
              Table[idcs + i*n, {dims-1}, {i,0,dims-1}]}],
          Flatten[hessians, 2]}],
      {indices, fns}],
    {n*dims, n*dims},
    0]];
Protect[CalculateGaussianGradient];

(* These functions are helpers for the Harmonic potential well functions *)
Options[HarmonicPotentialWellParseHarmonic] = {
  "Width" -> 1,
  "Shape" -> 2,
  "Weight" -> 1};
HarmonicPotentialWellParseHarmonic[x0_, OptionsPattern[]] := With[
  {width = OptionValue["Width"],
   shape = OptionValue["Shape"],
   weight = OptionValue["Weight"]},
  With[
    {const = weight / width^shape},
    Block[
      {X, U, distance},
      {(* Potential function *)
       Compile[
         {{X, _Real, 2}},
         const / shape * Sqrt[Total[(X - x0)^2]]^shape,
         RuntimeOptions -> {"Speed", "EvaluateSymbolically" -> False},
         Parallelization -> True],
       (* Gradient function *)
       Compile[
         {{X, _Real, 2}},
         With[
           {U = (X - x0)},
           With[
             {distances = (# + (1 - Unitize[#]))& @ Chop @ Sqrt @ Total[U^2]},
             With[
               {magnitudes = const * distances^(shape - 1),
                Unorm = U / ConstantArray[distances, Length[U]]},
               Unorm * ConstantArray[magnitudes, Length[U]]]]],
         RuntimeOptions -> {"Speed", "EvaluateSymbolically" -> False},
         Parallelization -> True],
       (* Hessian function *)
       Compile[
         {{X, _Real, 2}},
         With[
           {U = X - x0,
            dims = Length[X]},
           With[
             {r = (# + (1 - Unitize[#]))& @ Chop @ Sqrt @ Total[U^2]},
             With[
               {mag = ConstantArray[const / r^3, {dims, dims}],
                rdiag = With[
                  {r2 = r^2, zeros = ConstantArray[0, Length@First[X]]},
                  Table[If[i == j, r2, zeros], {i,1,dims}, {j,1,dims}]]},
               rdiag - Outer[Times, U, U, 1, 1]]]],
         RuntimeOptions -> {"Speed", "EvaluateSymbolically" -> False},
         Parallelization -> True]}]]];
Protect[HarmonicPotentialWellParseHarmonic];

HarmonicPotentialWellParseSpec[mesh_, Rule[idcsArg_, harmonicArg_]] := Check[
  {Which[
     IntegerQ[idcsArg], {VertexIndex[mesh, idcsArg]},
     ArrayQ[idcsArg, 1, IntegerQ], {VertexIndex[mesh, idcsArg]},
     All, Range[VertexCount[mesh]],
     True, Message[
       HarmonicPotentialWell::badarg,
       "Harmonic well spec's indices must be an integer or a list of integers"]],
   Which[
     ArrayQ[harmonicArg, 1, NumericQ], HarmonicPotentialWellParseHarmonic[harmonicArg],
     ListQ[harmonicArg], Apply[HarmonicPotentialWellParseHarmonic, harmonicArg],
     True, Message[
       HarmonicPotentialWell::badarg,
       "Harmonic well's spec must be a list x0 or a list of {x0} (possibly with options)"]]},
  $Failed];
Protect[HarmonicPotentialWellParseSpec];

ParseHarmonicPotentialWells[mesh_?CorticalObjectQ, spec_] := Check[
  Which[
    Head[spec] === Rule, {HarmonicPotentialWellParseSpec[mesh, spec]},
    ArrayQ[spec, 1, Head[#] === Rule&], Map[HarmonicPotentialWellParseSpec[mesh, #]&, spec],
    True, Message[HarmonicPotentialWell::badarg, "spec must be a rule or list of rules"]],
  $Failed];
Protect[ParseHarmonicPotentialWells];

(* #CalculateHarmonicPotential ********************************************************************)
CalculateHarmonicPotential[{indices_, fns_}, Xt_] := Total @ MapThread[
  Function[{idcs, f}, Total[(f[[1]])[Xt[[All, idcs]]]]],
  {indices, fns}];
Protect[CalculateHarmonicPotential];

(* #CalculateHarmonicGradient *********************************************************************)
CalculateHarmonicGradient[{indices_, fns_}, Xt_] := Module[
  {dXt = ConstantArray[0.0, Dimensions[Xt]]},
  MapThread[
    Function[{idcs, f}, dXt[[All, idcs]] += (f[[2]])[Xt[[All, idcs]]]],
    {indices, fns}];
  dXt];
Protect[CalculateHarmonicGradient];




(* ============================================================================================== *)
(* ====================================== Public Functions ====================================== *)
(* ============================================================================================== *)

(* #HarmonicEdgePotential *************************************************************************)
Options[HarmonicEdgePotential] = {MetaInformation -> {}};
HarmonicEdgePotential[mesh_?CorticalObjectQ] := With[
  {X0 = VertexCoordinatesTr[mesh],
   D0 = EdgeLengths[mesh],
   E = VertexIndex[mesh, EdgePairsTr[mesh]],
   m = EdgeCount[mesh],
   n = VertexCount[mesh],
   dims = If[CorticalMeshQ[mesh], 3, 2]},
  Module[
    {hessEdgeIdcs, hessEdgeToVertexMtcs},
    hessEdgeIdcs := hessEdgeIdcs = With[
      {rm = Range[m]},
      Transpose[
        {Flatten @ Table[rm + i*m, {i, 0, 2*dims - 1}, {2*dims}],
         Flatten @ Table[rm + i*m, {2*dims}, {i, 0, 2*dims - 1}]}]];
    hessEdgeToVertexMtcs := hessEdgeToVertexMtcs = {Transpose[#], #}& @ SparseArray[
      Rule[
        Transpose[
          {hessEdgeIdcs[[All, 1]],
           Join @@ Join @@ Table[
             If[i < dims, E[[1]], E[[2]]] + Mod[i, dims]*n, 
             {2*dims},
             {i, 0, 2*dims - 1}]}],
        1],
      {m*dims*2, n*dims},
      0];
    CorticalPotentialFunction[
      {Function[0.5 / m * Total[(EdgeLengthsTr[mesh, #] - D0)^2]],
       (* Gradient is easy: along the axis between the points *)
       Function@With[
         {x1 = #[[All, E[[1]]]], 
          x2 = #[[All, E[[2]]]]},
         With[
           {dX = CalculateDistanceGradient[x1, x2],
            r = CalculateDistance[x1, x2]},
           SumOverEdgesDirectedTr[
             mesh,
             ConstantArray[r - D0, dims] * dX[[1]] / m]]],
       (* Hessian is trickier *)
       Function@With[
         {x1 = #[[All, E[[1]]]], 
          x2 = #[[All, E[[2]]]]},
         With[
           {d2X = CalculateDistanceHessian[x1, x2],
            dX = Join @@ CalculateDistanceGradient[x1, x2],
            r = CalculateDistance[x1, x2]},
           With[
             {hess = Plus[
                d2X * ConstantArray[r, {2*dims, 2*dims}],
                Outer[Times, dX, dX, 1, 1]]},
             Dot[
               hessEdgeToVertexMtcs[[1]],
               SparseArray[hessEdgeIdcs -> Flatten[hess]],
               hessEdgeToVertexMtcs[[2]]]]]]},
      Print -> Subscript[Style["\[GothicCapitalH]",Bold], Row[{"Edges",",",Length@X0}]],
      CorticalMesh -> mesh,
      MetaInformation -> OptionValue[MetaInformation]]]];
Protect[HarmonicEdgePotential];

(* #VDWEdgePotential ******************************************************************************)
Options[VDWEdgePotential] = {MetaInformation -> {}, Order -> 2};
VDWEdgePotential[mesh_?CorticalObjectQ, OptionsPattern[]] := With[
  {X0 = VertexCoordinatesTr[mesh],
   R0 = EdgeLengths[mesh],
   E = VertexIndex[mesh, EdgePairsTr[mesh]],
   m = EdgeCount[mesh],
   n = VertexCount[mesh],
   dims = If[CorticalMeshQ[mesh], 3, 2],
   ord = OptionValue[Order] // Function[
     If[!NumericQ[#] || # <= 1, 
       Message[VDWEdgePotential::badarg, "VDWEdgePotential Order must be > 1"],
       #]]},
  With[
    {c1 = (ord - 1)/ord,
     c2 = ((ord - 1)/ord)^ord / (ord - 1)},
    CorticalPotentialFunction[
      {(* Potential... *)
       Function@With[
         {R0overR = c1 * R0 / EdgeLengthsTr[mesh, #]},
         Total[(R0overR - 1) * R0overR^(ord - 1) + c2] / m],
       (* Gradient... *)
       Function@With[
         {R = EdgeLengthsTr[mesh, #]},
         With[
           {dR = (ord * (R - R0) * (((ord - 1)*R0)/(ord*R))^ord)/(R*R0),
            dX = CalculateDistanceGradients[#[[All, E[[1]]]], #[[All, E[[2]]]]]},
           SumOverEdgesDirectedTr[
             mesh,
             ConstantArray[dR / m, dims] * dX[[1]] / m]]],
       (* No Hessian yet *)
       None},
      Print -> Subscript[Style["\[GothicCapitalV]",Bold], Row[{"Edges",",",Length@X0}]],
      CorticalMesh -> mesh,
      MetaInformation -> OptionValue[MetaInformation]]]];
Protect[VDWEdgePotential];

(* #HarmonicAnglePotential ************************************************************************)
Options[HarmonicAnglePotential] = {MetaInformation -> {}};
HarmonicAnglePotential[mesh_?CorticalObjectQ, OptionsPattern[]] := With[
  {X0 = VertexCoordinatesTr[mesh],
   Ft = VertexIndex[mesh, FaceListTr[mesh]],
   A0 = FaceAnglesTr[mesh],
   n = 3 * FaceCount[mesh],
   dims = If[CorticalMeshQ[mesh], 3, 2],
   pfun = CalculateHarmonicAnglePotential,
   gfun = CalculateHarmonicAngleGradient},
  With[
    {const = 1.0 / n},
    CorticalPotentialFunction[
      {Function@With[
         {X = #},
         With[
           {corners0 = X[[All, #]]& /@ Ft},
           const * Sum[
             pfun @@ Append[RotateLeft[corners0, i], A0[[i+1]]],
             {i, 0, 2}]]],
       Function@With[
         {X = #},
         With[
           {corners0 = X[[All, #]]& /@ Ft},
           (* corners0: {v1, v2, v3} (3xdxn); vi: {x, y, z} (3 x n) or {x, y} (2 x n) *)
           With[
             {facesGrad = Sum[
                RotateRight[
                  gfun @@ Append[RotateLeft[corners0, i], A0[[i + 1]]],
                  i],
                {i, 0, 2}]},
             (* facesGrad: same format as corners0 *)
             const * Table[SumOverFaceVerticesTr[mesh, facesGrad[[All, k]]], {k, 1, dims}]]]],
       (* No hessian yet *)
       None},
      Print -> Subscript[Style["\[GothicCapitalH]",Bold], Row[{"Angles",",",Length@X0}]],
      CorticalMesh -> mesh,
      MetaInformation -> OptionValue[MetaInformation]]]];
Protect[HarmonicAnglePotential];

(* #VWDAnglePotential *****************************************************************************)
Options[VDWAnglePotential] = {MetaInformation -> {}};
VDWAnglePotential[mesh_?CorticalObjectQ, OptionsPattern[]] := With[
  {X0 = VertexCoordinatesTr[mesh],
   Ft = VertexIndex[mesh, FaceListTr[mesh]],
   A0 = FaceAnglesTr[mesh],
   n = 3 * FaceCount[mesh],
   dims = If[CorticalMeshQ[mesh], 3, 2],
   pfun = CalculateVDWAnglePotential,
   gfun = CalculateVDWAngleGradient},
  With[
    {const = 1.0 / n},
    CorticalPotentialFunction[
      {Function@With[
         {X = #},
         With[
           {corners0 = X[[All, #]]& /@ Ft},
           const * Sum[
             pfun @@ Append[RotateLeft[corners0, i], A0[[i+1]]],
             {i, 0, 2}]]],
       Function@With[
         {X = #},
         With[
           {corners0 = X[[All, #]]& /@ Ft},
           (* corners0: {v1, v2, v3} (3xdxn); vi: {x, y, z} (3 x n) or {x, y} (2 x n) *)
           With[
             {facesGrad = Sum[
                RotateRight[
                  gfun @@ Append[RotateLeft[corners0, i], A0[[i + 1]]],
                  i],
                {i, 0, 2}]},
             (* facesGrad: same format as corners0 *)
             const * Table[SumOverFaceVerticesTr[mesh, facesGrad[[All, k]]], {k, 1, dims}]]]],
       (* No hessian yet *)
       None},
      Print -> Subscript[Style["\[GothicCapitalV]",Bold], Row[{"Angles",",",Length@X0}]],
      CorticalMesh -> mesh,
      MetaInformation -> OptionValue[MetaInformation]]]];
Protect[VDWAnglePotential];


(* #GaussianPotentialWell *************************************************************************)
Options[GaussianPotentialWell] = {MetaInformation -> {}};
GaussianPotentialWell[mesh_?CorticalObjectQ, spec_] := Check[
  With[
    {wells = Transpose @ ParseGaussianPotentialWells[mesh, spec]},
    CorticalPotentialFunction[
      {Function@CalculateGaussianPotential[wells, #],
       Function@CalculateGaussianGradient[wells, #],
       Function@CalculateGaussianHessian[wells, #]},
      Print -> Subscript[Style["\[GothicCapitalC]",Bold], Row[{"Angles",",",Length@X0}]],
      CorticalMesh -> mesh,
      MetaInformation -> OptionValue[MetaInformation]]],
  $Failed];
Protect[GaussianPotentialWell];

(* #HarmonicPotentialWell *************************************************************************)
Options[HarmonicPotentialWell] = {
  MetaInformation -> {}, 
  Print -> Automatic};
HarmonicPotentialWell[mesh_?CorticalObjectQ, spec_, OptionsPattern[]] := With[
  {wells = Transpose @ ParseHarmonicPotentialWells[mesh, spec],
   dims = If[CorticalMeshQ[mesh], 3, 2]},
  CorticalPotentialFunction[
    {Function@CalculateHarmonicPotential[wells, #],
     Function@CalculateHarmonicGradient[wells, #], 
     None},
    MetaInformation -> OptionValue[MetaInformation],
    CorticalMesh -> mesh,
    Print -> Replace[
      OptionValue[Print],
      Automatic :> Subscript[
        Style["\[GothicCapitalH]",Bold], 
        Row[{"Well","<",Length@wells,">",",",dims}]]]]];
Protect[HarmonicPotentialWell];

(* #RegionDistancePotential ***********************************************************************)
Options[RegionDistancePotential] = {
  Print -> Style["\[GothicCapitalR]", FontWeight -> Bold],
  MetaInformation -> {},
  VertexWeight -> Automatic};
RegionDistancePotential[mesh_?CorticalObjectQ, reg_?RegionQ, {F_, G_}, OptionsPattern[]] := With[
  {weight = Replace[
     OptionValue[VertexWeight],
     {list_ /; ArrayQ[list, 1, NumericQ] && Length[list] == VertexCount[mesh] :> list,
      list:{(_Integer -> _?NumericQ)..} /; Length[list] == VertexCount[mesh] :> SparseArray[
        VertexIndex[mesh, list[[All,1]]] -> list[[All, 2]], 
        VertexCount[mesh]],
      s_ /; ArrayQ[VertexPropertyValues[mesh, s], 1, NumericQ] :> VertexPropertyValues[mesh, s],
      Automatic :> ConstantArray[1, VertexCount[mesh]],
      _ :> Message[RegionDistancePotential::badarg, "Unrecognized VertexWeight option"]}]},
  With[
    {idcs = Complement[Range@Length[weight], Indices[Chop[weight], 0]],
     distFn = RegionDistance[reg],
     nearFn = RegionNearest[reg],
     W = Total[Abs@weight]},
    With[
      {w = Chop[weight], absw = Abs@Chop[weight[[idcs]]]},
      CorticalPotentialFunction[
        {Function[Dot[F @ distFn @ Transpose @ #[[All, idcs]], absw] / W],
         Function@With[
           {nears = nearFn @ Transpose @ #},
           With[
             {dX = # - Transpose[nears]},
             With[
               {dists = Chop@ColumnNorms[dX]},
               With[
                 {nulls = Unitize[dists]},
                 dX * ConstantArray[
                   nulls * w * G[dists] / (W * ((1 - nulls) + dists)),
                   Length[#]]]]]],
         None},
        Print -> OptionValue[Print],
        MetaInformation -> OptionValue[MetaInformation],
        CorticalMesh -> mesh]]]];
Protect[RegionDistancePotential];

(* #SignedRegionDistancePotential *****************************************************************)
Options[SignedRegionDistancePotential] = {
  Print -> Subscript[Style["\[GothicCapitalR]", FontWeight -> Bold], "Signed"],
  MetaInformation -> {},
  VertexWeight -> Automatic};
SignedRegionDistancePotential[mesh_?CorticalObjectQ,
                              reg_?BoundaryMeshRegionQ,
                              {F_, G_},
                              OptionsPattern[]] := With[
  {weight = Replace[
     OptionValue[VertexWeight],
     {list_ /; ArrayQ[list, 1, NumericQ] && Length[list] == VertexCount[mesh] :> list,
      list:{(_Integer -> _?NumericQ)..} /; Length[list] == VertexCount[mesh] :> SparseArray[
        VertexIndex[mesh, list[[All,1]]] -> list[[All, 2]], 
        VertexCount[mesh]],
      s_ /; ArrayQ[VertexPropertyValues[mesh, s], 1, Numeric] :> VertexPropertyValues[mesh, s],
      Automatic :> ConstantArray[1, VertexCount[mesh]],
      _ :> Message[SignedRegionDistancePotential::badarg, "Unrecognized VertexWeight option"]}]},
  With[
    {idcs = Indices[weight, Except[0|0.0]],
     distFn = SignedRegionDistance[reg],
     nearFn = RegionNearest@MeshRegion[
       MeshCoordinates[reg], 
       MeshCells[reg, If[CorticalMapQ[mesh], 1, 2]]],
     W = Total[Abs@weight]},
    With[
      {w = weight[[idcs]], absw = Abs[weight[[idcs]]]},
      CorticalPotentialFunction[
        {Function[Dot[F @ distFn @ Transpose @ #[[All, idcs]], absw] / W],
         Function@With[
           {nears = nearFn @ Transpose[#]},
           With[
             {dX = # - Transpose[nears]},
             With[
               {dists = distFn[Transpose @ #[[All, idcs]]]},
               dX * ConstantArray[
                 w * G[dists] / (# + (1 - Unitize[#]))&@(W * dists),
                 Length[#]]]]],
         None},
        Print -> OptionValue[Print],
        MetaInformation -> OptionValue[MetaInformation],
        CorticalMesh -> mesh]]]];
Protect[SignedRegionDistancePotential];

(* #HarmonicPerimeterPotential ********************************************************************)
Options[HarmonicPerimeterPotential] = {MetaInformation -> {}};
HarmonicPerimeterPotential[map_?CorticalMapQ, OptionsPattern[]] := With[
  {perimeter = MapBoundaryVertexList[map],
   X0 = VertexCoordinates[map]},
  HarmonicPotentialWell[
    map,
    Thread[perimeter -> X0[[VertexIndex[map, perimeter]]]],
    Print -> Subscript[Style["\[GothicCapitalH]", FontWeight -> Bold], "Perimeter"],
    MetaInformation -> OptionValue[MetaInformation]]];
Protect[HarmonicPerimeterPotential];

(* #MeshPotentialField ****************************************************************************)
Options[MeshPotentialField] = {
  Scale -> 1.0,
  ReferenceCoordinates -> Automatic,
  Order -> Automatic,
  StandardDeviation -> 1.0,
  Min -> 0.0,
  Max -> N[Pi],
  EdgePotentialScale -> Automatic,
  AnglePotentialScale -> Automatic};
DefineImmutable[
  MeshPotentialField[mesh_?CorticalObjectQ,
                     dir:(_String|(_String -> {_List, _List})|_List), 
                     shape_String,
                     opts:OptionsPattern[]
                    ] :> field,
  {SourceMesh[field] -> mesh,
   PotentialFieldQ[field] -> True,
   Direction[field] -> Switch[
     ToLowerCase@Which[
       StringQ[dir], dir, 
       Head[dir] === Rule, First[dir],
       VectorQ[dir, PotentialFieldQ], "sum",
       True, Message[
         PotentialField::badarg,
         "direction must be a string, rule, or list of potential fields"]],
     "edges"|"edge"|"edgelength", "Edges",
     "angles"|"angle", "Angles",
     "anchor"|"anchors", "Anchors",
     "perimeter", "Perimeter",
     "mesh"|"standard", "Mesh",
     "sum", "Sum",
     _, Message[PotentialField::badarg, "unrecognized potential direction"]],
   Shape[field] -> If[Direction[field] == "Mesh" || Direction[field] == "Sum",
     Automatic,
     Switch[
       ToLowerCase[shape],
       "harmonic", "Harmonic",
       "gaussian", "Gaussian",
       "lennardjones"|"lj"|"lennard-jones", "LennardJones",
       "well"|"infinitewell"|"infinite-well", "InfiniteWell",
       _, Message[PotentialField::badarg, "unrecognized potential shape"]]],
   Options[field] -> With[
     {d = Direction[field],
      sh = Shape[field]},
     Join[
       {Scale -> Replace[
          OptionValue[Scale],
          {s_?NumericQ :> s,
           s_ /; VectorQ[s, NumericQ] :> s,
           _ :> Message[
             PotentialField::badarg,
             "Scale must be a numeric quantity or a vector of numeric quantities"]}],
        ReferenceCoordinates -> Replace[
          OptionValue[ReferenceCoordinates],
          {Automatic :> VertexCoordinatesTr[mesh],
           X0_ /; And[MatrixQ[X0, NumericQ],
                      Dimensions[X0] == Dimensions@VertexCoordinatesTr[mesh]] :> X0,
           X0_ /; And[MatrixQ[X0, NumericQ],
                      Dimensions[X0] == Dimensions@VertexCoordinates[mesh]] :> Transpose[X0],
           _ :> Message[
             PotentialField::badarg,
             "ReferenceCoordinates was not Automatic or validly-sized matrix"]}]},
       If[d == "Anchors",
         With[
           {vtcs = dir[[2,1]],
            pts = dir[[2,2]]},
           Which[
             !VectorQ[vtcs, IntegerQ], Message[PotentialField::badarg, "incorrect anchor specification"],
             !MatrixQ[pts, NumericQ], Message[PotentialField::badarg, "Points must be a numeric matrix"],
             Length[pts] != Length[vtcs] && Length@First[pts] != Length[vtcs], Message[
               Potentia::badarg,
               "Anchors points specification must be the same size as the number of vertices"],
             And[Length[pts] != Length@VertexCoordinatesTr[mesh],
                 Length@First[pts] != Length@VertexCoordinatesTr[mesh]], Message[
               Potentia::badarg,
               "Anchors points specification must have the same dimensionality as the mesh"],
             True, {
               VertexList -> vtcs,
               Indices -> VertexIndex[mesh, vtcs],
               Points -> If[Length[pts] != Length[vtcs], pts, Transpose[pts]]}]],
         {}],
       {Order -> Replace[
          OptionValue[Order],
          {q_?NumericQ /; q > 0 :> q,
           q_ /; VectorQ[q, NumericQ] :> q,
           Automatic :> If[StringQ[sh] && sh == "InfiniteWell", 0.5, 2.0],
           _ :> Message[PotentialField::badarg, "Invalid order parameter"]}]},
       If[StringQ[sh] && sh == "InfiniteWell",
         With[
           {min = OptionValue[Min], max = OptionValue[Max]},
           Which[
             !NumericQ[min] && !VectorQ[min, NumericQ], Message[
               PotentialField::badarg,
               "Min must be a numeric quantity"],
             !NumericQ[max] && !VectorQ[max, NumericQ], Message[
               PotentialField::badarg, 
               "Max must be a numeric quantity"],
             True, {Min -> min, Max -> max}]],
         {}],
       If[StringQ[sh] && sh == "Gaussian",
         {StandardDeviation -> Replace[
            OptionValue[StandardDeviation],
            {s_?NumericQ /; s > 0 :> s,
             _ :> Message[PotentialField::badarg, "standard deviation must be a positive number"]}]},
         {}],
       If[d == "Mesh",
         With[
           {Se = Replace[OptionValue[EdgePotentialScale], Automatic -> 500],
            Sa = Replace[OptionValue[AnglePotentialScale], Automatic -> 1]},
           If[!NumericQ[Se] || Se <= 0,
             Message[PotentialField::badarg, "EdgePotentialScale must be a number > 0"]];
           If[!NumericQ[Sa] || Sa <= 0,
             Message[PotentialField::badarg, "AnglePotentialScale must be a number > 0"]];
           {EdgePotentialScale -> Se, AnglePotentialScale -> Sa}],
         {}]]],
   Options[field, name_] := Replace[name, Append[Options[field], name :> $Failed]],
   JavaObject[field] -> With[
     {ops = Options[field],
      faces = VertexIndex[mesh, FaceListTr[mesh]] - 1,
      X0 = ReferenceCoordinates /. Options[field]},
     Switch[
       {Direction[field], Shape[field]},
       {"Sum", _}, With[
          {PE = nben`mesh`registration`Fields`newSum[]},
          Do[PE@addField[JavaObject[f]], {f, dir}];
          PE],
       {"Edges", "Harmonic"}, nben`mesh`registration`Fields`newHarmonicEdgePotential[
          N[Scale /. ops], 
          N[Order /. ops],
          faces, X0],
       {"Edges", "Gaussian"}, nben`mesh`registration`Fields`newGaussianEdgePotential[
          N[Scale /. ops], 
          N[StandardDeviation /. ops],
          N[Order /. ops],
          faces, X0],
       {"Edges", "LennardJones"}, nben`mesh`registration`Fields`newLJEdgePotential[
          N[Scale /. ops], 
          N[Order /. ops],
          faces, X0],
       {"Edges", "InfiniteWell"}, nben`mesh`registration`Fields`newWellEdgePotential[
          N[Scale /. ops],
          N[Order /. ops],
          N[Min /. ops],
          N[Max /. ops],
          faces, X0],
       {"Angles", "Harmonic"}, nben`mesh`registration`Fields`newHarmonicAnglePotential[
          N[Scale /. ops], 
          N[Order /. ops],
          faces, X0],
       {"Angles", "Gaussian"}, nben`mesh`registration`Fields`newGaussianAnglePotential[
          N[Scale /. ops], 
          N[StandardDeviation /. ops],
          N[Order /. ops],
          faces, X0],
       {"Angles", "LennardJones"}, nben`mesh`registration`Fields`newLJAnglePotential[
          N[Scale /. ops], 
          N[Order /. ops],
          faces, X0],
       {"Angles", "InfiniteWell"}, nben`mesh`registration`Fields`newWellAnglePotential[
          N[Scale /. ops],
          N[Order /. ops],
          N[Min /. ops],
          N[Max /. ops],
          faces, X0],
       {"Anchors", "Harmonic"}, nben`mesh`registration`Fields`newHarmonicAnchorPotential[
          N[Scale /. ops], 
          N[Order /. ops],
          (Indices /. ops) - 1,
          (Points /. ops),
          X0],
       {"Anchors", "Gaussian"}, nben`mesh`registration`Fields`newGaussianAnchorPotential@@(Global`dbg=List[
          N[Scale /. ops], 
          N[StandardDeviation /. ops],
          N[Order /. ops],
          (Indices /. ops) - 1,
          (Points /. ops),
          X0]),
       {"Anchors", "LennardJones"}, nben`mesh`registration`Fields`newLJAnchorPotential[
          N[Scale /. ops], 
          N[Order /. ops],
          (Indices /. ops) - 1,
          (Points /. ops),
          X0],
       {"Anchors", "InfiniteWell"}, nben`mesh`registration`Fields`newWellAnchorPotential[
          N[Scale /. ops],
          N[Order /. ops],
          N[Min /. ops],
          N[Max /. ops],
          (Indices /. ops) - 1,
          (Points /. ops),
          X0],
       {"Perimeter", "Harmonic"}, nben`mesh`registration`Fields`newHarmonicPerimeterPotential[
          N[Scale /. ops], 
          N[Order /. ops],
          faces, X0],
       {"Perimeter", "Gaussian"}, nben`mesh`registration`Fields`newGaussianPerimeterPotential[
          N[Scale /. ops], 
          N[StandardDeviation /. ops],
          N[Order /. ops],
          faces, X0],
       {"Perimeter", "LennardJones"}, nben`mesh`registration`Fields`newLJPerimeterPotential[
          N[Scale /. ops], 
          N[Order /. ops],
          faces, X0],
       {"Perimeter", "InfiniteWell"}, nben`mesh`registration`Fields`newWellPerimeterPotential[
          N[Scale /. ops],
          N[Order /. ops],
          N[Min /. ops],
          N[Max /. ops],
          faces, X0],
       {"Mesh", _}, nben`mesh`registration`Fields`newStandardMeshPotential[
          N[EdgePotentialScale /. ops],
          N[AnglePotentialScale /. ops],
          faces, X0]]]},
  SetSafe -> True,
  Symbol -> MeshPotentialField];
Unprotect[MeshPotentialField];

PotentialFieldQ[_] := False;

(* MakeBoxes for potential fields... *)
MakeBoxes[p_MeshPotentialField, form_] := MakeBoxes[#, form]&@With[
  {dir = Direction[p], sh = Shape[p]},
  With[
    {name = Which[
       dir === "Mesh", "StandardMeshPotentialField",
       dir === "Sum", "PotenfialFieldSum", 
       True, Subscript[dir <> "PotentialField", sh]]},
    Row[{name["..."]}]]];

MeshPotentialField /: Plus[f___MeshPotentialField] := With[
  {fields = {f}},
  If[Length@Union[SourceMesh /@ fields] != 1,
    Message[PotentialField::badarg, "cannot add potential fields with different meshes"],
    MeshPotentialField[
      SourceMesh@First[fields],
      fields,
      "-"]]];

Protect[MeshPotentialField, ReferenceCoordinates, PotentialFieldQ];

(* #PotentialField ********************************************************************************)
Options[PotentialField] = Options[MeshPotentialField];
PotentialField[mesh_?CorticalObjectQ,
               dir_String, shape_String,
               opts:OptionsPattern[]] := MeshPotentialField[mesh, dir, shape, opts];
PotentialField[mesh_?CorticalObjectQ, 
               dir:(_String -> {_,_}), shape_String,
               opts:OptionsPattern[]] := MeshPotentialField[mesh, dir, shape, opts];
PotentialField[mesh_?CorticalObjectQ, 
               str_String /; StringMatchQ[ToLowerCase[str], "mesh"|"standard"],
               opts:OptionsPattern[]] := MeshPotentialField[mesh, str, "-", opts];
Protect[PotentialField, EdgePotentialField, AnglePotentialField];


(* #MeshRegister **********************************************************************************)
Options[MeshRegister] = {
  MaxSteps -> 10000,
  MaxStepSize -> 0.1,
  MaxPotentialChange -> 1.0,
  VertexCoordinates -> Automatic};
MeshRegister[field_?PotentialFieldQ, opts:OptionsPattern[]] := With[
  {steps = Replace[
     OptionValue[MaxSteps],
     {x_Integer?Positive :> x,
      _ :> Message[MeshRegister::badarg, "MaxSteps must be a postive integer"]}],
   stepsz = Replace[
     OptionValue[MaxStepSize],
     {x_?NumericQ /; Positive[x] :> x,
      _ :> Message[MeshRegister::badarg, "MaxStepSize must be a postive number"]}],
   pchange = Replace[
     OptionValue[MaxPotentialChange],
     {x_?NumericQ /; 0 <= x <= 1 :> x,
      _ :> Message[MeshRegister::badarg, "MaxPotentialChange must be a number between 0 and 1"]}],
   X = Replace[
     OptionValue[VertexCoordinates],
     {x_ /; MatrixQ[x, NumericQ] :> If[Length[x] > Length@First[x], Transpose[x], x],
      Automatic :> VertexCoordinatesTr@SourceMesh[field],
      _ :> Message[MeshRegister::badarg, "VertexCoordinates must be a numerical matrix"]}],
   mesh = SourceMesh[field],
   PE = JavaObject[field]},
  With[
    {min = JavaNew["nben.mesh.registration.Minimizer", PE, X]},
    min@step[pchange, steps, stepsz];
    Clone[mesh, VertexCoordinatesTr -> (min@getX[])]]];
MeshRegister[mesh_?CorticalObjectQ, instr_List, opts:OptionsPattern[]] := MeshRegister[
  Plus@@Table[
    PotentialField[mesh, Sequence@@q],
    {q, instr}],
  opts];
Protect[MeshRegister, MaxPotentialChange];

(* #CortexGradientPlot ****************************************************************************)
Options[CortexGradientPlot] = Join[
  {Arrowheads -> Small,
   PlotStyle -> Automatic,
   Scaled -> Automatic},
  Options[CortexPlot]];
CortexGradientPlot[mesh_?CorticalMapQ, functions_List, opts : OptionsPattern[]] := With[
  {arrowheads = OptionValue[Arrowheads],
   scale = Replace[OptionValue[Scaled], Automatic :> 0.75*Max[EdgeLengths[mesh]]],
   styles = OptionValue[PlotStyle],
   cortexPlotOpts = FilterRules[{opts}, Options[CortexPlot]],
   showOpts = FilterRules[{opts}, Options[Show]],
   X = VertexCoordinates[mesh],
   n = Length[functions]},
  Show[
    {CortexPlot[
       mesh,
       Sequence @@ cortexPlotOpts,
       EdgeRenderingFunction -> ({Thin, Line[#Coordinates, VertexColors -> #VertexColors]}&),
       FaceRenderingFunction -> None],
     With[
       {grads = Partition[Grad[#, X], 2] & /@ functions},
       With[
         {coef = scale/Max[Join @@ Map[RowNorms, grads]]},
         Graphics[
           {Arrowheads[arrowheads],
            MapIndexed[
              Function@With[
                {style = Which[
                   styles === Automatic, ColorData["Rainbow"][
                     If[n == 1, 0.5, (#2[[1]] - 1)/(n - 1)]],
                   !ListQ[styles], styles,
                   True, styles[[Mod[#2[[1]] - 1, n] + 1]]]},
                Join[
                  If[ListQ[style], style, {style}],
                  MapThread[Arrow[{#1, #1 + coef*#2}] &, {X, #1}]]],
              grads]}]]]},
    Sequence @@ showOpts]];
CortexGradientPlot[Rule[X_ /; ArrayQ[X, 2, NumericQ], mesh_?CorticalMapQ],
                   functions_List,
                   opts:OptionsPattern[]] := CortexGradientPlot[
  CorticalMap[mesh, VertexCoordinates -> X],
  functions,
  opts];
Protect[CortexGradientPlot];

(* #MapTangledQ ***********************************************************************************)
(* This compiled function does the work, but expects the triangle coordinates to be rows of the
   given matrix FX, which should have the triangles in either clockwise or counter-clockwise order;
   in other words, this is to be called by MapTangledQ privately. *)
TrianglesTangledQ = Compile[
  {{FX, _Real, 3}},
  (* The basic check here:
     We know that all triangles should be expressed in the same order (clockwise or
     counterclockwise); accordingly, all angles should be either positive or negative (when the
     smallest angle is taken in the ordered direction). *)
  With[
    {th1 = ArcTan[FX[[2, 1]] - FX[[1, 1]], FX[[2, 2]] - FX[[1, 2]]],
     th2 = ArcTan[FX[[3, 1]] - FX[[1, 1]], FX[[3, 2]] - FX[[1, 2]]]},
    With[
      {th = th2 - th1},
      (* For angles that are over pi or under -pi, we want to invert the sign; then we just check
         to make sure all signes are the same. *)
      Length@Union@Sign[Sign[Pi - Abs[th]] * th] != 1]],
  RuntimeOptions -> {"Speed", "EvaluateSymbolically" -> False}];
(* Here are the MapTangledQ wrappers: *)
MapTangledQ[map_?CorticalMapQ, Xarg_List, vtxIdcs_List] := With[
  {X = If[Length[Xarg] == 2, Xarg, Transpose[Xarg]],
   idcs = Union@Apply[Join, VertexFaceList[map][[vtxIdcs]]]},
  TrianglesTangledQ[X[[All, #]]& /@ VertexIndex[map, FaceListTr[map][[All, idcs]]]]];
MapTangledQ[map_?CorticalMapQ, Xarg_List] := With[
  {X = If[Length[Xarg] == 2, Xarg, Transpose[Xarg]]},
  TrianglesTangledQ[X[[All, #]]& /@ VertexIndex[map, FaceListTr[map]]]];
MapTangledQ[map_?CorticalMapQ, Automatic, is_] := MapTangledQ[map, VertexCoordinatesTr[map], is];
MapTangledQ[map_?CorticalMapQ, Automatic] := MapTangledQ[map, VertexCoordinatesTr[map]];
MapTangledQ[map_?CorticalMapQ] := MapTangledQ[map, VertexCoordinatesTr[map]];
Protect[MapTangledQ];

(* #MapTangles ************************************************************************************)
MapTangles[map_?CorticalMapQ, Xarg_] := With[
  {X = If[Length[Xarg] == 2, Xarg, Transpose[Xarg]],
   Ft = VertexIndex[map, FaceListTr[map]]},
  With[
    {FX = X[[All, #]]& /@ Ft},
    With[
      {th1 = ArcTan[FX[[2, 1]] - FX[[1, 1]], FX[[2, 2]] - FX[[1, 2]]],
       th2 = ArcTan[FX[[3, 1]] - FX[[1, 1]], FX[[3, 2]] - FX[[1, 2]]]},
      With[
        {th = Sign[Sign[Pi - Abs[#]] * #]&[th2 - th1]},
        With[
          {signs = Tally[th]},
          (* if signs has 1 value, there are no tangles *)
          If[Length[signs] == 1, 
            {},
            (* otherwise, we take the more common value to be the 'correct' one *)
            With[
              {bad = SortBy[signs, Last][[1 ;; -2, 1]]},
              With[
                {angIDs = Union@Indices[th, If[Length[bad] == 1, bad[[1]], Alternatives@@bad]]},
                Union@Flatten@Part[Ft, All, angIDs]]]]]]]]];                  
(*
MapTangles[map_?CorticalMapQ, X_] := With[
  {nei = VertexIndex[map, NeighborhoodList[map]]},
  Flatten@Last@Reap[
    MapThread[
      Function @ With[
        {x = #1, n = #2, k = #3},
        If[Length[n] > 2,
          With[
            {xnei = X[[#]] & /@ n},
            With[
             {ord = Ordering[ArcTan[#[[1]] - x[[1]], #[[2]] - x[[2]]] & /@ xnei]},
             With[
               {rot = RotateLeft[ord, Position[ord, 1, {1}][[1, 1]] - 1]},
               If[Range[Length[rot]] != rot || !Graphics`Mesh`InPolygonQ[xnei, x],
                 Sow[k]]]]]]],
      {X, nei, Range[Length@X]}]]];
*)
MapTangles[map_?CorticalMapQ] := MapTangles[map, VertexCoordinates[map]];
Protect[MapTangles];

(* #MapUntangle ***********************************************************************************)
Options[MapUntangle] = {MaxIterations -> 50};
MapUntangle[map_?CorticalMapQ, Xtang_?MatrixQ, OptionsPattern[]] := With[
  {nei = VertexIndex[map, NeighborhoodList[map]],
   max = OptionValue[MaxIterations]},
  NestWhile[
   Function@With[
     {tangles = Flatten[Join[#, nei[[#]] & /@ #]]& @ MapTangles[map, #],
      X = #},
     ReplacePart[X, Thread[tangles -> Map[Mean[X[[nei[[#]]]]]&, tangles]]]],
    Xtang,
    MapTangledQ[map, #]&,
    1,
    max]];
MapUntangle[map_?CorticalMapQ, opts:OptionsPattern[]] := MapUntangle[
  map,
  VertexCoordinates[map],
  opts];
Protect[MapUntangle];

(* #RegistrationChooseStepSize ********************************************************************)
$RegistrationChooseStepSizeCutoff = 0.05;
RegistrationChooseStepSize[
  potential_, (* Potential Function *)
  map_, (* The original map *)
  X0_List, (* The starting coordinates for this step *)
  PE0_, (* initial potential energy *)
  J0_ (* Initial Jacobian/Gradient *)
 ] := With[
  {nJ0 = ColumnNorms[J0]},
   With[
    {maxNormJ0 = Max[nJ0]},
    With[
      {scaledJ0 = J0 / maxNormJ0,
       scnormJ0 = nJ0 / maxNormJ0},
      (* okay, partition the list into high and low and grab the vtcs in the faces for the highs *)
      With[
        {faceIDs = Union@@Part[
             VertexFaceList[map],
             Pick[Range@Length[scnormJ0], Sign[scnormJ0 - 0.05], 1]]},
        With[
          {vtxIDs = VertexIndex[map, FaceListTr[map][[All, faceIDs]]]},
          (* we now minimize with these, starting with stepsize = 1...
             We are primarily looking for a decrease in the potential function... *)
          NestWhile[
            (* So if the step-size is two high, halve it... *)
            Function@If[Chop[#] == 0, 0, 0.5 * #],
            (* We look for the starting step size ss0 such that size ss0 does not tangle the map,
               but size 2*ss does; this should be a good place to start the search, and guarantees
               that we will never tangle the map. *)
            With[
              {grad = J0[[All, #]]& /@ vtxIDs,
               x = X0[[All, #]]& /@ vtxIDs},
              With[
                {ss0 = NestWhile[
                   (* We either increase or decrease stepsize, depending on last result *)
                   Function@With[
                     {ss = #[[1]] * If[#[[2]], 0.5, 2.0]},
                     If[Chop[ss] == 0, 0, {ss, TrianglesTangledQ[x - ss*grad]}]],
                   (* start with a stepsize of 1 *)
                   {1/maxNormJ0, TrianglesTangledQ[x - (1/maxNormJ0)*grad]},
                   (* run until one step is true and one is false *)
                   Function@Which[
                     #1 === 0, False,
                     Xor[#1[[2]], #2[[2]]], False,
                     True, True],
                   (* provide two arguments to the test *)
                   2]},
                If[ss0[[2]], 0.5*ss0[[1]], ss0[[1]]]]],
            (* okay, we basically just want to check if the function has decreased *)
            Function@If[# == 0, 
              False,
              With[
                {X1 = X0 - #*J0},
                potential[X1] >= PE0 || MapTangledQ[map, X1]]]]]]]]];
Protect[RegistrationChooseStepSize];

(* #RegistrationFrame *****************************************************************************)
Options[RegistrationFrame] = {
  Method -> "GradientDescent",
  MaxVertexChange -> 0.1,
  MinStepSize -> 10^-5,
  CorticalMap -> None};
DefineImmutable[
  RegistrationFrame[P_CorticalPotentialFunctionInstance, X0_, OptionsPattern[]] :> frame,
  {VertexCount[frame] -> Length@First[X0],
   Dimensions[frame] -> Dimensions[X0],

   PotentialFunction[frame] = P,
   VertexCoordinatesTr[frame] = If[Length[X0] < Length[X0[[1]]], X0, Transpose[X0]],
   MaxVertexChange[frame] = OptionValue[MaxVertexChange],
   MinStepSize[frame] = OptionValue[MinStepSize],
   Method[frame] = OptionValue[Method],
   CorticalMap[frame] -> Replace[
     OptionValue[CorticalMap],
     Except[None | _?CorticalMapQ] :> Message[
       RegistrationFrame::badarg,
       "CorticalMap option to RegistrationFrame must be a cortical map or None"]],
   StepNumber[frame] = 0,
   
   VertexCoordinates[frame] := Transpose@VertexCoordinatesTr[frame],
   
   Value[frame] :> PotentialFunction[frame][VertexCoordinatesTr[frame]],
   Gradient[frame] :> Partition[
     Grad[PotentialFunction[frame], VertexCoordinatesTr[frame]],
     VertexCount[frame]],
   (*
   StepSize[frame] :> With[
     {grad = Gradient[frame],
      x0 = VertexCoordinatesTr[frame],
      pe0 = Value[frame],
      pf = PotentialFunction[frame],
      minss = MinStepSize[frame],
      map = CorticalMap[frame]},
     With[
       {gradNorms = ColumnNorms[grad]},
       With[
         {maxVtxNorm = Max[gradNorms]},
         NestWhile[
           (0.5*#) &,
           MaxVertexChange[frame]/maxVtxNorm,
           If[map === None,
             # > minss && pe0 < pf[x0 - #*grad] &,
             Function@With[
               {x1 = x0 - #*grad},
               # > minss && pe0 < pf[x1] && !MapTangledQ[map, x1]]]]]]],
      *)
   StepSize[frame] :> RegistrationChooseStepSize[
     PotentialFunction[frame],
     CorticalMap[frame],
     VertexCoordinatesTr[frame],
     Value[frame],
     Gradient[frame]],
   Next[frame] := If[StepSize[frame] <= MinStepSize[frame],
     None,
     Clone[
       frame,
       VertexCoordinatesTr -> (VertexCoordinatesTr[frame] - StepSize[frame]*Gradient[frame]),
       StepNumber -> (StepNumber[frame] + 1)]],
   
   RegistrationFrameQ[frame] -> Which[
     Dimensions[VertexCoordinatesTr[frame]] != Dimensions[frame], Message[
       RegistrationFrame::badarg,
       "coordinates have the wrong dimensions: " <> ToString[Dimensions@VertexCoordinatesTr[frame]]
         <> " and " <> ToString@Dimensions[frame]],
     Method[frame] == "ConjugateGradient" || Method[frame] == "MonteCarlo", Message[
       RegistrationFrame::badarg,
       "Method option must be \"ConjucateGradient\" or \"MonteCarlo\""],
     MaxVertexChange[frame] <= 0, Message[
       RegistrationFrame::badarg,
       "MaxVertexChange must be > 0"],
     MinStepSize[frame] <= 0, Message[
       RegistrationFrame::badarg,
       "MinSteoSize must be > 0"],
     Head[PotentialFunction[frame]] =!= 
     CorticalPotentialFunctionInstance, Message[
       RegistrationFrame::badarg,
       "Invalid cortical potential function"],
     True, True]},
  SetSafe -> True,
  Symbol -> RegistrationFrameData];
MakeBoxes[frame:RegistrationFrameData[___], form_] := MakeBoxes[#]& @ With[
  {style = {
     FontSize -> 11,
     FontColor -> Gray,
     FontFamily -> "Arial",
     FontWeight -> "Thin"}},
  Row[
    {"RegistrationFrame"[
       Panel[
         Grid[
           MapThread[
             Function[{Spacer[4], Style[#1, Sequence @@ style], Spacer[2], #2, Spacer[4]}],
             {{"Step Number:", "Potential:", "Max Vertex Gradient:"},
              {StepNumber[frame], Value[frame], Max@ColumnNorms[Gradient[frame]]}}],
           Alignment -> Table[{Right, Right, Center, Left, Left}, {3}]]]]},
    BaseStyle -> Darker[Gray]]];
Protect[RegistrationFrame, RegistrationFrameData, RegistrationFrameQ, StepNumber, StepSize, 
        MinStepSize, MaxVertexChange];

(* #RegistrationTrajectory ************************************************************************)
RegistrationTrajectorySymbolLookup[sym_, k_] := With[
  {f = sym["CacheFrequency"],
   autoCacheFn = sym["AutoCache"],
   F0 = sym["InitialFrame"],
   cmesh = sym["Mesh"],
   ss0 = sym["MinStepSize"]},
  Check[
    If[k > Last[sym],
      sym[Last[sym]],
      With[
        {cache = If[ac =!= None && Mod[k, f] == 0,
           AutoCache[ac[k], $Failed],
           $Failed],
         prev = f*Floor[(k - 1)/f]},
        With[
          {rule = If[cache =!= $Failed,
             cache,
             With[
               {res = Check[
                  NestWhile[
                    Function@With[
                      {nexts = Rest@NestWhileList[Next, #, (Next[#] =!= None)&, 1, k - prev]},
                      With[
                        {ln = If[nexts == {}, None, Last[nexts]]},
                        Which[
                          nexts == {}, prev -> VertexCoordinatesTr[#],
                          !RegistrationFrameQ[ln], $Failed,
                          !ArrayQ[VertexCoordinatesTr[ln], 2, NumericQ], $Failed,
                          MapTangledQ[cmesh, VertexCoordinates[ln]], $Failed,
                          Length[nexts] < k - prev, Rule[
                            prev + Length[nexts],
                            VertexCoordinatesTr[ln]],
                          True, k -> VertexCoordinatesTr[ln]]]],
                    sym[prev],
                    Function@Which[
                      Head[#] === Rule, False,
                      # === $Failed, (
                        Message[RegistrationTrajectory::nocnv, "could not take valid step"];
                        False),
                      MaxVertexChange[#] <= minss, (
                        Message[RegistrationTrajectory::nocnv, "max stepsize reduced to 0"];
                        False),
                      True, True]],
                  $Failed]},
               If[And[ac =!= None, Mod[k, f] == 0, prev < Last[sym], 
                      res =!= $Failed, cache === $Failed],
                 AutoCache[ac[k], res],
                 res]]]},
          If[rule === $Failed,
            $Failed,
            With[
              {frame = Clone[
                 F0,
                 StepNumber -> rule[[1]],
                 VertexCoordinatesTr -> rule[[2]]]},
              If[rule[[1]] > Max[sym], sym /: Max[sym] = rule[[1]]];
              If[k > rule[[1]], sym /: Last[sym] = rule[[1]]];
              If[Mod[k, f] == 0 && k == rule[[1]], sym[k] = frame];
              frame]]]]],
    $Failed]];
Options[RegistrationTrajectory] = {
  MaxVertexChange -> 0.1,
  MinStepSize -> 10^-5,
  CacheFrequency -> 50,
  InitialVertexCoordinates -> Automatic,
  AutoCache -> None};
DefineImmutable[
  RegistrationTrajectory[
    mesh_?CorticalObjectQ,
    P_CorticalPotentialFunctionInstance,
    opts:OptionsPattern[]] :> traj,
  {MaxVertexChange[traj] = OptionValue[MaxVertexChange],
   MinStepSize[traj] = OptionValue[MinStepSize],

   CacheFrequency[traj] -> OptionValue[CacheFrequency],
   AutoCache[traj] -> With[
     {tmp = OptionValue[AutoCache]},
     Which[
       tmp === None || tmp === False, None,
       !StringQ[tmp] || !StringContainsQ[tmp, "XXX"~~"X"..], Message[
         RegistrationTrajectory::badarg,
         "AutoCache must be a string with at least 4 X's"],
       True, Function[
         Which[
           # === Normal, tmp,
           True, StringReplace[
             tmp,
             s:("XXX"~~"X"..) :> If[IntegerQ[#], 
               IntegerString[#, 10, StringLength[s]],
               ToString[#]]]]]]],
   PotentialFunction[traj] -> P,
   CorticalMesh[traj] -> mesh,
   InitialFrame[traj] -> With[
     {f0 = RegistrationFrame[
        PotentialFunction[traj],
        Replace[
          OptionValue[InitialVertexCoordinates],
          {Automatic :> VertexCoordinatesTr[mesh],
           X_ /; Dimensions[X] == Dimensions@VertexCoordinates[mesh] :> Transpose[X],
           X_ /; Dimensions[X] == Dimensions@VertexCoordinatesTr[mesh] :> X,
           _ :> Message[
             RegistrationTrajectory::badarg,
             "InitialVertexCoordinaets must be Automatic or an appropriately sized matrix"]}],
        Sequence @@ Join[
          FilterRules[{opts}, Options[RegistrationFrame]],
          If[CorticalMapQ[mesh], {CorticalMap -> mesh}, {}]]]},
     Which[
       !RegistrationFrameQ[f0], Message[
         RegistrationTrajectory::badarg,
         "arguments did not produce a valid initial frame"],
       Dimensions@VertexCoordinatesTr@f0 != Dimensions@VertexCoordinatesTr[mesh], Message[
         RegistrationTrajectory::badarg,
         "X0 must be appropriately sized for VertexCoordinatesTr[mesh]"],
       True, f0]],
   Symbol[traj] -> With[
     {f = CacheFrequency[traj],
      F0 = InitialFrame[traj],
      ac = AutoCache[traj],
      cmesh = CorticalMesh[traj],
      sym = TemporarySymbol["trajectory"],
      minss = MinStepSize[traj]},
     If[!IntegerQ[f] || f < 1,
       Message[RegistrationTrajectory::badarg, "CacheFrequency must be an integer > 0"]];
     sym /: Max[sym] = 0;
     sym /: Last[sym] = Infinity;
     sym[0] = F0;
     sym["CacheFrequency"] = f;
     sym["AutoCache"] = ac;
     sym["InitialFrame"] = F0;
     sym["Mesh"] = cmesh;
     sym["MinStepSize"] = minss;
     sym[k_Integer /; k > 0] := RegistrationTrajectorySymbolLookup[sym, k];
     sym],
   
   Frame[traj, k_Integer /; k >= 0] := Symbol[traj][k],
   Frame[traj, ks_ /; ArrayQ[ks, 1, IntegerQ[#] && # >= 0 &]] := Map[Symbol[traj], ks],
   FinalFrame[traj] :> With[
     {sym = Symbol[traj], f = CacheFrequency[traj]},
     NestWhile[
       (sym[#]; # + f)&,
       f * Floor[(Max[sym] + f) / f],
       Last[sym] === Infinity &];
     sym[Last[sym]]],

   RegistrationTrajectoryQ[traj] -> With[
     {mvc = MaxVertexChange[traj],
      mss = MinStepSize[traj],
      cfreq = CacheFrequency[traj]},
     Which[
       mvc <= 0, Message[
         RegistrationTrajectory::badarg,
         "MaxVertexChange must be > 0"],
       mss <= 0, Message[
         RegistrationTrajectory::badarg,
         "MinSteoSize must be > 0"],
       cfreq =!= None && (! IntegerQ[cfreq] || cfreq < 1), Message[
         RegistrationTrajectory::badarg,
         "CacheFrequency must be an integer > 0 or None"],
       True, True]]},
  Symbol -> RegistrationTrajectoryData,
  SetSafe -> True];
MakeBoxes[frame:RegistrationTrajectoryData[___], form_] := MakeBoxes[#]& @ With[
  {style = {
     FontSize -> 11,
     FontColor -> Gray,
     FontFamily -> "Arial",
     FontWeight -> "Thin"}},
  Row[
    {"RegistrationTrajectory"[
       Panel[
         Grid[
           MapThread[
             {Spacer[4], Style[#1, Sequence @@ style], Spacer[2], #2, Spacer[4]}&,
             {{"Dimensions:", "Potential Function:"},
              {Dimensions@VertexCoordinates@InitialFrame[frame], PotentialFunction[frame]}}],
         Alignment -> Table[{Right, Right, Center, Left, Left}, {3}]]]]},
     BaseStyle -> Darker[Gray]]];
Protect[RegistrationTrajectory, RegistrationTrajectoryData, RegistrationTrajectoryQ, CacheFrequency,
        InitialFrame, InitialVertexCoordinates];

(* #MeshRegister **********************************************************************************)

(*
Options[MeshRegister] = Join[
  Options[RegistrationTrajectory],
  {MaxIterations -> 10000}];
MeshRegister[mesh_?CorticalObjectQ,
             PF_CorticalPotentialFunctionInstance,
             opts:OptionsPattern[]] := With[
  {maxiter = Replace[
     Option},
  With[
    {},
    ]];
*)

End[];
EndPackage[];

