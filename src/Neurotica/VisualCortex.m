(* VisualCortex.m
 *
 * Functions and data structures for efficiently implementing retinotopy on the cortical surface,
 * including the Banded Double-Sech model of Schira, Tyler, Spehar, and Breakspear (2010), in
 * Mathematica.
 * References:
 *   [Retinotopic Templates]
 *   Benson NC, Butt OH, Brainard DH, Aguirre GK (2014) Correction of Distortion in Flattened 
 *     Representations of the Cortical Surface Allows Prediction of V1-V3 Functional Organization
 *     from Anatomy. PLoS Comput Biol 10(3): e1003538. doi: 10.1371/journal.pcbi.1003538
 *   [Schira Model of Retinotopy]
 *   Schira MM, Tyler CW, Spehar B, Breakspear M (2010) Modeling Magnification and Anisotropy in the
 *     Primate Foveal Confluence. PLoS Comput Biol 6(1): e1000651. doi: 10.1371/journal.pcbi.1000651
 *   [V1 Masks and Hulls]
 *   Hinds OP, Rajendran N, Polimeni JR, Augustinack JC, Wiggins G, Wald LL, Diana Rosas H,
 *     Potthast A, Schwartz EL, Fischl B (2008) Accurate prediction of V1 location from cortical
 *     folds in a surface coordinate system. Neuroimage. 2008 Feb 15;39(4):1585-99.
 *     doi: 10.1016/j.neuroimage.2007.10.033
 *
 * Copyright (C) 2015 by Noah C. Benson.
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
  "Neurotica`VisualCortex`",
  {"NDSolve`FEM`",
   "Neurotica`Global`", 
   "Neurotica`Util`",
   "Neurotica`Coordinates`",
   "Neurotica`Mesh`",
   "Neurotica`FreeSurfer`",
   "Neurotica`Registration`"}];
Unprotect["Neurotica`VisualCortex`*", "Neurotica`VisualCortex`Private`*"];
ClearAll[ "Neurotica`VisualCortex`*", "Neurotica`VisualCortex`Private`*"];

(* Complex Representation of the Visual Field  ****************************************************)
VisualAngleToComplex::usage = "VisualAngleToComplex[polarAngle, eccentricity] yields an imaginary number that represents the visual field coordinate. This complex number is always of the form r * Exp[I * t] where t is between -Pi/2 and Pi/2.";
ComplexToVisualAngle::usage = "ComplexToVisualAngle[z] yields a {polarAngle, eccentricity} pair that represents the visual angle coordinates represented by the complex number z. This complex  number should be of the form r * Exp[I * t] where t is between -Pi/2 and Pi/2.";
ComplexToPolarAngle::usage = "ComplexToPolarAngle[z] yields the polar angle value that is represented by the complex number z. This complex number should be of the form r * Exp[I * t] where t is between -Pi/2 and Pi/2.";
ComplexToEccentricity::usage = "ComplexToEccentricity[z] yields the eccentricity value that is represented by the complex number z. This complex number should be of the form r * Exp[I * t] where t is between -Pi/2 and Pi/2.";
ComplexToCoordinate::usage = "ComplexToCoordinate[z] yields {Re[z], Im[z]} for a complex z and a list of such for a list of complex z.";
CoordinateToComplex::usage = "CoordinateToComplex[{x,y}] yields x + I*y for real x and y.
CoordinateToComplex[{{x1,y1},{x2,y2}...}] yields a list of zi = xi + yi*I.
CoordinateToComplex[x, y] is equivalent to CoordinateToComplex[Thread[{x,y}]].";
CoordinateToVisualAngle::usage = "CoordinateToVisualAngle[x, y] yields the visual angle coordinate {\[Theta], \[Rho]} of the point given in the visual field by {x,y}. The polar angle is \[Theta] and is between -180 and 180 (degrees) with negative values indicating the left visual field, and the eccentricity value is \[Rho] and is between 0 and \[Infinity]. If the x and y arguments are equally sized lists, the funtion is automatically threaded over the arguments.
CoordinateToVisualAngle[{x,y}] for atomic x and y is equivalent to CoordinateToVisualAngle[x,y].
CoordinateToVisualAngle[{{x1, y1}, {x2, y2}, ...}] automatically threads the CoordinateToVisualAngle function across all given points. ";
VisualAngleToCoordinate::usage = "VisualAngleToCoordinate[\[Theta], \[Rho]] yields the cartesian coordinate {x,y} of the point in the visual field corresponding the polar angle \[Theta] and eccentricity \[Rho] where \[Theta] must be between -180 and 180 (degrees) with negative degrees representing the left visual field and the eccentricity \[Rho] should be between 0 and \[Infinity] (though negative values are accepted). If the \[Theta] and \[Rho] arguments are equally sized lists, then the result is automatically threaded over them.
VisualAngleToCoordinate[{\[Theta], \[Rho]] is equivalent to VisualAngleToCoordinate[\[Theta], \[Rho]].
VisualAngleToCoordinate[{{\[Theta]1, \[Rho]1}, {\[Theta]2, \[Rho]2}, ...}] is automatically threaded over all the individual \[Theta] and \[Rho] values.";

(* Visual Areas ***********************************************************************************)
VisualAreaQ::usage = "VisualAreaQ[area] yields true if and only if area is a valid visual area id.";
VisualAreaData::usage = "VisualAreaData[id] yields a list of data regarding the visual area whose ID is given. VisualAreaData[] yields a list of all known IDs and their data.";
VisualAreaName::usage = "VisualAreaName[id] yields the visual area name of the given visual area ID.";
VisualAreaSimplify::usage = "VisualAreaSimplify[id] yields a simplified version of id such that border IDs (e.g., V1/V2 border) are converted to the lower visual area, ventral and dorsal signatures are removed (all made dorsal), and anything outside of V1, V2, and V3 is converted to 0. Accordingly, this function always yields 0, 1, 2, or 3 for any valid visual area id.";
VisualAreaData::badarea = "An unknown visual area was given: `1`";

(* Retinotopic Templates **************************************************************************)
FSAverageSymPolarAngle::usage = "FSAverageSymPolarAngle yields the polar angle data for the fsaverage_sym subject as a field on the fsaverage_sym hemisphere.";
FSAverageSymEccentricity::usage = "FSAverageSymEccentricity yields the eccentricity data for the  fsaverage_sym subject as a field on the fsaverage_sym hemisphere.";
FSAverageSymVisualArea::usage = "FSAverageSymVisualArea yields the visual area data for the  fsaverage_sym subject as a field on the fsaverage_sym hemisphere. Each vertex will either be labeled as 1, 2, 3 (for dorsal regions), or -1, -2, -3 (for ventral regions) or None.";
FSAverageSymRetinotopy::usage = "FSAverageSymRetinotopy yields the retinotopy field for the fsaverage_sym subject. The field for each vertex is {PA, E, A} where PA is the polar angle, E is the eccentricity, and A is the area.";

PredictRetinotopy::usae = "PredictRetinotopy[subject, meshName, hemi] is equivalent to Cortex[subject, meshName, hemi] but yields a cortical surface with a prediction of the subject's retinotopy, based on the subject's fsaverage_sym-aligned spherical hemisphere. If the subject does not have a sym hemisphere, an error is raised. The retinotopy data is included in the surface's properties with the names \"PolarAngle\", \"Eccentricity\", and \"VisualAngle\".";

(* The Schira Model *******************************************************************************)
$DefaultSchiraA::usage = "The default value of the Schira model A parameter.";
$DefaultSchiraB::usage = "The default value of the Schira model B parameter.";
$DefaultSchira\[CapitalLambda]::usage = "The default value of the Schira model lambda parameter.";
$DefaultSchira\[CapitalPsi]::usage = "The default rotation of the Schira model in radians.";
$DefaultSchira\[CapitalRho]90::usage = "The default degree position for 90 degrees in the Schira model.";
$DefaultSchiraV1Size::usage = "The default size of V1 in the Schira model.";
$DefaultSchiraV2Size::usage = "The default size of V2 in the Schira model.";
$DefaultSchiraV3Size::usage = "The default size of V3 in the Schira model.";
$DefaultSchiraHV4Size::usage = "The default size of the hV4 pseudo-area in the Schira model.";
$DefaultSchiraV3ASize::usage = "The default size of V3A-like pseudo-area in the Schira model.";
$DefaultSchiraFC::usage = "The default position of the foveal confluence in the Schira model.";
$DefaultSchiraScale::usage = "The default {x,y} scale of the Schira model.";
$DefaultSchiraShear::usage = "The default shear matrix of the Schira model.";
$SchiraParameters::usage = "$SchiraParameters yields a list of the default parameters for the Schira model. The list rules are delayed so any Block redefining the default Schira parameters will affect this list's values.";

A::usage="The A parameter of the Schira model.";
B::usage="The B parameter of the Schira model.";
\[CapitalPsi]::usage="The rotation parameter of the Schira model.";
\[CapitalLambda]::usage="The rotation parameter of the Schira model.";
\[CapitalRho]90::usage = "The degree position of 90 degrees in the Schira model.";
V1Size::usage="The parameter of the Schira model that determines the size of V1.";
V2Size::usage="The parameter of the Schira model that determines the size of V2.";
V3Size::usage="The parameter of the Schira model that determines the size of V3.";
HV4Size::usage="The parameter of the Schira model that determines the size of the hV4 pseudo-area.";
V3ASize::usage="The parameter of the Schira model that determines the size of the V3A pseudo-area.";
FC::usage="The foveal confluence position parameter of the Schira model.";
Shear::usage="The shear matrix parameter of the Schira model.";

SchiraModel::usage = "SchiraModel[parameters...] yields a SchiraModelObject using the given
 parameters when possible and using the default Schira parameters otherwise.";
SchiraModel::badarg = "Bad argument(s); all arguments must be rules: `1`";
SchiraModel::badarea = "Bad area given to Schira function; areas must be 1, 2, 3, or 4.";
SchiraModelObject::usage = "A SchiraModelObject form stores the data for a Schira model.";
SchiraModelObject::badarg = "Unrecognized SchiraModelObject argument: `1`";
SchiraFunction::usage = "SchiraFunction[mdl] yields the forward tranformation function for the given Schira model mdl. This is equivalent to mdl[Function].";
SchiraFunction::badarg = "Bad argument given to Schira function: `1`";
SchiraInverse::usage = "SchiraFunction[mdl] yields the inverse tranformation function for the given Schira model mdl. This is equivalent to mdl[Inverse].";
SchiraInverse::badarg = "Bad argument given to Schira inverse function: `1`";
CorticalMapToVisualField::usage = "CorticalMapToVisualField[model, map] yields a list of predictions, one per vertex in map, of the {polar angle, eccentricity, visual area} for the vertex in the given SchiraModelObject model. For any vertex that lies outside of the model's bounds, the absolute value of the visual area will be greater than 4. A list of vertices or a single vertex may also be substituted for map.
CorticalMapToVisualField[model, X, Y] is equivalent to CorticalMapToVisualField[model, Thread[{X,Y}]].
CorticalMapToVisualField[model] yields a pure curried function that may be called with map directly.";
VisualFieldToCorticalMap::usage = "VisualFieldToCorticalMap[model, retinotopy] yields a list of 5 x 2 matrices, each row of which gives the {x, y} coordinate predictions for one of the visual areas. The result contains one such matrix for each retinotopic coordinate given. The retinotopy argument must be a list of {polarAngle, eccentricity}. The rows of each coordinate matrix returned represent, in order, the V1, V2, V3, HV4, and V3A predictions for the given retinotopic coordinate. A single retinotopy coordinate may be given, in which case, a single coordinate matrix is returned.
VisualFieldToCorticalMap[model, polarAngles, eccentricities] is equivalent to VisualFieldToCorticalMap[model, Thread[{polarAngles, eccentricitie}]].
VisualFieldToCorticalMap[model] yields a pure curried function that may be called with retinotopy directly.";

PolarAngleLegend::usage = "PolarAngleLegend[hemi] yields a graphic that is appropriate for a polar angle legend. All options that are valid for DensityPlot can be passed.";
EccentricityLegend::usage = "EccentricityLegend[hemi,max] yields a graphic that is appropriate for an eccentricity legend. All options that are valid for DensityPlot can be passed. The range is an optional argument that specifies the max eccentricity for this legend (default: 20)";

SchiraParametricPlot::usage = "SchiraParametricPlot[mdl, options...] is equivalent to ParametricPlot[f[th,r], {th, thMin, thMax}, {r, rMin, rMax}, options] where f is the Schira function for the given SchiraModelObject mdl, and thMin, thMax, rMin, and rMax can be controlled via the additional option Range. The visual areas specified by the VisualAreas argument are plotted. The Range argument may be of the following forms: {thRange, rRange} or rRange where thRange must be {thMin, thMax} and rRange may be rMax (with rMin 0) or {rMin, rMax}.";
SchiraParametricPlot::badarg = "Bad argument to SchiraParametricPlot: `1`";
V123ParametricPlot::usage = "V123ParametricPlot[mdl, options...] is equivalent to ParametricPlot[f[th,r], {th, thMin, thMax}, {r, rMin, rMax}, options] where f is the function for the given V123Model object mdl, and thMin, thMax, rMin, and rMax can be controlled via the additional option Range. The visual areas specified by the VisualAreas argument are plotted. The Range argument may be of the following forms: {thRange, rRange} or rRange where thRange must be {thMin, thMax} and rRange may be rMax (with rMin 0) or {rMin, rMax}.";
V123ParametricPlot::badarg = "Bad argument to V123ParametricPlot: `1`";
VisualAreas::usage = "VisualAreas is an optional keyword argument to the SchiraParametricPlot and V123ParametriclPlot functions. VisualAreas may be a list of any of {-4,-3,-2,-1,1,2,3,4}, which represent areas hV4, V3ventral, V2ventral, V1Ventral, V1Dorsal, V2Dorsal, V3Dorsal, V3A, respectively. The default value or Automatic will yield {-3,-2,-1,1,2,3}.";

SchiraLinePlot::usage = "SchiraLinePlot[mdl,options...] is equivalent to calling ParametricPlot[] with the given options, but such that the first arguments are optimized to draw polar angle and eccentricity lines from the given SchiraModelObject mdl. The following additional parameters may be given:
PolarAngleLines gives the polar angle values at which to daw iso-eccentric lines.
EccentricityLines gives the eccentricity values at which to draw iso-angular lines.
PolarAngleStyleFunction gives a function that, given a polar angle value, yields the style instruction for that polar angle line.
EccentricityFunction gives a function that, given a polar angle value, yields the style instruction for that polar angle line.
VisualAreas is used as in SchiraParametricPlot.
Range is used as in SchiraParametricPlot.";
SchiraLinePlot::badarg = "Bad argument to SchiraLinePlot: `1`";
V123LinePlot::usage = "V123LinePlot[mdl,options...] is equivalent to calling ParametricPlot[] with the given options, but such that the first arguments are optimized to draw polar angle and eccentricity lines from the given V123Model association mdl. The following additional parameters may be given:
PolarAngleLines gives the polar angle values at which to daw iso-eccentric lines.
EccentricityLines gives the eccentricity values at which to draw iso-angular lines.
PolarAngleStyleFunction gives a function that, given a polar angle value, yields the style instruction for that polar angle line.
EccentricityFunction gives a function that, given a polar angle value, yields the style instruction for that polar angle line.
VisualAreas is used as in V123ParametricPlot.
Range is used as in V123ParametricPlot.";
V123LinePlot::badarg = "Bad argument to V123LinePlot: `1`";
PolarAngleLines::usage = "PolarAngleLines is an option to SchiraLinePlot and V123LinePlot that specifies the polar angle values at which to draw the iso-eccentric lines.";
EccentricityLines::usage = "EccentricityLines is an option to SchiraLinePlot and V123LinePlot that specifies the eccentricity values at which to draw the iso-angular lines.";
PolarAngleStyleFunction::usage = "PolarAngleStyleFunction is an option to SchiraLinePlot and V123LinePlot that must accept a polar angle value and yield the style directive for plotting that angle's iso-eccentric line. The value Automatic results in ColorCortex[PolarAngle] beign used. The values Thick, Dashed, Dotted, etc. will result in the same coloring schele with the option directive appended.";
EccentricityStyleFunction::usage = "EccentricityStyleFunction is an option to SchiraLinePlot and V123LinePlot that must accept an eccentricity value and yield the style directive for plotting that eccentricity's iso-angular line. The value Automatic results in ColorCortex[Eccentricity] beign used. The values Thick, Dashed, Dotted, etc. will result in the same coloring schele with the option directive appended.";

ManipulateSchiraModel::usage = "ManipulateSchiraModel[map] yields a manipulate form that allows one to edit the parameters of the Schira model which is projected over a cortex plot of the given map.
ManipulateSchiraModel[map, model] uses the given model as the starting parameterization of the model.
ManipulateSchiraModel[map, model :> var] assigns the updated model value to the given var every time a change is entered, assuming that the \"Save on change?\" option is set to True.";

(* Registration potential functions ***************************************************************)
GaussianSchiraPotential::usage = "GaussianSchiraPotential[map, model] yields a potential function that describes the agreement between the given SchiraModel object, model, and the retinotopy data found in the properties \"PolarAngle\" and \"Eccentricity\" of map. If the property \"VertexWeight\" is also present, then all vertices with a value above 0 are included and they are weighted by the given weights. Otherwise, any vertex with missing polar angle or eccentricity values ($Failed, None, or Indeterminate) are excluded.";
GaussianSchiraPotential::badarg = "Bad argument given to GaussianSchiraPotential: `1`";

SchiraAnchors::usage = "SchiraAnchors[map, model] yields a list appropriate for specification of a Schira potential of the given model over the given mesh via the \"Anchors\" potential type of the PotentialField function. Generally, this would be used as PotentialField[mesh, \"Anchors\" -> SchiraAnchors[mesh, model]]. To construct these anchors, the function uses the retinotopy data found in the mesh properties \"PolarAngle\" and \"Eccentricity\" of the given mesh. If the property \"VertexWeight\" or \"Weight\" is also present, then all vertices with a value above 0 are included and they are weighted by the given weights. Otherwise, any vertex with missing polar angle or eccentricity values ($Failed, None, or Indeterminate) are excluded.";
SchiraAnchors::badarg = "Bad argument given to SchiraAnchors: `1`";
SchiraAnchorsParameter::usage = "SchiraAnchorsParameter[map, param] yields a list appropriate for specification of a Schira potential parameters via the \"Anchors\" potential type of the PotentialField function. Generally, this would be used as PotentialField[mesh, \"Anchors\" -> SchiraAnchors[mesh, model], Scale -> SchiraAnchorsParameter[map, scaleData]].";

V123Anchors::usage = "V123Anchors[map, model] yields a list appropriate for specification of a V1-V3 model potential of the given model over the given mesh via the \"Anchors\" potential type of the PotentialField function. Generally, this would be used as PotentialField[mesh, \"Anchors\" -> V123Anchors[mesh, model]]. To construct these anchors, the function uses the retinotopy data found in the mesh properties \"PolarAngle\" and \"Eccentricity\" of the given mesh. If the property \"VertexWeight\" or \"Weight\" is also present, then all vertices with a value above 0 are included and they are weighted by the given weights. Otherwise, any vertex with missing polar angle or eccentricity values ($Failed, None, or Indeterminate) are excluded.";
V123Anchors::badarg = "Bad argument given to V123Anchors: `1`";
V123AnchorsParameter::usage = "V123AnchorsParameter[map, param] yields a list appropriate for specification of a V1-V3 model potential parameters via the \"Anchors\" potential type of the PotentialField function. Generally, this would be used as PotentialField[mesh, \"Anchors\" -> V123Anchors[mesh, model], Scale -> V123AnchorsParameter[map, scaleData]].";

MeshRegionOrthogonalFields::usage = "MeshRegionOrthogonalFields[mesh, init1, init2] yields a pair of orthogonal fields at the vertex coordinates in the given mesh region; the fields are discovered by minimization in which the smoothness and the orthogonality of the fields are minimized such that the given init values are held constant. The init values should be lists of (vertexIndex -> value) rules. All options that can be given the FindArgMin may be passed to this function.";

V123Model::usage = "V123Model[] yields a data structure (an Association) of relevant data for a model of V1-V3 on the flattened cortical surface. The data structure includes the following elements:
    * \"Mesh\" holds the MeshRegion of the model;
    * \"BoundaryMesh\" holds a BoundaryMeshRegion of the model;
    * \"PolarAngle\" and \"Eccentricity\" hold the data at the mesh points;
    * \"PolarAngleConstraints\" and \"EccentricityConstraints\" hold the data used to initialize the mesh;
    * \"V1BoundaryMesh\", \"V2BoundaryMesh\", \"V3BoundaryMesh\", and \"V4BoundaryMesh\" hold boundary mesh regions for the individual visual areas.
  
  The following options may be given:
    * Scale (default: {15,5,5,3}) must be a 4-element list of positive numbers corresponding to the widths of each visual area (V1, V2, V3, V4);
    * AspectRatio (default: 0.32) specifies the height/width ratio of the visual areas;
    * Intersection (default: Pi * 3/5) specifies the angle at which the dorsal and ventral arms of the model intersect on the X-axis;
    * Exponent (default: 3.5) specifies the eccentricity magnification m used in the eccentricity function (M * x^m);
    * MaxValue (default: 90.0) specifies the maximum eccentricity value, M;
    * CellSize (default: 0.1) specifies the approximate spacing of coordinates in the mesh;
    * MaxIterations (default: 10,000) specifies the maximum number of iterations when finding the polar angle and eccentricity values;
    * AffineTransform (default: a shear of {{1,-0.2},{0,1}} followed by a 5\[Degree] rotation and a translation of {-7,-1}) specifies the transform to be applied to the mesh after it has been solved.";
V123ModelQ::usage = "V123ModelQ[mdl] yields true if the given obejct, mdl, is a valid V123 model, otherwise False.";

(**************************************************************************************************)
(**************************************************************************************************)
Begin["`Private`"];

(* Visual Angles and Visual Areas *****************************************************************)
VisualAngleToComplex[th_, r_] := r * Exp[I * (Pi/2 - th * Pi / 180)];
ComplexToVisualAngle[z_] := {180/Pi * (Pi/2 - Arg[z]), Abs[z]};
ComplexToPolarAngle[z_] := 180/Pi * (Pi/2 - Arg[z]);
ComplexToEccentricity[z_] := Abs[z];
ComplexToCoordinate[z_] := {Re[z],Im[z]};
SetAttributes[ComplexToCoordinate, Listable];
CoordinateToComplex[x_,y_] := x + I*y;
CoordinateToComplex[xy_List] := Which[
  Length[xy] == 0, {},
  ListQ[First[xy]], With[{tr = Transpose[xy]}, tr[[1]] + I*tr[[2]]],
  True, xy[[1]] + I * xy[[2]]];
CoordinateToVisualAngle[x_, y_] := {90 - 180/Pi*ArcTan[x,y], Sqrt[x^2 + y^2]};
CoordinateToVisualAngle[{x:Except[_List], y:Except[_List]}] := {
  90 - 180/Pi*ArcTan[x,y],
  Sqrt[x^2 + y^2]};
CoordinateToVisualAngle[l_List /; MatchQ[Dimensions[l], {_,2}]] := With[
  {tr = Transpose[l]},
  Transpose[
    {90 - 180/Pi*ArcTan[tr[[1]], tr[[2]]],
     Sqrt[tr[[1]]^2 + tr[[2]]^2]}]];
VisualAngleToCoordinate[t_, r_] := {r * Cos[Pi/180*(90 - t)], r * Sin[Pi/180*(90 - t)]};
VisualAngleToCoordinate[{t:Except[_List], r:Except[_List]}] := r * {
  Cos[Pi/180*(90 - t)],
  Sin[Pi/180*(90 - t)]};
VisualAngleToCoordinate[l_List /; MatchQ[Dimensions[l], {_,2}]] := With[
  {tr = Transpose[l]},
  Transpose[
    tr[[2]] * {Cos[Pi/180*(90 - tr[[1]])], Sin[Pi/180*(90 - tr[[1]])]}]];
Protect[VisualAngleToComplex, ComplexToVisualAngle, ComplexToPolarAngle, ComplexToEccentricity,
        CoordinateToComplex, ComplexToCoordinate, CoordinateToVisualAngle, VisualAngleToCoordinate];

$VisualAreasData = {
  1  -> {"Name" -> "V1 Dorsal",  "Areas" -> {"V1"}, "Simple" -> 1, "Stream" -> "Dorsal"},
  -1 -> {"Name" -> "V1 Ventral", "Areas" -> {"V1"}, "Simple" -> 1, "Stream" -> "Ventral"},
  2  -> {"Name" -> "V2 Dorsal",  "Areas" -> {"V2"}, "Simple" -> 2, "Stream" -> "Dorsal"},
  -2 -> {"Name" -> "V2 Ventral", "Areas" -> {"V2"}, "Simple" -> 2, "Stream" -> "Ventral"},
  3  -> {"Name" -> "V3 Dorsal",  "Areas" -> {"V3"}, "Simple" -> 3, "Stream" -> "Dorsal"},
  -3 -> {"Name" -> "V3 Ventral", "Areas" -> {"V3"}, "Simple" -> 3, "Stream" -> "Ventral"},
  4  -> {"Name" -> "V3A", "Areas" -> {"V3A"}, "Simple" -> 0, "Stream" -> "Dorsal"},
  -4 -> {"Name" -> "HV4", "Areas" -> {"HV4"}, "Simple" -> 0, "Stream" -> "Ventral"},
  5  -> {"Name" -> "Unknown Dorsal Area",  "Areas" -> {}, "Simple" -> 0, "Stream" -> "Dorsal"},
  -5 -> {"Name" -> "Unknown Ventral Area", "Areas" -> {}, "Simple" -> 0, "Stream" -> "Ventral"},
  0  -> {"Name" -> "V1 Horizontal Meridian", "Areas" -> {"V1"}, "Simple" -> 1, "Stream" -> None},
  (3/2)  -> {"Name" -> "V1/V2 Lower Vertical Meridian", "Areas" -> {"V1", "V2"}, "Simple" -> 1, "Stream" -> "Dorsal"},
  (-3/2) -> {"Name" -> "V1/V2 Upper Vertical Meridian", "Areas" -> {"V1", "V2"}, "Simple" -> 1, "Stream" -> "Ventral"},
  (5/2)  -> {"Name" -> "V2/V3 Lower Vertical Meridian", "Areas" -> {"V2", "V3"}, "Simple" -> 2, "Stream" -> "Dorsal"},
  (-5/2) -> {"Name" -> "V2/V3 Upper Vertical Meridian", "Areas" -> {"V2", "V3"}, "Simple" -> 2, "Stream" -> "Ventral"},
  (7/2)  -> {"Name" -> "V3/V3A Lower Vertical Meridian", "Areas" -> {"V3", "V3A"}, "Simple" -> 3, "Stream" -> "Dorsal"},
  (-7/2) -> {"Name" -> "V3/HV4 Upper Vertical Meridian", "Areas" -> {"V3", "HV4"}, "Simple" -> 3, "Stream" -> "Ventral"},
  (9/2)  -> {"Name" -> "V3A Upper Vertical Meridian", "Areas" -> {"V3A"}, "Simple" -> 0, "Stream" -> "Dorsal"},
  (-9/2) -> {"Name" -> "HV4 Lower Vertical Meridian", "Areas" -> {"HV4"}, "Simple" -> 0, "Stream" -> "Ventral"}};
$VisualAreasDispatch = Dispatch[
  Append[
    $VisualAreasData,
    a_ :> (Message[VisualAreaData::badarea, a]; Indeterminate)]];

VisualAreaData[] := $VisualAreasData;
VisualAreaData[id_] := Replace[id, $VisualAreasDispatch];
VisualAreaName[id_] := Check[Replace["Name", Replace[id, $VisualAreasDispatch]], Indeterminate];
VisualAreaSimplify[id_] := Check[Replace["Simple", Replace[id, $VisualAreasDispatch]], Indeterminate];
SetAttributes[VisualAreaData, Listable];
SetAttributes[VisualAreaName, Listable];
SetAttributes[VisualAreaSimplify, Listable];
Protect[$VisualAreasData, $VisualAreasDispatch, VisualAreaData, VisualAreaName, VisualAreaSimplify];

(* Retinotopic Templates **************************************************************************)
$FSAverageSymRetinotopicTemplateAddresses = Association[
  {Rule[
     PolarAngle,
     "https://cfn.upenn.edu/aguirreg/public/ES_template/mgh_files/angle-template-2.5.sym.mgh"],
   Rule[
     Eccentricity,
     "https://cfn.upenn.edu/aguirreg/public/ES_template/mgh_files/eccen-template-2.5.sym.mgh"],
   Rule[
     VisualArea,
     "https://cfn.upenn.edu/aguirreg/public/ES_template/mgh_files/areas-template-2.5.sym.mgh"]}];
Protect[$FSAverageSymRetinotopicTemplateAddresses];
FSAverageSymRetinotopy := With[
  {field = Check[
     Association @ MapThread[
       (* #here -- why does this need to be reversed?! *)
       Rule[#1, Reverse[#2]]&,
       {{PolarAngle, Eccentricity, VisualArea},
        Transpose @ MapThread[
          Function[{pa, e, a},
            Which[
              Abs[a] < 1 || Abs[a] > 3,  {None, None, None},
              NumberQ[pa] && pa > 90.0,  {pa, e, -Abs[a]},
              NumberQ[pa] && pa <= 90.0, {pa, e, Abs[a]},
              True, {None, None, None}]],
          Import[$FSAverageSymRetinotopicTemplateAddresses[#], "MGH"]& /@ {
            PolarAngle, Eccentricity, VisualArea}]}],
     $Failed]},
  If[field === $Failed,
    $Failed,
    (Unprotect[FSAverageSymRetinotopy];
     Set[FSAverageSymRetinotopy, field];
     Protect[FSAverageSymRetinotopy];
     FSAverageSymRetinotopy)]];
FSAverageSymPolarAngle := FSAverageSymRetinotopy[PolarAngle];
FSAverageSymEccentricity := FSAverageSymRetinotopy[Eccentricity];
FSAverageSymVisualArea := FSAverageSymRetinotopy[VisualArea];

PredictRetinotopy[sub_, meshName_, hemi_] := Check[
  If[sub === FSAverageSymSubject,
    SetProperty[
      {Cortex[sub, meshName, hemi], VertexList},
      Map[ToString[#] -> FSAverageSymRetinotopy[#]&, Keys[FSAverageSymRetinotopy]]],
    With[
      {sym = Cortex[sub, "Sym", hemi],
       props = {"PolarAngle", "Eccentricity", "VisualArea"},
       fsaverageSym = SetProperty[
         Cortex[FSAverageSymSubject, "Sphere", LH],
         Map[ToString[#] -> FSAverageSymRetinotopy[#]&, {PolarAngle, Eccentricity, VisualArea}]]},
      With[
        {resamp = CortexResample[
           fsaverageSym -> sym,
           Properties -> props]},
        Fold[
          SetProperty[{#1, VertexList}, #2 -> PropertyValue[{resamp, VertexList}, #2]]&,
          Cortex[sub, meshName, hemi],
          props]]]],
  $Failed];

Protect[
  FSAverageSymRetinotopy, FSAverageSymPolarAngle, FSAverageSymEccentricity, FSAverageSymVisualArea,
  PredictRetinotopy];


(* Default Schira Parameters **********************************************************************)
$DefaultSchiraA = 1.5;
$DefaultSchiraB = 45.0;
$DefaultSchira\[CapitalLambda] = 3.0;
$DefaultSchira\[CapitalPsi] = 0.70;
$DefaultSchiraV1Size = 1.10;
$DefaultSchiraV2Size = 0.45;
$DefaultSchiraV3Size = 0.35;
$DefaultSchiraHV4Size = 0.3;
$DefaultSchiraV3ASize = 0.3;
$DefaultSchiraFC = {-0.15, -0.4};
$DefaultSchiraScale = {40.0, 40.0};
$DefaultSchiraShear = {{1, 0}, {-0.2, 1}};
$DefaultSchira\[CapitalRho]90 = 90;
$SchiraParameters = List[
   A :> $DefaultSchiraA,
   B :> $DefaultSchiraB,
   \[CapitalLambda] :> $DefaultSchira\[CapitalLambda],
   \[CapitalPsi] :> $DefaultSchira\[CapitalPsi],
   V1Size :> $DefaultSchiraV1Size,
   V2Size :> $DefaultSchiraV2Size,
   V3Size :> $DefaultSchiraV3Size,
   HV4Size :> $DefaultSchiraHV4Size,
   V3ASize :> $DefaultSchiraV3ASize,
   FC :> $DefaultSchiraFC,
   Scale :> $DefaultSchiraScale,
   Shear :> $DefaultSchiraShear,
   \[CapitalRho]90 :> $DefaultSchira\[CapitalRho]90];

Protect[ 
  $DefaultSchiraA, $DefaultSchiraB,
  $DefaultSchira\[CapitalLambda], $DefaultSchira\[CapitalPsi],
  $DefaultSchira\[CapitalRho]90, $DefaultSchiraV1Size,
  $DefaultSchiraV2Size, $DefaultSchiraV3Size, $DefaultSchiraHV4Size,
  $DefaultSchiraV3ASize, $DefaultSchiraFC, $DefaultSchiraScale,
  $DefaultSchiraShear, $SchiraParameters];

(* Schira Parameters Keys *************************************************************************)
A = A; B = B; \[CapitalLambda] = \[CapitalLambda];
V1Size = V1Size; V2Size = V2Size; V3Size = V3Size; HV4Size = HV4Size; V3ASize = V3ASize; 
\[CapitalPsi] = \[CapitalPsi]; FC = FC; Shear = Shear;
\[CapitalRho]90 = \[CapitalRho]90;
Protect[A, B, \[CapitalLambda], 
       V1Size, V2Size, V3Size, V3ASize, HV4Size,
       \[CapitalPsi], \[CapitalRho]90, FC, Shear];

(* The Schira Model Objects ***********************************************************************)

(* private variable used below *)
Unprotect[$SchiraParameterPositions];
ClearAll[$SchiraParameterPositions];
$SchiraParameterPositions = Dispatch[
  MapThread[
    Rule,
    {$SchiraParameters[[All, 1]], Range[Length[$SchiraParameters]]}]];
Protect[$SchiraParameterPositions];

(* These actually compile the low-level Schira calculation functions *)
Unprotect[CompileSchiraFunction, CompileSchiraInverse];
ClearAll[CompileSchiraFunction, CompileSchiraInverse];
CompileSchiraFunction[a_, b_, lambda_, psi_, shearMtx_, scale_, fc_, areas_] := Check[
  With[
    {v1b = N[areas[[1]]],
     v2b = N[areas[[2]]],
     v3b = N[areas[[3]]],
     hv4b = N[areas[[4]]],
     v3ab = N[areas[[5]]],
     dsech1 = 0.1821,
     dsech2 = 0.76,
     fcx0 = N[Log[(a + lambda) / (b + lambda)]],
     xscale = If[NumberQ[scale], scale, scale[[1]]],
     yscale = If[NumberQ[scale], scale, scale[[2]]],
     mtx = Dot[
       RotationMatrix[N[psi]],
       N[shearMtx]]},
    (* sanity checks should already be done: just compile the function *)
    Compile[
      {{z, _Complex}},
      (* Layerless Transorm:
         First, set zLayered to the layerless tranform for each of the 4 areas *)
      With[
        {zLayered = With[
          {zz = Arg[z]*2.0/Pi},
          With[
            {sgn = If[zz == 0, 1.0, Sign[zz]]},
            Abs[z] * Exp[I * {
              v1b*zz,
              sgn * (v1b + v2b*(1.0 - Abs[zz])),
              sgn * (v1b + v2b + v3b*Abs[zz]),
              v1b + v2b + v3b + hv4b * (1.0 - 0.5 * (zz + 1.0)),
              -(v1b + v2b + v3b + 0.5 * v3ab * (zz + 1.0))}]]]},
        (* Log-Polar Transform:
           Now, do the log-polar part of the transform *)
        With[
          {zLogPolar = Table[
             With[
               {argz = Arg[
                  If[Re[zz] >= 0, zz + lambda, zz + 2.0 * lambda * (1.0 - Abs[Arg[zz]]/Pi)]],
                absz = Abs[
                  If[Re[zz] >= 0, zz + lambda, zz + 2.0 * lambda * (1.0 - Abs[Arg[zz]]/Pi)]]},
               If[absz == 0,
                 Log[a/b] + 0.0 * I,
                 Log@Divide[
                   a + absz * Exp[I * argz * Sech[argz]^(dsech1 * Sech[dsech2*Log[absz/a]])],
                   b + absz * Exp[I * argz * Sech[argz]^(dsech1 * Sech[dsech2*Log[absz/b]])]]]],
             {zz, zLayered}]},
          (* We center the FC on zero to start, by subtracting fcx0, then we scale and shear, and,
             last, we push things back to the specified FC *)
          fc[[1]] + I*fc[[2]] + Flatten@Dot[
            {{xscale, I*yscale}},
            mtx,
            (* Note that we flip the z here so that the arrangement matchis the LH *)
            {Re[zLogPolar] - fcx0, -Im[zLogPolar]}]]],
      RuntimeOptions -> {"Speed", "EvaluateSymbolically" -> True},
      Parallelization -> True,
      RuntimeAttributes -> {Listable}]],
  $Failed];

CompileSchiraInverse[a_, b_, lambda_, psi_, shearMtx_, scale_, fc_, areas_] := Check[
  With[
    {v1b = N[areas[[1]]],
     v2b = N[areas[[2]]],
     v3b = N[areas[[3]]],
     hv4b = N[areas[[4]]],
     v3ab = N[areas[[5]]],
     maxAngle = Total[N[areas[[{1,2,3,5}]]]],
     minAngle = Total[N[areas[[{1,2,3,4}]]]],
     dsech1 = 0.1821,
     dsech2 = 0.76,
     tol = 0.000001,
     fcx0 = N[Log[(a + lambda) / (b + lambda)]],
     xscale = If[NumberQ[scale], scale, scale[[1]]],
     yscale = If[NumberQ[scale], scale, scale[[2]]],
     mtx = Dot[
       RotationMatrix[N[psi]],
       N[shearMtx]](*,
     sideFn = With[
       {u = If[ListQ[scale] && Times@@scale < 0,
          {Cos[psi + 0.5*Pi], Sin[psi + 0.5*Pi]},
          {Cos[psi - 0.5*Pi], Sin[psi - 0.5*Pi]}]},
       Function[Sign[Dot[{Re@#, Im@#} - fc, u]]]]*)},
    (* sanity checks should already be done: we compile the function next:
       this function is basically identical to the compiled function just above, but it does not
       perform the layered transform (ie, translating area V1-V4 into an imaginary number. This
       allows us to make a simple inverse function *)
    With[
      {forward = Compile[
         {{z, _Complex}},
         (* Log-Polar Transform: Do the log-polar part of the transform *)
         With[
           {ztr = With[
              {argz = Arg[
                 If[Re[z] >= 0, z + lambda, z + 2.0 * lambda * (1.0 - Abs[Arg[z]]/Pi)]],
               absz = Abs[
                 If[Re[z] >= 0, z + lambda, z + 2.0 * lambda * (1.0 - Abs[Arg[z]]/Pi)]]},
              If[absz == 0,
                Log[a/b] + 0.0 * I,
                Log@Divide[
                  a + absz * Exp[I * argz * Sech[argz]^(dsech1 * Sech[dsech2*Log[absz/a]])],
                  b + absz * Exp[I * argz * Sech[argz]^(dsech1 * Sech[dsech2*Log[absz/b]])]]]]},
           (* We center the FC on zero to start, by subtracting fcx0, then we scale and shear, 
              and, last, we push things back to the specified FC *)
           fc[[1]] + I*fc[[2]] + First@Dot[
             {{xscale, I*yscale}},
             mtx,
             {Re[ztr] - fcx0, Im[ztr]}]],
         RuntimeOptions -> {"Speed", "EvaluateSymbolically" -> False},
         Parallelization -> True,
         RuntimeAttributes -> {Listable}]},
      With[
        {samples = Nearest[
           Apply[
             Rule,
             Transpose@Flatten[
               Table[
                 With[
                   {z = (90.0*r^3.5)*Exp[I*t]},
                   {{Re@#, Im@#}, z} & @ forward[z]],
                 {t, -minAngle, maxAngle, 0.01 * (maxAngle + minAngle)},
                 {r, 0, 1, 0.01}],
               1]]],
         distTol = EuclideanDistance[
           {Re@#, Im@#}&@forward[10.0],
           {Re@#, Im@#}&@forward[15.0]],
         distTolSmall = EuclideanDistance[
           {Re@#, Im@#}&@forward[10.0],
           {Re@#, Im@#}&@forward[11.0]]},
        (* Now, we create an inverse function for the forward function *)
        With[
          {inverse = Function[
             With[
               {close = samples[{Re@#, Im@#}, {1, distTol}]},
               If[Length[close] == 0,
                 $Failed,
                 Conjugate@FindRoot[
                   # == forward[w],
                   {w, close[[1]] - 0.001 - 0.001*I, close[[1]] + 0.001 + 0.001*I}
                  ][[1,2]]]]]},
          (* And wrap this inverse in a translator for the areas *)
          Function[{z},
            With[
              {w = inverse[z]},
              If[
                !And[
                  NumberQ[w],
                  Norm[{Re@z, Im@z} - {Re@#, Im@#}]&@forward[Conjugate@w] <= distTolSmall],
                {0.0, 0},
                With[
                  {argw = Arg[w],
                   absw = Abs[w],
                   side = Sign[Arg[w]]},
                  Which[
                    Round[Abs@argw, tol] == 0, {absw + 0.0*I, 0},
                    Round[v1b - Abs@argw, tol] >= 0, {
                      absw*Exp[I * Pi/2 * argw / v1b], 
                      If[side == 0, 1, side] * If[Round[Abs@argw - v1b, tol] == 0, 3/2, 1]},
                    Round[v2b - (Abs@argw - v1b), tol] >= 0, {
                      absw*Exp[I * Pi/2 * Sign[argw] * (1.0 - (Abs[argw] - v1b)/v2b)], 
                      side * If[Round[Abs@argw - v1b - v2b, tol] == 0, 5/2, 2]},
                    Round[v3b - (Abs@argw - v1b - v2b), tol] >= 0, {
                      absw*Exp[I * Pi/2 * Sign[argw] * (Abs[argw] - v1b - v2b)/v3b],
                      side * If[Round[Abs@argw - v1b - v2b - v3b, tol] == 0, 7/2, 3]},
                    side < 0 && Round[hv4b - (Abs@argw - v1b - v2b - v3b), tol] >= 0, {
                      absw*Exp[I * Pi * ((Abs@argw - v1b - v2b - v3b)/hv4b - 0.5)],
                      If[Round[Abs@argw - v1b - v2b - v3b - hv4b, tol] == 0, -9/2, -4]},
                    side > 0 && Round[v3ab - (Abs@argw - v1b - v2b - v3b), tol] >= 0, {
                      absw*Exp[I * Pi * (0.5 - (Abs@argw - v1b - v2b - v3b)/v3ab)],
                      If[Round[Abs@argw - v1b - v2b - v3b - v3ab, tol] == 0, 9/2, 4]},
                    side < 0, {
                      absw*Exp[I * Pi * (
                        0.5 - (Abs@argw - v1b - v2b - v3b - hv4b)/(Pi - v1b - v2b - v3b - hv4b))],
                      -5},
                    side > 0, {
                      absw*Exp[I * Pi * (
                        (Abs@argw - v1b - v2b - v3b - v3ab)/(Pi - v1b - v2b - v3b - v3ab) - 0.5)],
                      5},
                    True, {-absw, Infinity}]]]],
            {Listable}]]]]],
  $Failed];
Protect[CompileSchiraFunction, CompileSchiraInverse];


(* This is used to force preparation of the function and inverse below *)
Unprotect[SchiraModelObjectPrep];
ClearAll[SchiraModelObjectPrep];
SchiraModelObjectPrep[params_List] := With[
  {ff = Unique["fun"],
   if = Unique["inv"],
   prf = Unique["predictRetinotopy"],
   prc = Unique["predictCoordinates"],
   a = A /. params,
   b = B /. params,
   lambda = \[CapitalLambda] /. params,
   areas = {V1Size, V2Size, V3Size, HV4Size, V3ASize} /. params,
   (* These we may change so have to update in params *)
   psi = (\[CapitalPsi] /. params) /. None -> 0,
   shearMtx = (Shear /. params) /. None -> {{1,0},{0,1}},
   scale = Replace[Scale /. params, x_?NumericQ :> {x,x}],
   fc = Replace[FC /. params, z_Complex :> {Re[x], Im[z]}]},
  (* sanity checking *)
  Which[
    !NumericQ[a] || a <= 0, Message[SchiraModelObject::badarg, "A must be numeric and > 0"],
    !NumericQ[b] || b <= 0, Message[SchiraModelObject::badarg, "B must be numeric and > 0"],
    !NumericQ[lambda] || lambda < 0, Message[
      SchiraModelObject::badarg,
      "\[CapitalLambda] must be numeric and >= 0"],
    !NumericQ[psi], Message[
      SchiraModelObject::badarg,
      "\[CapitalPsi] must be numeric and >= 0"],
    Dimensions[shearMtx] != {2,2} || shearMtx[[1,1]] != 1 || shearMtx[[2,2]] != 1, Message[
      SchiraModelObject::badarg,
      "Shear must be a 2 x 2 matrix with ones on the diagonal or None"],
    !NumericQ[shearMtx[[1,2]]] || !NumericQ[shearMtx[[2,1]]], Message[
      SchiraModelObject::badarg,
      "Shear matrix must have numeric off-diagonal elements"],
    !ListQ[scale] || Length[scale] != 2 || !NumericQ[scale[[1]]] || !NumericQ[scale[[2]]], Message[
      SchiraModelObject::badarg,
      "Scale must be a single numeric quantity or a pair of {x-scale, y-scale} numeric quantities"],
    !ListQ[fc] || Length[scale] != 2 || !NumericQ[fc[[1]]] || !NumericQ[fc[[2]]], Message[
      SchiraModelObject::badarg,
      "Foveal convluence (FC) must be a either a coordinate {x, y} with numeric quantities or a"
       <> " single complex number with numeric real and imaginary parts"],
    Not[And@@Map[NumericQ, areas]] || Not[And@@Map[(#>0)&, areas]], Message[
      SchiraModelObject::badarg,
      "V1Size, V2Size, V3Size, HV4Size, and V3ASize must all be numeric quantities > 0"],
    Total[Most[areas]] > Pi || Total[areas[[{1,2,3,5}]]] > Pi, Message[
      SchiraModelObject::badarg,
      "V1Size, V2Size, V3Size, and neither HV4Size nor V3ASize may sum to be >= Pi"]];
  (* Note that these are temporary variables *)
  SetAttributes[Evaluate[ff], Temporary];
  SetAttributes[Evaluate[if], Temporary];
  SetAttributes[Evaluate[prf], Temporary];
  SetAttributes[Evaluate[prc], Temporary];
  (* set these to auto-memoize themselves if requested *)
  ff := With[
    {f = Check[
       CompileSchiraFunction[a, b, lambda, psi, shearMtx, scale, fc, areas],
       $Failed]},
    If[f === $Failed, $Failed, (ff = f)]];
  if := With[
    {f = Check[
       CompileSchiraInverse[a, b, lambda, psi, shearMtx, scale, fc, areas],
       $Failed]},
    If[f === $Failed, $Failed, (if = f)]];
  (* And make a dispatch for this model *)
  With[
    {mdl = SchiraModelObject[
      Dispatch[
        Join[
          Cases[
            params,
            Except[
              Alternatives[
                Rule[(\[CapitalPsi])|Shear|Scale|FC, _],
                RuleDelayed[(\[CapitalPsi])|Shear|Scale|FC, _]]]],
          {\[CapitalPsi] -> psi,
           Shear -> shearMtx,
           Scale -> scale,
           FC -> fc,
           Function :> ff,
           Inverse :> if,
           CorticalMapToVisualField :> prf,
           VisualFieldToCorticalMap :> prc,
           All :> params,
           x_ :> Message[SchiraModelObject::badarg, x]}]]]},
    prf := With[
      {res = Check[CorticalMapToVisualField[mdl], $Failed]},
      If[res === $Failed, res, (prf = res)]];
    prc := With[
      {res = Check[VisualFieldToCorticalMap[mdl], $Failed]},
      If[res === $Failed, res, (prc = res)]];
    mdl]];
Protect[SchiraModelObjectPrep];

SchiraModel[
   opts : Evaluate[
     Repeated[
      (Rule | RuleDelayed)[
       Apply[Alternatives, First /@ $SchiraParameters],
       _]]]
 ] := SchiraModelObjectPrep[
  ReplaceAll[
    $SchiraParameters,
    Map[
      Function[{rule},
        Rule[
          (First[rule] :> _),
          Head[rule][First[rule], Last[rule]]]],
        {opts}]]];
SchiraModel[
   SchiraModelObject[disp_],
   opts : Evaluate[
     Repeated[
      (Rule | RuleDelayed)[
       Apply[Alternatives, First /@ $SchiraParameters],
       _]]]
 ] := SchiraModelObjectPrep[
  ReplaceAll[
    Select[
      disp[[1]],
      And[Head[#[[1]]] =!= Pattern, MemberQ[$SchiraParameters[[All,1]], #[[1]]]]&],
    Map[
      Function[{rule},
        Rule[
          (Rule|RuleDelayed)[First[rule], _],
          Head[rule][First[rule], Last[rule]]]],
        {opts}]]];
SchiraModel[] := SchiraModelObjectPrep[$SchiraParameters];
 
SchiraModelObject[disp_][x_] := Replace[x, disp];
SchiraFunction[SchiraModelObject[disp_]] := Replace[Function, disp];
SchiraInverse[SchiraModelObject[disp_]] := Replace[Inverse, disp];

(* #CorticalMapToVisualField ***********************************************************************)
CorticalMapToVisualField[SchiraModelObject[disp_], map_?CorticalMapQ] := With[
  {inv = Replace[Inverse, disp],
   Z = VertexCoordinatesTr[map],
   r90 = Replace[\[CapitalRho]90, disp]},
  Quiet[
    Map[
      Append[ComplexToVisualAngle[#[[1]]], #[[2]]]&,
      inv[(Z[[1]] + I * Z[[2]])] * (90.0 / r90)],
    {FindRoot::cvmit, FindRoot::frmp, FindRoot::lstol}]];
CorticalMapToVisualField[SchiraModelObject[disp_], {x:Except[_List], y:Except[_List]}] := With[
  {inv = Replace[Inverse, disp],
   r90 = Replace[\[CapitalRho]90, disp]},
  With[
   {z = Quiet[inv[x + I*y] * (90.0 / r90), {FindRoot::cvmit, FindRoot::frmp, FindRoot::lstol}]},
   Append[ComplexToVisualAngle[z[[1]]], z[[2]]]]];
CorticalMapToVisualField[SchiraModelObject[disp_],
                        {x_List, y_List} /; Length[x] == Length[y] && Length[x] != 2] := With[
  {inv = Replace[Inverse, disp],
   r90 = Replace[\[CapitalRho]90, disp]},
  With[
    {res = Quiet[inv[X + I*Y] * (90.0 / r90), {FindRoot::cvmit, FindRoot::frmp, FindRoot::lstol}]},
    If[ListQ[First@res],
      Append[ComplexToVisualAngle[#[[1]]], #[[2]]]& /@ res,
      Append[ComplexToVisualAngle[res[[1]]], res[[2]]]]]];
CorticalMapToVisualField[SchiraModelObject[disp_], coords:{{_,_}..}] := With[
  {inv = Replace[Inverse, disp],
   Z = Transpose[coords],
   r90 = Replace[\[CapitalRho]90, disp]},
  Map[
    Append[ComplexToVisualAngle[#[[1]]], #[[2]]]&,
    Quiet[inv[Z[[1]] + I * Z[[2]]] * (90.0 / r90), {FindRoot::cvmit, FindRoot::frmp, FindRoot::lstol}]]];
CorticalMapToVisualField[SchiraModelObject[disp_], X_, Y_] := With[
  {inv = Replace[Inverse, disp],
   r90 = Replace[\[CapitalRho]90, disp]},
  With[
    {res = Quiet[inv[X + I*Y] * (90.0 / r90), {FindRoot::cvmit, FindRoot::frmp, FindRoot::lstol}]},
    If[ListQ[First@res],
      Append[ComplexToVisualAngle[#[[1]]], #[[2]]]& /@ res,
      Append[ComplexToVisualAngle[res[[1]]], res[[2]]]]]];
CorticalMapToVisualField[mdl_SchiraModelObject] := Function[CorticalMapToVisualField[mdl, ##]];

VisualFieldToCorticalMap[SchiraModelObject[disp_], retinotopy:{{_,_}..}] := With[
  {fun = Replace[Function, disp],
   tr = Transpose[retinotopy],
   r90 = Replace[\[CapitalRho]90, disp]},
  ComplexToCoordinate[fun[r90 / 90.0 * VisualAngleToComplex[tr[[1]], tr[[2]]]]]];
VisualFieldToCorticalMap[
  SchiraModelObject[disp_],
  {polarAngle:Except[_List], eccentricity:Except[_List]}
 ] := With[
  {fun = Replace[Function, disp],
   r90 = Replace[\[CapitalRho]90, disp]},
  ComplexToCoordinate[fun[r90 / 90.0 * VisualAngleToComplex[polarAngle, eccentricity]]]];
VisualFieldToCorticalMap[SchiraModelObject[disp_], polarAngles_, eccentricities_] := With[
  {fun = Replace[Function, disp],
   r90 = Replace[\[CapitalRho]90, disp]},
  ComplexToCoordinate[fun[r90 / 90.0 * VisualAngleToComplex[polarAngles, eccentricities]]]];
VisualFieldToCorticalMap[mdl_SchiraModelObject] := Function[VisualFieldToCorticalMap[mdl, ##]];

(* We also have a generic handler for associations that have the right keys: 
 * This expects the CorticalMapToVisualField key to contain a function that translates a point
 * x + I*y to x+I*y in the visual field.
 *)
CorticalMapToVisualField[assc_?AssociationQ /; KeyExistsQ[assc, "CorticalMapToVisualField"],
                         coords_?MatrixQ] := With[
  {f = assc["CorticalMapToVisualField"],
   regf = If[KeyExistsQ[assc, "VisualAreaFunction"],
     assc["VisualAreaFunction"],
     ConstantArray[Indeterminate, Length[#]]&]},
  With[
    {zz = f[coords], areas = regf[coords]},
    Transpose[{90 - Arg[zz]*180/Pi, Abs[zz], areas}]]];
CorticalMapToVisualField[
  assc_?AssociationQ /; KeyExistsQ[assc, "CorticalMapToVisualField"],
  map_?CorticalMapQ
  ] := CorticalMapToVisualField[assc, VertexCoordinates[map]];
CorticalMapToVisualField[
  assc_?AssociationQ /; KeyExistsQ[assc, "CorticalMapToVisualField"],
  x0_?VectorQ
  ] := First@CorticalMapToVisualField[assc, {x0}];
CorticalMapToVisualField[
  assc_?AssociationQ /; KeyExistsQ[assc, "CorticalMapToVisualField"],
  x_?VectorQ, y_?VectorQ
  ] := CorticalMapToVisualField[assc, Transpose[{x, y}]];
CorticalMapToVisualField[
  assc_?AssociationQ /; KeyExistsQ[assc, "CorticalMapToVisualField"],
  x_?NumericQ, y_?NumericQ
  ] := First@CorticalMapToVisualField[assc, {{x, y}}];

VisualFieldToCorticalMap[assc_?AssociationQ /; KeyExistsQ[assc, "VisualFieldToCorticalMap"],
                         angles:{{_,_}..}] := With[
  {f = assc["VisualFieldToCorticalMap"],
   a = Transpose[angles]},
  f@VisuaAngleToComplex[a[[1]], a[[2]]]];
VisualFieldToCorticalMap[assc_?AssociationQ /; KeyExistsQ[assc, "VisualFieldToCorticalMap"],
                         ang_List, ecc_List] := With[
  {f = assc["VisualFieldToCorticalMap"]},
  f@VisualAngleToComplex[ang, ecc]];
VisualFieldToCorticalMap[assc_?AssociationQ /; KeyExistsQ[assc, "VisualFieldToCorticalMap"],
                         ang:Except[_List], ecc:Except[_List]] := With[
  {f = assc["VisualFieldToCorticalMap"]},
  First@f@VisualAngleToComplex[{ang}, {ecc}]];


Protect[SchiraModel, SchiraModelObject, SchiraFunction, SchiraInverse,
        CorticalMapToVisualField, VisualFieldToCorticalMap];

(* Plotting Data **********************************************************************************)
PolarAngleLegend[hemi : (LH|RH), opts___Rule] := DensityPlot[
  If[hemi === LH, ArcTan[x, y], ArcTan[90 - x, y]],
  {x, 0, 90},
  {y, -90, 90},
  opts,
  RegionFunction -> If[hemi === LH, (Norm[{#1, #2}] < 90 &), (Norm[{90 - #1, #2}] < 90 &)],
  ColorFunctionScaling -> False,
  ColorFunction -> Function[CorticalColorData["PolarAngle"][<|"PolarAngle" -> 90.0 - #*180.0/Pi|>]],
  Frame -> False,
  Axes -> False,
  BaseStyle -> Directive[10, FontFamily -> "Arial"],
  ImageSize -> 1.25*72,
  AspectRatio -> 2,
  Background -> White];

EccentricityLegend[hemi : (LH | RH), max_?NumericQ /; 0 < max <= 90, opts___Rule] := DensityPlot[
  If[hemi === LH, Norm[{x, y}], Norm[{max - x, y}]],
  {x, 0, max},
  {y, -max, max},
  opts,
  RegionFunction -> If[hemi === LH, (Norm[{#1, #2}] < max &), (Norm[{max - #1, #2}] < max &)],
  ColorFunctionScaling -> False,
  ColorFunction -> Function[CorticalColorData["Eccentricity"][<|"Eccentricity" -> #|>]],
  Frame -> False,
  Axes -> False,
  BaseStyle -> Directive[10, FontFamily -> "Arial"],
  ImageSize -> 1.25*72,
  AspectRatio -> 2,
  Background -> White];


(* #SchiraParametricPlot **************************************************************************)
Options[SchiraParametricPlot] = Join[
   Options[ParametricPlot],
   {VisualAreas -> Automatic,
    Range -> Full}];
SchiraParametricPlot[mdl_SchiraModelObject, opts:OptionsPattern[]] := Catch[
  With[
    {epsilon = 0.000001,
     plotRangeArg = OptionValue[PlotRange],
     colorFun = OptionValue[ColorFunction],
     colorFunSc = OptionValue[ColorFunctionScaling],
     areas = Union@Replace[
        OptionValue[VisualAreas],
        {All -> {-4, -3, -2, -1, 1, 2, 3, 4},
         Automatic -> {-3, -2, -1, 1, 2, 3},
         i_Integer /; -5 < i < 5 && i != 0 :> {i},
         l_List /; Length[l] ==  Count[l, i_Integer /; -5 < i < 5 && i != 0, {1}] :> Union[l],
         _ :> Message[
           SchiraParametricPlot::badarg,
           "VisualAreas must be All, one of +/- {1,2,3,4}, or a list of such integers"]}],
     f = mdl[VisualFieldToCorticalMap],
     range = Replace[
       OptionValue[Range],
       {(All | Full | Automatic) -> {{0, 180}, {0, 90}},
        r : {{_, _}, {_, _}} :> r,
        {t : {_, _}, r : Except[{_, _}]} :> {t, {0, r}},
        {t : {_, _}, (All | Full | Automatic)} :> {t, {0, 90}},
        {(All | Full | Automatic), r : {_, _}} :> {{0, 180}, r},
        {(All | Full | Automatic), 
          r : Except[{_, _}]} :> {{0, 180}, {0, r}},
        r_ :> {{0, 180}, {0, r}}}]},
    With[
      {msg = Which[
         ! NumericQ[range[[1, 1]]], "theta-min must be a number",
         ! NumericQ[range[[1, 2]]], "theta-max must be a number",
         ! NumericQ[range[[2, 1]]], "rho-min must be a number",
         ! NumericQ[range[[2, 2]]], "rho-max must be a number",
         ! (0 <= range[[1, 1]] < range[[1, 2]]), "theta-min must be in [0,theta-max)",
         range[[1, 2]] > 180, "theta-max must be in (theta-min,180]",
         ! (0 <= range[[2, 1]] < range[[2, 2]]), "rho-min must be in [0,rho-max)",
         range[[2, 2]] > 90, "rho-max must be in (rho-min,90]",
         True, None]},
      If[StringQ[msg], (Message[SchiraParametricPlot::badarg, msg]; Throw[$Failed])]];
    With[
     {optseq = Sequence @@ FilterRules[{opts}, Options[ParametricPlot][[All,1]]],
      optseqShow = Sequence @@ FilterRules[{opts}, Options[Show][[All,1]]],
      rhoTrans = Function[90.0*#^3.5],
      rhoTransInv = Function[(#/90.0)^(1/3.5)],
      thetaLower = If[range[[1, 1]] > 90,
        None,
        {range[[1, 1]], Min[{90, range[[1, 2]]}]}],
      thetaUpper = If[range[[1, 2]] < 90,
        None,
        {Max[{90, range[[1, 1]]}], range[[1, 2]]}]},
     With[
       {rMin = rhoTransInv[range[[2,1]]],
        rMax = rhoTransInv[range[[2,2]]]},
       With[
         {graphics = Map[
            Function[
              With[
                {k = Which[# == -4, 4, # == 4, 5, True, Abs[#]],
                 thetaMinIdeal = If[# == 4 || # < 0, thetaLower[[1]], thetaUpper[[1]]],
                 thetaMaxIdeal = If[# == -4 || # > 0, thetaUpper[[2]], thetaLower[[2]]]},
                With[
                  {thetaMin = If[# == 2 || # == 3, 
                     Max[{thetaMinIdeal, 90.0 + epsilon}],
                     thetaMinIdeal],
                   thetaMax = If[# == -3 || # == -2,
                     Min[{thetaMaxIdeal, 90.0 - epsilon}],
                     thetaMaxIdeal]},
                  Quiet[
                    ParametricPlot[
                      Part[f[theta, rhoTrans[rho]], k],
                      {theta, thetaMin, thetaMax},
                      {rho, rMin, rMax},
                      PlotRange -> plotRangeArg,
                      ColorFunction -> If[
                        Or[colorFun === Automatic,
                           colorFun === None,
                           colorFun === False,
                           StringQ[colorFun]],
                        colorFun,
                        If[colorFunSc === False,
                          Function[colorFun[#1,#2,#3,rhoTrans[#4]]],
                          Function[colorFun[#1,#2,#3,rhoTrans[#4]/90.0]]]],
                      optseq],
                    {CompiledFunction::cfsa}]]]],
            areas]},
         With[
           {plotRange = With[
              {ranges = Cases[
                 AbsoluteOptions[#,PlotRange]& /@ graphics,
                 (PlotRange -> r_) :> r,
                 {2}]},
              {{Min[ranges[[All, 1, 1]]], Max[ranges[[All, 1, 2]]]},
               {Min[ranges[[All, 2, 1]]], Max[ranges[[All, 2, 2]]]}}]},
           Show[
             graphics,
             PlotRange -> plotRange, 
             optseqShow]]]]]]];

(* #V123ParametricPlot ****************************************************************************)
Options[V123ParametricPlot] = Join[
   Options[ParametricPlot],
   {VisualAreas -> Automatic,
    Range -> Full,
    Hemi -> LH}];
V123ParametricPlot[mdl_?AssociationQ, opts:OptionsPattern[]] := Catch[
  With[
    {epsilon = 0.000001,
     hemi = OptionValue[Hemi] /. {Left -> LH, Right -> RH, Automatic -> LH},
     plotRangeArg = OptionValue[PlotRange],
     colorFun = OptionValue[ColorFunction],
     colorFunSc = OptionValue[ColorFunctionScaling],
     areas = Union@Replace[
        OptionValue[VisualAreas],
        {All -> {-4, -3, -2, -1, 1, 2, 3, 4},
         Automatic -> {-3, -2, -1, 1, 2, 3},
         i_Integer /; -5 < i < 5 && i != 0 :> {i},
         l_List /; Length[l] ==  Count[l, i_Integer /; -5 < i < 5 && i != 0, {1}] :> Union[l],
         _ :> Message[
           V123ParametricPlot::badarg,
           "VisualAreas must be All, one of +/- {1,2,3,4}, or a list of such integers"]}],
     f0 = Function@VisualFieldToCorticalMap[mdl, ##],
     range = Replace[
       OptionValue[Range],
       {(All | Full | Automatic) -> {{0, 180}, {0, 90}},
        r : {{_, _}, {_, _}} :> r,
        {t : {_, _}, r : Except[{_, _}]} :> {t, {0, r}},
        {t : {_, _}, (All | Full | Automatic)} :> {t, {0, 90}},
        {(All | Full | Automatic), r : {_, _}} :> {{0, 180}, r},
        {(All | Full | Automatic), 
          r : Except[{_, _}]} :> {{0, 180}, {0, r}},
        r_ :> {{0, 180}, {0, r}}}]},
    With[
      {msg = Which[
         ! NumericQ[range[[1, 1]]], "theta-min must be a number",
         ! NumericQ[range[[1, 2]]], "theta-max must be a number",
         ! NumericQ[range[[2, 1]]], "rho-min must be a number",
         ! NumericQ[range[[2, 2]]], "rho-max must be a number",
         ! (0 <= range[[1, 1]] < range[[1, 2]]), "theta-min must be in [0,theta-max)",
         range[[1, 2]] > 180, "theta-max must be in (theta-min,180]",
         ! (0 <= range[[2, 1]] < range[[2, 2]]), "rho-min must be in [0,rho-max)",
         range[[2, 2]] > 90, "rho-max must be in (rho-min,90]",
         True, None]},
      If[StringQ[msg], (Message[V123ParametricPlot::badarg, msg]; Throw[$Failed])]];
    With[
     {optseq = Sequence @@ FilterRules[{opts}, Options[ParametricPlot][[All,1]]],
      optseqShow = Sequence @@ FilterRules[{opts}, Options[Show][[All,1]]],
      rhoTrans = Function[90.0*#^3.5],
      rhoTransInv = Function[(#/90.0)^(1/3.5)],
      thetaLower = If[range[[1, 1]] > 90,
        None,
        {range[[1, 1]], Min[{90, range[[1, 2]]}]}],
      thetaUpper = If[range[[1, 2]] < 90,
        None,
        {Max[{90, range[[1, 1]]}], range[[1, 2]]}],
      f = If[hemi === LH, f0, Function@ReplaceAll[f0[##], {x_,y_}:>{-x,y}]]},
     With[
       {rMin = rhoTransInv[range[[2,1]]],
        rMax = rhoTransInv[range[[2,2]]]},
       With[
         {graphics = Map[
            Function@With[
              {k = Which[# == -4, 4, # == 4, 5, True, Abs[#]],
               thetaMinIdeal = If[# == 4 || # < 0, thetaLower[[1]], thetaUpper[[1]]],
               thetaMaxIdeal = If[# == -4 || # > 0, thetaUpper[[2]], thetaLower[[2]]]},
              With[
                {thetaMin = If[# == 2 || # == 3, 
                   Max[{thetaMinIdeal, 90.0 + epsilon}],
                   thetaMinIdeal],
                 thetaMax = If[# == -3 || # == -2,
                   Min[{thetaMaxIdeal, 90.0 - epsilon}],
                   thetaMaxIdeal]},
                Quiet[
                  ParametricPlot[
                    Part[f[theta, rhoTrans[rho]], k],
                    {theta, thetaMin, thetaMax},
                    {rho, rMin, rMax},
                    PlotRange -> plotRangeArg,
                    ColorFunction -> If[
                      Or[colorFun === Automatic,
                         colorFun === None,
                         colorFun === False,
                         StringQ[colorFun]],
                      colorFun,
                      If[colorFunSc === False,
                        Function[colorFun[#1,#2,#3,rhoTrans[#4]]],
                        Function[colorFun[#1,#2,#3,rhoTrans[#4]/90.0]]]],
                    optseq],
                  {CompiledFunction::cfsa}]]],
            areas]},
         With[
           {plotRange = With[
              {ranges = Cases[
                 AbsoluteOptions[#,PlotRange]& /@ graphics,
                 (PlotRange -> r_) :> r,
                 {2}]},
              {{Min[ranges[[All, 1, 1]]], Max[ranges[[All, 1, 2]]]},
               {Min[ranges[[All, 2, 1]]], Max[ranges[[All, 2, 2]]]}}]},
           Show[
             graphics,
             PlotRange -> plotRange, 
             optseqShow]]]]]]];


(* #SchiraLinePlot ********************************************************************************)
Options[SchiraLinePlot] = Join[
  FilterRules[Options[ParametricPlot], Except[PlotStyle | ColorFunction | ColorFunctionScaling]],
  {EccentricityStyleFunction -> Automatic,
   PolarAngleStyleFunction -> Automatic,
   EccentricityLines -> Automatic,
   PolarAngleLines -> Automatic,
   VisualAreas -> Automatic,
   Range -> Full,
   Hemi -> LH}];
SchiraLinePlot[mdl_SchiraModelObject, opts : OptionsPattern[]] := Catch[
  With[
    {epsilon = 0.000001,
     hemi = OptionValue[Hemi] /. {Left -> LH, Right -> RH, Automatic -> LH},
     plotRangeArg = OptionValue[PlotRange],
     areas = Union@Replace[
        OptionValue[VisualAreas],
        {All -> {-4, -3, -2, -1, 1, 2, 3, 4},
         Automatic -> {-3, -2, -1, 1, 2, 3},
         i_Integer /; -5 < i < 5 && i != 0 :> {i},
         l_List /; Length[l] == Count[l, i_Integer /; -5 < i < 5 && i != 0, {1}] :> Union[l],
         _ :> Message[
           SchiraLinePlot::badarg, 
           "VisualAreas must be All, one of +/- {1,2,3,4}, or a list of such integers"]}],
     fn0 = mdl[VisualFieldToCorticalMap],
     f = TemporarySymbol["f"],
     range = Replace[
       OptionValue[Range],
       {(All | Full | Automatic) -> {{0, 180}, {0, 90}},
        r : {{_, _}, {_, _}} :> r,
        {t : {_, _}, r : Except[{_, _}]} :> {t, {0, r}},
        {t : {_, _}, (All | Full | Automatic)} :> {t, {0, 90}},
        {(All | Full | Automatic), r : {_, _}} :> {{0, 180}, r},
        {(All | Full | Automatic), 
        r : Except[{_, _}]} :> {{0, 180}, {0, r}},
        r_ :> {{0, 180}, {0, r}}}],
     esf = Replace[
       OptionValue[EccentricityStyleFunction],
       x : (Automatic | Thick | Thin | Dotted | Dashed) :> With[
         {clr = With[{f = CorticalColorData["Eccentricity"]}, f[<|"Eccentricity" -> #|>]&]},
         If[x === Automatic, clr, {x, clr[#]} &]]],
     psf = Replace[
       OptionValue[PolarAngleStyleFunction],
       x : (Automatic | Thick | Thin | Dotted | Dashed) :> With[
         {clr = With[{f = CorticalColorData["PolarAngle"]}, f[<|"PolarAngle" -> #|>]&]},
         If[x === Automatic, clr, {x, clr[#]} &]]],
     eclines = Replace[
       OptionValue[EccentricityLines],
       {None -> {},
        x_?NumberQ :> {x},
        Automatic -> {1.25, 2.5, 5.0, 10.0, 20.0, 40.0, 90.0}}],
     palines = Replace[
       OptionValue[PolarAngleLines],
       {None -> {},
        x_?NumberQ :> {x},
        Automatic -> {0.0, 45.0, 90.0, 135.0, 180.0}}]},
    If[hemi === LH,
      f[t_?NumericQ, r_?NumericQ, k_Integer] := Part[fn[t, r], k],
      f[t_?NumericQ, r_?NumericQ, k_Integer] := ReplaceAll[Part[fn[t, r], k], {x_,y_}:>{-x,y}]];
    With[
      {msg = Which[
         ! NumericQ[range[[1, 1]]], "theta-min must be a number",
         ! NumericQ[range[[1, 2]]], "theta-max must be a number",
         ! NumericQ[range[[2, 1]]], "rho-min must be a number",
         ! NumericQ[range[[2, 2]]], "rho-max must be a number",
         ! (0 <= range[[1, 1]] < range[[1, 2]]), "theta-min must be in [0,theta-max)",
         range[[1, 2]] > 180, "theta-max must be in (theta-min,180]",
         ! (0 <= range[[2, 1]] < range[[2, 2]]), "rho-min must be in [0,rho-max)",
         range[[2, 2]] > 90, "rho-max must be in (rho-min,90]",
         True, None]},
      If[StringQ[msg], (Message[SchiraLinePlot::badarg, msg]; Throw[$Failed])]];
    With[
     {optseq = Sequence @@ FilterRules[
        {opts},
        ReplacePart[Options[ParametricPlot], {_, 2} -> _]],
      optseqShow = Sequence @@ FilterRules[
        {opts},
        ReplacePart[Options[Show], {_, 2} -> _]],
      rhoTrans = Function[{rho}, range[[2, 1]] + (range[[2, 2]] - range[[2, 1]])*rho^3.5], 
      thetaLower = If[range[[1, 1]] > 90, None, {range[[1, 1]], Min[{90, range[[1, 2]]}]}], 
      thetaUpper = If[range[[1, 2]] < 90, None, {Max[{90, range[[1, 1]]}], range[[1, 2]]}],
      angLines = Replace[
        areas,
        {(-4 | 4) :> palines,
         (-3 | -2) :> Select[palines, # <= 90 &] /. (90|90.0 -> 90.0 - epsilon),
         -1 :> Select[palines, # <= 90 &],
         (2 | 3) :> Select[palines, # >= 90 &] /. (90|90.0 -> 90.0 + epsilon),
         1 -> Select[palines, # >= 90 &]},
        {1}],
      eccLines = Table[eclines, {Length@areas}]},
     With[
       {areaIdcs = Which[# == -4, 4, # == 4, 5, True, Abs[#]] & /@ areas,
        thetaMinIdeals = If[# == 4 || # < 0, thetaLower[[1]], thetaUpper[[1]]] & /@ areas,
        thetaMaxIdeals = If[# == -4 || # > 0, thetaUpper[[2]], thetaLower[[2]]] & /@ areas},
       With[
         {thetaMins = MapThread[
            Function[
              If[#1 == 2 || #1 == 3, 
                Max[{#2, 90.0 + epsilon}], 
                #2]], 
            {areas, thetaMinIdeals}], 
          thetaMaxs = MapThread[
            Function[
              If[#1 == -3 || #1 == -2, 
                Min[{#2, 90.0 - epsilon}], 
                #2]],
            {areas, thetaMaxIdeals}]},
         With[
           {angPrep = Flatten[
              MapThread[
                Function[{lines, idx},
                  Map[
                    Function[
                      {Hold[
                         Evaluate[#1],
                         Evaluate[Block[{rho}, rhoTrans[rho]]],
                         Evaluate[idx]],
                       #1}],
                    lines]],
                {angLines, areaIdcs}],
              1],
            eccPrep = Flatten[
              MapThread[
                Function[{lines, min, max, idx},
                  Map[
                    Function[
                      {Hold[
                         Evaluate[Block[{theta}, min + (max - min)*theta]],
                         Evaluate[#1],
                         Evaluate[idx]],
                       #1}],
                    lines]],
                {eccLines, thetaMins, thetaMaxs, areaIdcs}],
              1]},
           With[
             {graphics = Flatten[
                {If[Length[angPrep] == 0,
                   {},
                   ReplacePart[
                     Hold[
                       Evaluate[ReplacePart[angPrep[[All, 1]], {_, 0} -> f]],
                       {rho, 0, 1},
                       Evaluate[PlotStyle -> psf /@ angPrep[[All, 2]]],
                       Evaluate[PlotRange -> plotRangeArg],
                       Evaluate[optseq]],
                     {0 -> ParametricPlot}]],
                 If[Length[eccPrep] == 0,
                   {},
                   ReplacePart[
                     Hold[
                       Evaluate[ReplacePart[eccPrep[[All, 1]], {_, 0} -> f]],
                       {theta, 0, 1},
                       Evaluate[PlotStyle -> esf /@ eccPrep[[All, 2]]],
                       Evaluate[PlotRange -> plotRangeArg],
                       Evaluate[optseq]],
                     {0 -> ParametricPlot}]]}]},
             With[
               {plotRange = If[plotRangeArg =!= Automatic && plotRangeArg =!= Full,
                  plotRangeArg,
                  With[
                    {ranges = Cases[Options /@ graphics, (PlotRange -> r_) :>r, {2}]},
                    {{Min[ranges[[All, 1, 1]]], Max[ranges[[All, 1, 2]]]},
                     {Min[ranges[[All, 2, 1]]], Max[ranges[[All, 2, 2]]]}}]]}, 
               Show[graphics, PlotRange -> plotRange, optseqShow]]]]]]]]];

(* #InterpretRangeOption (Private) ****************************************************************)
InterpretRangeOption[opt_] := Replace[
  opt,
  {(All | Full | Automatic) -> {{0, 180}, {0, 90}},
   r : {{_, _}, {_, _}} :> r,
   r : {_?NumericQ, _?NumericQ} :> {{0, 180}, r},
   {t : {_, _}, r : Except[{_, _}]} :> {t, {0, r}},
   {t : {_, _}, (All | Full | Automatic)} :> {t, {0, 90}},
   {(All | Full | Automatic), r : {_, _}} :> {{0, 180}, r},
   {(All | Full | Automatic), 
    r : Except[{_, _}]} :> {{0, 180}, {0, r}},
   r_ :> {{0, 180}, {0, r}}}];
Protect[InterpretRangeOption];

(* #InterpretVisualAreasOption (Private) **********************************************************)
InterpretVisualAreasOption[opt_] := With[
  {r = Union@Replace[
     opt,
     {All -> {1, 2, 3, 4, 5},
      Automatic -> {1, 2, 3},
      i_Integer :> {i},
      s_String :> List@Replace[
        s,
        {"V1" -> 1, "V2" -> 2, "V3" -> 3, "V4"|"hV4"|"HV4" -> 4, "V3a"|"V3A"|"V3AB"|"V3ab" -> 5}],
      l_List :> Replace[
        l,
        {"V1" -> 1, "V2" -> 2, "V3" -> 3, "V4"|"hV4"|"HV4" -> 4, "V3a"|"V3A"|"V3AB"|"V3ab" -> 5},
        {1}]}]},
  If[!ListQ[r] || Length@Complement[r, {1,2,3,4,5}] != 0,
    (Message[
       V123LinePlot::badarg, 
       "VisualAreas must be All, one of {V1, V2, V3, hV4, V4a} or a list of such ids"];
     $Failed),
    r]];
Protect[InterpretVisualAreasOption];

(* #V123IsoEccentricityLine ***********************************************************************)
Options[V123IsoEccentricityLine] = {
  Spacings -> Automatic,
  Range -> Full,
  VisualAreas -> Automatic};
V123IsoEccentricityLine[mdl_?V123ModelQ, ecc_?NumericQ, OptionsPattern[]] := With[
  {space = Replace[
     OptionValue[Spacings],
     Automatic :> Mean@PropertyValue[{mdl["Mesh"], 1}, MeshCellMeasure]],
   range = InterpretRangeOption@OptionValue[Range],
   areas = InterpretVisualAreasOption@OptionValue[VisualAreas]},
  With[
    {amin = range[[1,1]], amax = range[[1,2]]},
    With[
      {angs = If[Last[#] == amax, #, Append[#, amax]]&@Range[amin, amax, space],
       boolAreas = Table[MemberQ[areas, k], {k, 1, 5}]},
      With[
        {dat = VisualFieldToCorticalMap[mdl, angs, ConstantArray[ecc, Length[angs]]],
         pre90 = First@FirstPosition[angs, x_ /; x >= 90, {1}] - 1},
        Pick[
          {Reverse[dat[[ All,              4]]],
           dat[[         1 ;; pre90-1,     3]],
           Reverse[dat[[ 1 ;; pre90-1,     2]]],
           dat[[         All,              1]],
           Reverse[dat[[ pre90 + 2 ;; All, 2]]],
           dat[[         pre90 + 2 ;; All,  3]],
           Reverse[dat[[ All,              5]]]},
          boolAreas[[{4,3,2,1,2,3,5}]],
          True]]]]];
V123IsoEccentricityLine[mdl_?V123ModelQ, e_ /; VectorQ[e, NumericQ], opts:OptionsPattern[]] := Map[
  V123IsoEccentricityLine[mdl, #, opts]&,
  e];  
Protect[V123IsoEccentricityLine];

(* #V123IsoAngleLine ******************************************************************************)
Options[V123IsoAngleLine] = {Spacings -> Automatic, Range -> Full, VisualAreas -> Automatic};
V123IsoAngleLine[mdl_?V123ModelQ, ang_?NumericQ, OptionsPattern[]] := With[
  {space = Replace[
     OptionValue[Spacings],
     Automatic :> Mean@PropertyValue[{mdl["Mesh"], 1}, MeshCellMeasure]],
   range = InterpretRangeOption@OptionValue[Range],
   areas = InterpretVisualAreasOption@OptionValue[VisualAreas],
   exp = 9.0},
  With[
    {emin = range[[2,1]],
     emax = range[[2,2]]},
    With[
      {eccs = Join[Most[#], Union[{Last[#], emax}]]&@Range[emin, emax, space],
       boolAreas = Table[MemberQ[areas, k], {k, 1, 5}]},
      With[
        {dat = VisualFieldToCorticalMap[
           mdl,
           ConstantArray[ang, Length[eccs]],
           (((eccs - emin)/(emax - emin))^exp * (emax - emin) + emin)]},
        Pick[
          Transpose[dat],
          boolAreas[[{1,2,3,4,5}]],
          True]]]]];
V123IsoAngleLine[mdl_?V123ModelQ, a_Rule, opts:OptionsPattern[]] := V123IsoAngleLine[
  mdl, a[[1]], VisualAreas -> a[[2]], opts];
V123IsoAngleLine[mdl_?V123ModelQ, a_List, opts:OptionsPattern[]] := Map[
  Function@If[Head[#] === Rule,
    With[
      {r = V123IsoAngleLine[mdl, #[[1]], VisualAreas -> #[[2]], opts]},
      If[StringQ[#[[2]]] || IntegerQ[#[[2]]], r[[1]], r]],
    V123IsoAngleLine[mdl, #, opts]],
  a];
Protect[V123IsoAngleLine];

(* #V123LinePlot **********************************************************************************)
Options[V123LinePlot] = Join[
  FilterRules[Options[ParametricPlot], Except[PlotStyle | ColorFunction | ColorFunctionScaling]],
  {EccentricityStyleFunction -> Automatic,
   PolarAngleStyleFunction -> Automatic,
   EccentricityLines -> Automatic,
   PolarAngleLines -> Automatic,
   VisualAreas -> Automatic,
   Range -> Full}];
V123LinePlot[mdl_?AssociationQ, opts : OptionsPattern[]] := Catch[
  With[
    {epsilon = 0.01,
     areas = InterpretVisualAreasOption@OptionValue[VisualAreas],
     range = InterpretRangeOption@OptionValue[Range],
     esf = Replace[
       OptionValue[EccentricityStyleFunction],
       x : (Automatic | Thick | Thin | Dotted | Dashed) :> With[
         {clr = With[{f = CorticalColorData["Eccentricity"]}, f[<|"Eccentricity" -> #|>]&]},
         If[x === Automatic, clr, {x, clr[#]} &]]],
     psf = Replace[
       OptionValue[PolarAngleStyleFunction],
       x : (Automatic | Thick | Thin | Dotted | Dashed) :> With[
         {clr = With[{f = CorticalColorData["PolarAngle"]}, f[<|"PolarAngle" -> #|>]&]},
         If[x === Automatic, clr, {x, clr[#]} &]]],
     eclines = Replace[
       OptionValue[EccentricityLines],
       {None -> {},
        x_?NumberQ :> {x},
        Automatic -> {1.25, 2.5, 5.0, 10.0, 20.0, 40.0, 88.0}}],
     palines = Replace[
       OptionValue[PolarAngleLines],
       {None -> {},
        x_?NumberQ :> {x},
        Automatic -> {0.0, 90.0, 180.0}}]},
    With[
      {msg = Which[
         !NumericQ[range[[1, 1]]], "theta-min must be a number",
         !NumericQ[range[[1, 2]]], "theta-max must be a number",
         !NumericQ[range[[2, 1]]], "rho-min must be a number",
         !NumericQ[range[[2, 2]]], "rho-max must be a number",
         !(0 <= range[[1, 1]] < range[[1, 2]]), "theta-min must be in [0,theta-max)",
         range[[1, 2]] > 180, "theta-max must be in (theta-min,180]",
         !(0 <= range[[2, 1]] < range[[2, 2]]), "rho-min must be in [0,rho-max)",
         range[[2, 2]] > 90, "rho-max must be in (rho-min,90]",
         True, None]},
      If[StringQ[msg], (Message[V123LinePlot::badarg, msg]; Throw[$Failed])]];
    With[
      {rhoTrans = Function[{rho}, range[[2, 1]] + (range[[2, 2]] - range[[2, 1]])*rho^9.0], 
       thetaLower = If[range[[1, 1]] > 90, None, {range[[1, 1]], Min[{90, range[[1, 2]]}]}], 
       thetaUpper = If[range[[1, 2]] < 90, None, {Max[{90, range[[1, 1]]}], range[[1, 2]]}],
       angLines = Join@@Replace[
         areas,
         {a:1|4|5 :> Thread[palines -> a],
          2 :> Thread[# -> 2]&@Union@Fold[
            Function@If[MemberQ[areas, #2[[1]]], Pick[#1, Unitize[#1 - #2[[2]]], 1]],
            Join[
              palines /. (90|90.0 -> (90.0 - epsilon)),
              palines /. (90|90.0 -> (90.0 + epsilon))],
            {{1, 0}, {1, 180}}],
          3 :> Thread[# -> 3]&@Union@Fold[
            Function@If[MemberQ[areas, #2[[1]]], Pick[#1, Unitize[#1 - #2[[2]]], 1]],
            If[MemberQ[areas, 2],
              palines,
              Join[
                palines /. (90|90.0 -> (90.0 - epsilon)),
                palines /. (90|90.0 -> (90.0 + epsilon))]],
            {{2, 90}}]},
         {1}],
       eccLines = eclines,
       erange = ReplacePart[range, {2,2} -> Min[{range[[2,2]], Max[eclines]}]]},
      WithOptions[
        Graphics@Join[
          Map[
            {esf[#], Line@WithOptions[V123IsoEccentricityLine[mdl, #], opts]}&,
            eccLines],
          Map[
            {psf[#[[1]]], Line@WithOptions[V123IsoAngleLine[mdl, #], Range -> erange, opts]}&,
            angLines]],
        opts]]]];

Protect[SchiraParametricPlot, SchiraLinePlot,
        V123ParametricPlot, V123LinePlot,
        PolarAngleLegend, EccentricityLegend, VisualAreas,
        SchiraLinePlot, EccentricityStyleFunction, PolarAngleStyleFunction,
        EccentricityLines, PolarAngleLines];

(* #ManipulateSchiraModel *************************************************************************)
(* #here *)
ManipulateSchiraModel[map_?CorticalMapQ, model_SchiraModelObject :> var_Symbol] := $Failed;
Protect[ManipulateSchiraModel];

(* #GaussianSchiraPotential ***********************************************************************)
Options[GaussianSchiraPotential] = {
  VertexWeight -> Automatic,
  StandardDeviation -> 1.0,
  Print -> Subscript[Style["\[GothicCapitalG]", FontWeight -> Bold], "Schira"],
  MetaInformation -> {}};
GaussianSchiraPotential[map_?CorticalMapQ, model_SchiraModelObject, OptionsPattern[]] := With[
  {angle = VertexPropertyValues[map, "PolarAngle"],
   eccen = VertexPropertyValues[map, "Eccentricity"],
   weights = Replace[
     OptionValue[VertexWeight],
     {list_ /; VectorQ[list] && Length[list] == VertexCount[map] :> list,
      list:{(_Integer -> _)..} /; Length[list] == VertexCount[map] :> SparseArray[
        Select[VertexIndex[map, list[[All,1]]] -> list[[All, 2]], NumericQ[#[[2]]]&],
        VertexCount[map]],
      s_ /; VectorQ[VertexPropertyValues[map, s]] :> VertexPropertyValues[map, s],
      Automatic :> ConstantArray[1, VertexCount[map]],
      _ :> Message[GaussianSchiraPotential::badarg, "Unrecognized VertexWeight option"]}]},
  With[
    {idcs = Indices[
       Thread[{angle, eccen, weights}],
       {a_?NumericQ, e_?NumericQ /; e >= 0, w_ /; NumericQ[w] && Positive[w]}]},
    With[
      {indexedWeights = (#/Total[#])& @ weights[[idcs]],
       stddev = Replace[
         OptionValue[StandardDeviation],
         {x_?NumericQ :> ConstantArray[x, Length[idcs]],
          x_ /; VectorQ[x] && Length[x] == VertexCount[map] :> x[[idcs]],
          x_ /; VectorQ[x] && Length[x] == Length[idcs] :> x,
          s_ /; VectorQ[VertexPropertyValues[map, s]] :> Part[
            VertexPropertyValues[map, s]
            idcs],
          _ :> Message[
            GaussianSchiraPotential::badarg,
            StringJoin[
              "StandardDeviation option should be a number, a list of numbers, or a property name",
              " that is a list of numbers for all vertices"]]}],
       preds = Transpose[VisualFieldToCorticalMap[model, angle[[idcs]], eccen[[idcs]]], {3, 1, 2}]},
      Which[
        Length[idcs] == 0, Message[
          GaussianSchiraPotential::badarg,
          "No vertices selected"]; $Failed,
        !VectorQ[stddev, NumericQ[#] && Positive[#]&], Message[
          GaussianSchiraPotential::badarg,
          "standard deviations contains non-positive or non-numeric quantities"]; $Failed,
        True, With[
          {ks = indexedWeights / (Sqrt[2.0 * Pi] * stddev),
           vardenom = -1.0 / (2.0 * stddev^2),
           var = stddev^2,
           zeros = ConstantArray[0, {2, VertexCount[map]}]},
          CorticalPotentialFunction[
            {Function@Dot[
               ks,
               With[
                 {Xidcs = #[[All, idcs]]},
                 Total@Map[
                   (1.0 - Exp[Total[(Xidcs - #)^2] * vardenom]) &,
                   preds]]],
             Function@With[
               {Xidcs = #[[All, idcs]],
                dims = Length[#]},
               Total@Map[
                 Function@With[
                   {dX = Xidcs - #},
                   With[
                     {d = ColumnNorms[dX]},
                     With[
                       {grad = Times[
                          ConstantArray[ks / var * Exp[d^2 * vardenom], dims],
                          dX]},
                       Module[
                         {res = zeros},
                         res[[All, idcs]] = grad;
                         res]]]],
                 preds]],
             (* No hessian yet *)
             None},
            Print -> OptionValue[Print],
            CorticalMesh -> map,
            MetaInformation -> OptionValue[MetaInformation]]]]]]];
Protect[GaussianSchiraPotential];

(* #SchiraAnchors *********************************************************************************)
Options[SchiraAnchorsIndices] = {VertexWeight -> Automatic};
SchiraAnchorsIndices[map_?CorticalMapQ, OptionsPattern[]] := With[
  {angle = VertexPropertyValues[map, "PolarAngle"],
   eccen = VertexPropertyValues[map, "Eccentricity"],
   weights = Replace[
     OptionValue[VertexWeight],
     {list_ /; VectorQ[list] && Length[list] == VertexCount[map] :> list,
      list:{(_Integer -> _)..} /; Length[list] == VertexCount[map] :> SparseArray[
        Select[VertexIndex[map, list[[All,1]]] -> list[[All, 2]], NumericQ[#[[2]]]&],
        VertexCount[map]],
      s_ /; VectorQ[VertexPropertyValues[map, s]] :> VertexPropertyValues[map, s],
      Automatic :> Which[
        ListQ@VertexPropertyValues[map, "VertexWeight"], VertexPropertyValues[map, "VertexWeight"],
        ListQ@VertexPropertyValues[map, "Weight"], VertexPropertyValues[map, "Weight"],
        True, ConstantArray[1, VertexCount[map]]],
      _ :> Message[SchiraAnchors::badarg, "Unrecognized VertexWeight option"]}]},
  With[
    {idcs = Indices[
       Thread[{angle, eccen, weights}],
       {a_?NumericQ, e_?NumericQ /; e >= 0, w_ /; NumericQ[w] && Positive[w]}]},
    If[Length[idcs] == 0,
      (Message[SchiraAnchors::badarg, "No vertices selected"]; $Failed),
      idcs]]];

Options[SchiraAnchors] = Options[SchiraAnchorsIndices];
SchiraAnchors[map_?CorticalMapQ, model_SchiraModelObject, opts:OptionsPattern[]] := With[
  {angle = VertexPropertyValues[map, "PolarAngle"],
   eccen = VertexPropertyValues[map, "Eccentricity"],
   idcs = SchiraAnchorsIndices[map, opts]},
  With[
    {preds = VisualFieldToCorticalMap[model, angle[[idcs]], eccen[[idcs]]]},
    {Join@@ConstantArray[VertexList[map][[idcs]], 5],
     Transpose[Join @@ Transpose[preds]]}]];

Options[SchiraAnchorsParameter] = Options[SchiraAnchors];
SchiraAnchorsParameter[map_?CorticalMapQ, paramArg_, opts:OptionsPattern[]] := With[
  {param = Which[
     StringQ[paramArg], VertexPropertyValues[map, paramArg],
     ListQ[paramArg], paramArg,
     True, Message[SchiraAnchors::badarg, "parameter must be a list or a property name"]],
   idcs = SchiraAnchorsIndices[map, opts]},
  Join@@ConstantArray[
    Which[
      Length[param] == VertexCount[map], param[[idcs]],
      Length[param] == Length[idcs], param,
      True, Message[SchiraAnchors::badarg, "given parameter is wrong size"]],
    5]];

Protect[SchiraAnchors, SchiraAnchorsParameter, SchiraAnchorsIndices];

(* #V123Anchors ***********************************************************************************)
Options[V123AnchorsIndices] = {VertexWeight -> Automatic};
V123AnchorsIndices[map_?CorticalMapQ, OptionsPattern[]] := With[
  {angle = VertexPropertyValues[map, "PolarAngle"],
   eccen = VertexPropertyValues[map, "Eccentricity"],
   weights = Replace[
     OptionValue[VertexWeight],
     {list_ /; VectorQ[list] && Length[list] == VertexCount[map] :> list,
      list:{(_Integer -> _)..} /; Length[list] == VertexCount[map] :> SparseArray[
        Select[VertexIndex[map, list[[All,1]]] -> list[[All, 2]], NumericQ[#[[2]]]&],
        VertexCount[map]],
      s_ /; VectorQ[VertexPropertyValues[map, s]] :> VertexPropertyValues[map, s],
      Automatic :> Which[
        ListQ@VertexPropertyValues[map, "VertexWeight"], VertexPropertyValues[map, "VertexWeight"],
        ListQ@VertexPropertyValues[map, "Weight"], VertexPropertyValues[map, "Weight"],
        True, ConstantArray[1, VertexCount[map]]],
      _ :> Message[V123Anchors::badarg, "Unrecognized VertexWeight option"]}]},
  With[
    {idcs = Indices[
       Thread[{angle, eccen, weights}],
       {a_?NumericQ, e_?NumericQ /; e >= 0, w_ /; NumericQ[w] && Positive[w]}]},
    If[Length[idcs] == 0,
      (Message[V123Anchors::badarg, "No vertices selected"]; $Failed),
      idcs]]];

Options[V123Anchors] = Options[V123AnchorsIndices];
V123Anchors[map_?CorticalMapQ, model_?AssociationQ, opts:OptionsPattern[]] := With[
  {angle = VertexPropertyValues[map, "PolarAngle"],
   eccen = VertexPropertyValues[map, "Eccentricity"],
   idcs = V123AnchorsIndices[map, opts]},
  With[
    {preds = VisualFieldToCorticalMap[model, angle[[idcs]], eccen[[idcs]]]},
    {Join@@ConstantArray[VertexList[map][[idcs]], 5],
     Transpose[Join @@ Transpose[preds]]}]];

Options[V123AnchorsParameter] = Options[V123Anchors];
V123AnchorsParameter[map_?CorticalMapQ, paramArg_, opts:OptionsPattern[]] := With[
  {param = Which[
     StringQ[paramArg], VertexPropertyValues[map, paramArg],
     ListQ[paramArg], paramArg,
     True, Message[SchiraAnchors::badarg, "parameter must be a list or a property name"]],
   idcs = V123AnchorsIndices[map, opts]},
  Join@@ConstantArray[
    Which[
      Length[param] == VertexCount[map], param[[idcs]],
      Length[param] == Length[idcs], param,
      True, Message[V123Anchors::badarg, "given parameter is wrong size"]],
    5]];

Protect[V123Anchors, V123AnchorsParameter, V123AnchorsIndices];


(* #V123MeshFunctions (Private) *******************************************************************)
V123MeshFunctions[data_, tx0_:Automatic] := With[
  {mesh = data["Mesh"],
   coords = MeshCoordinates@data["Mesh"],
   cells = MeshCells[data["Mesh"], 2][[All, 1]],
   angle = data["PolarAngle"],
   eccen = data["Eccentricity"],
   iregions = {"V1", "V2", "V3", "V4Ventral", "V4Dorsal"},
   tx = Which[
     tx0 === None, Identity,
     tx0 === Automatic, With[
       {tx1 = If[AssociationQ[#], #["AffineTransform"], None]&@data["Options"]},
       If[tx1 === None, Identity, tx1]],
     True, tx0]},
  With[
    {field = eccen*Exp[I*angle],
     ifieldAll = Flatten[coords.{{1}, {I}}],
     itx = InverseFunction[tx]},
    With[
      {xy = Transpose[{Re[field], Im[field]}]},
      With[
        {imesh = With[
           {F = Transpose@Select[cells, Length@Union[xy[[#]]] == 3 &]},
           Association@Table[
             id -> With[
               {bound = data[id <> "BoundaryMesh"]},
               MeshRegion[
                 xy,
                 Polygon /@ Pick[
                   Transpose[F],
                   Total@Unitize[RegionDistance[bound, coords[[#]]] & /@ F],
                   0]]],
             {id, iregions}]],
         near = Association@Table[
           id -> With[
             {idcs = Pick[
                Range@Length[coords],
                Unitize@RegionDistance[data[id<>"BoundaryMesh"], coords],
                0]},
             Nearest[xy[[idcs]] -> idcs]],
           {id, iregions}]},
        With[
          {ifields = Association@Table[
             reg -> ifieldAll[[near[reg][MeshCoordinates@imesh[reg], 1][[All, 1]]]],
             {reg, iregions}]},
          Join[
            data,
            <|"CorticalMapToVisualField" -> Function@Which[
                MatrixQ[#], MeshRegionInterpolate[mesh, field, itx[#], Chop -> 0.1],
                VectorQ[#], First@MeshRegionInterpolate[mesh, field, itx[{#}], Chop -> 0.1],
                True,       $Failed],
              "VisualFieldToCorticalMap" -> Function@With[
                {dat = If[VectorQ[#], #, {#}],
                 post = If[VectorQ[#], Identity, #[[All,1]]&]},
                post@Transpose@Table[
                  tx@ReIm@MeshRegionInterpolate[
                    imesh[reg], ifields[reg], ReIm[dat],
                    Chop -> Infinity],
                  {reg, iregions}]],
              "VisualFieldMeshes" -> imesh,
              "VisualFieldFields" -> ifield|>]]]]]];
Protect[V123MeshFunctions];

(* MeshRegionOrthogonalFields *********************************************************************)
Options[MeshRegionFields] = Options[FindArgMin];
MeshRegionFields[region_, initF1Arg_List, initF2Arg_List, opts:OptionsPattern[]] := With[
  {initF1 = Select[Union /@ GatherBy[initF1Arg, First], Length[#] == 1 &][[All, 1]],
   initF2 = Select[Union /@ GatherBy[initF2Arg, First], Length[#] == 1 &][[All, 1]]},
  With[
    {f1Vertices = initF1[[All, 1]],
     f2Vertices = initF2[[All, 1]],
     xT = Transpose@MeshCoordinates[region],
     edgesT = Transpose[MeshCells[region, 1] /. Line[idcs_] :> idcs],
     n = MeshCellCount[region, 0],
     m = MeshCellCount[region, 1]},
    With[
      {sumOverEdgesT = Transpose@SparseArray@MapThread[
         {#1, #2} -> #3 &,
         {Join @@ edgesT,
          Join[#, #] &@Range[m],
          Join[ConstantArray[1, m], ConstantArray[-1, m]]}],
       zeroTh = ReplacePart[initF1, {_, 2} -> 0],
       zeroR = ReplacePart[initF2, {_, 2} -> 0],
       x0 = Join @@ MapThread[
         ReplacePart,
         {{ConstantArray[Mean@initF1[[All, 2]], n],
           ConstantArray[Mean@initF2[[All, 2]], n]},
          {initF1, initF2}}]},
      Block[
        {f, x},
        (* The Cost Function *)
        f[z0_ /; VectorQ[z0, NumericQ]] := With[
          {z = MapThread[ReplacePart, {Partition[z0, n], {initF1, initF2}}]},
          Plus[
            (* The Squared-Difference between Edge-Values: *)
            0.5*Total@Total[(z[[All, edgesT[[1]]]] - z[[All, edgesT[[2]]]])^2],
            (* The Orthogonality of the Fields: *)
            0.5*(Dot @@ (z - Mean /@ z))^2]];
        (* The Gradient Function *)
        f /: Grad[f, z0_ /; VectorQ[z0, NumericQ]] := Join @@ With[
          {z = Partition[z0, n]},
          MapThread[ReplacePart, {#, {zeroTh, zeroR}}] &@Plus[
            (* The z-values want to approach each other like springs: *)
            (z[[All, edgesT[[1]]]] - z[[All, edgesT[[2]]]]).sumOverEdgesT,
            (* The z-values want to push away from each other *)
            (#*(Dot@@#))&@Reverse[z - Mean /@ z]]];
        Partition[First@FindArgMin[f[x], {x, x0}, Gradient :> Grad[f, x], opts], n]]]]];
Protect[MeshRegionOrthogonalFields];

(* #V123Model *************************************************************************************)
Options[V123Model] = {
  Scale -> {15, 5, 5, 3},
  AspectRatio -> 0.32,
  Intersection -> 0.6*Pi,
  Function -> Function[0.06*# + #^14],
  MaxValue -> 90.0,
  AffineTransform -> AffineTransform[{RotationMatrix[5 Degree].{{1, -0.2}, {0, 1}}, {-7, -1}}],
  CellSize -> 0.1,
  MaxIterations -> 10000};
V123Model[OptionsPattern[]] := With[
  {assc = Association@Table[ToString[k] -> OptionValue[k], {k, Options[V123Model][[All, 1]]}]},
  assc // Function@With[
    {ellipseCenter = {0.5*#Scale[[1]]/#AspectRatio, 0},
     ellipseAxes = {#Scale[[1]]/#AspectRatio, #Scale[[1]]}/2,
     v1w = #Scale[[1]],
     v2w = #Scale[[2]],
     v3w = #Scale[[3]],
     v4w = #Scale[[4]],
     tx = If[#AffineTransform === None, Identity, #AffineTransform],
     itx = If[#AffineTransform === None, Identity, InverseFunction[#AffineTransform]],
     cellsz = #CellSize},
    With[
      {ipoint = ellipseCenter[[1]] + 0.5*v1w*Tan[Pi/2 - #Intersection/2],
       cpoint = ellipseCenter[[1]] - 0.5*v1w*Tan[#Intersection/2],
       rad = v1w/(2*Cos[#Intersection/2]),
       v2axes = ellipseAxes + v2w*{0.5, 1},
       v3axes = ellipseAxes + (v2w + v3w)*{0.5, 1},
       v4axes = ellipseAxes + (v2w + v3w + v4w)*{0.5, 1},
       xmin = ellipseCenter[[1]] - (ellipseAxes[[1]] + 0.5*(v2w + v3w + v4w))},
      With[
        {tri = Triangle[{{ipoint,0}, {xmin,#}, {xmin,-#}}]&[(ipoint - xmin)*Tan[#Intersection/2]],
         circle = Disk[{cpoint, 0}, rad],
         v1ellipse = Ellipsoid[ellipseCenter, ellipseAxes],
         v2ellipse = Ellipsoid[ellipseCenter, v2axes],
         v3ellipse = Ellipsoid[ellipseCenter, v3axes],
         v4ellipse = Ellipsoid[ellipseCenter, v4axes]},
        With[
          {utri = Normalize[tri[[1, 2]] - tri[[1, 1]]],
           fixfn = Function@With[
             {x = #1, v = #2},
             Pick[
               Join[Reverse[#], ({#[[1, 1]], -#[[1, 2]]} -> v) & /@ x],
               Join[Reverse[#], #]&@Chop@RegionDistance[tri, x[[All, 1]]],
               0]]},
          With[
            {xpts = Pick[#, Chop@SignedRegionDistance[tri, #[[All, 1]]], x_ /; x <= 0]&@Join[
               Table[{-x, 0} -> {If[x == 0, 0, _], 0}, {x, 0, -xmin, cellsz}],
               Table[{x, 0} -> {0, _}, {x, cellsz, cpoint + rad, cellsz}]],
             v1v2pts = fixfn[#, {Pi/2, _}] &@Table[
               (ellipseCenter + {-Cos[t], Sin[t]}*ellipseAxes) -> {-Pi/2, _},
               {t, 0, Pi/2, cellsz/Mean[ellipseAxes]}],
             v2v3pts = fixfn[#, {0, _}] &@Table[
               (ellipseCenter + {-Cos[t], Sin[t]}*(v2axes)) -> {0, _},
               {t, 0, Pi/2, cellsz/Mean[v2axes]}],
             v3v4pts = fixfn[#, {Pi/2, _}] &@Table[
               (ellipseCenter + {-Cos[t], Sin[t]}*(v3axes)) -> {-Pi/2, _},
               {t, 0, Pi/2, cellsz/Mean[v3axes]}],
             v4pts = fixfn[#, {-Pi/2, _}] &@Table[
               (ellipseCenter + {-Cos[t], Sin[t]}*(v4axes)) -> {Pi/2, _},
               {t, 0, Pi/2, cellsz/Mean[v4axes]}],
             peripts = Union@With[
               {possibles = Select[
                  Join[
                    Table[
                      ({ipoint, 0} + k*utri) -> {_, 1},
                      {k, 0, Norm[{(ipoint - xmin), tri[[1, 2, 2]]}], #CellSize}],
                    Table[
                      ({ipoint, 0} + k*{utri[[1]], -utri[[2]]}) -> {_, 1},
                      {k, 0, Norm[{(ipoint - xmin), tri[[1, 2, 2]]}], #CellSize}]],
                  Abs[#[[1, 2]]] > v1w/2 &]},
               Join[
                 Pick[possibles, Chop@RegionDistance[v4ellipse, possibles[[All, 1]]], 0],
                 Join[# -> {_, 1} & /@ #, {#[[1]], -#[[2]]} -> {_, 1} & /@ #]&@Table[
                   {cpoint, 0} + {Cos[t], Sin[t]}*rad,
                   {t, 0, (Pi - #Intersection)/2, #CellSize/rad}]]]},
            With[
              {b = BoundaryMesh@DelaunayMesh[Join[v4pts, peripts][[All, 1]]],
               v1peripts = Pick[
                 peripts[[All, 1]],
                 Chop@RegionDistance[v1ellipse, peripts[[All, 1]]],
                 0],
               v2peripts0 = Pick[
                 peripts[[All, 1]],
                 Chop@RegionDistance[v2ellipse, peripts[[All, 1]]],
                 0],
               v3peripts0 = Pick[
                 peripts[[All, 1]],
                 Chop@RegionDistance[v2ellipse, peripts[[All, 1]]],
                 0],
               pickfn = Function@Pick[#1, Chop@RegionDistance[#2, #1], 1],
               yExtreme = (v1w + v2w + v3w + v4w)},
              With[
                {points = With[
                   {fillpts = Join @@ MapIndexed[
                      Function@With[
                        {ymod = If[Mod[#2[[1]], 2] == 0, cellsz/2, 0]},
                        Table[
                          {#1, y} -> {_, _},
                          {y, -yExtreme - ymod, yExtreme + ymod, cellsz}]],
                      Range[xmin, cpoint + rad, cellsz]],
                    boundpts = Union[xpts, v1v2pts, v2v3pts, v3v4pts, v4pts, peripts]},
                   Map[
                     Function@With[
                       {pa = Union@Select[#[[All, 2, 1]], # =!= _ &],
                        ec = Union@Select[#[[All, 2, 2]], # =!= _ &]},
                       #[[1, 1]] -> If[Length[#] == 1,
                         #[[1, 2]],
                         {If[pa == {}, _, Mean[pa]],
                          If[ec == {}, _, Mean[ec]]}]],
                     GatherBy[#, First] &@SortBy[#, First] &@N@Union[
                       boundpts,
                       Intersection[
                         Delete[
                           fillpts,
                           Transpose[
                             {Union@@Nearest[
                                fillpts[[All, 1]] -> Automatic,
                                boundpts[[All, 1]],
                                {Infinity, cellsz}]}]],
                         Pick[fillpts, Chop@RegionDistance[b, fillpts[[All, 1]]], 0]]]]],
                 v2peripts = pickfn[v2peripts0, v1ellipse],
                 v3peripts = pickfn[v3peripts0, v2ellipse],
                 v4peripts = pickfn[peripts[[All, 1]], v3ellipse]},
                With[
                  {mesh = With[
                     {mesh0 = Quiet[ToElementMesh[points[[All, 1]]], {ToElementMesh::femimq}]},
                     With[
                       {X = mesh0["Coordinates"],
                        F = mesh0["MeshElements"][[1, 1]],
                        q = mesh0["Quality"]},
                       MeshRegion[X, Polygon /@ Delete[F, Transpose@List@Indices[Sign[q], -1|0]]]]],
                   angc = MapIndexed[
                     If[NumericQ[#1[[2, 1]]], #2[[1]] -> #1[[2, 1]], Nothing] &,
                     points],
                   eccc = MapIndexed[
                     If[NumericQ[#1[[2, 2]]], #2[[1]] -> #1[[2, 2]], Nothing] &,
                     points],
                   v2peri = Pick[v2peripts, Sign[v2peripts[[All, 2]]], #] & /@ {1, -1},
                   v3peri = Pick[v3peripts, Sign[v3peripts[[All, 2]]], #] & /@ {1, -1},
                   v4peri = Pick[v4peripts, Sign[v4peripts[[All, 2]]], #] & /@ {1, -1},
                   bmeshfn = Function@With[
                     {x = Keys@GroupBy[#, Identity]},
                     BoundaryMeshRegion[x, Line@Append[Range@Length[x], 1]]]},
                  With[
                    {near = Nearest[MeshCoordinates[mesh] -> Automatic],
                     Xang = points[[angc[[All, 1]], 1]],
                     Xecc = points[[eccc[[All, 1]], 1]],
                     v1b = BoundaryMesh@DelaunayMesh@Union[v1v2pts[[All, 1]], v1peripts],
                     v2b = bmeshfn@Join[
                       v1v2pts[[All, 1]], Reverse[v2peri[[2]]],
                       Reverse[v2v3pts[[All, 1]]], v2peri[[1]]],
                     v3b = bmeshfn@Join[
                       v2v3pts[[All, 1]], Reverse[v3peri[[2]]],
                       Reverse[v3v4pts[[All, 1]]], v3peri[[1]]],
                     v4b = bmeshfn@Join[
                       v3v4pts[[All, 1]], Reverse[v4peri[[2]]],
                       Reverse[v4pts[[All, 1]]], v4peri[[1]]]},
                    With[
                      {angConst = Thread[# -> angc[[All, 2]]]&[near[Xang, 1][[All, 1]]],
                       eccConst = Thread[# -> eccc[[All, 2]]]&[near[Xecc, 1][[All, 1]]],
                       v4vb = bmeshfn[Gather[N@#][[All, 1]]] &@Join[
                         Select[v3v4pts[[All, 1]], #[[2]] <= 0 &],
                         Reverse[v4peri[[2]]],
                         Select[Reverse[v4pts[[All, 1]]], #[[2]] <= 0 &],
                         Reverse@Pick[
                           xpts[[All, 1]],
                           Sign@RegionDistance[v4b, xpts[[All, 1]]],
                           0]],
                       v4db = bmeshfn[Gather[N@#][[All, 1]]]&@Join[
                         Select[v3v4pts[[All, 1]], #[[2]] >= 0 &],
                         Reverse[v4peri[[1]]],
                         Select[Reverse[v4pts[[All, 1]]], #[[2]] >= 0 &]]},
                      With[
                        {fields = MeshRegionFields[
                           mesh,
                           angConst, eccConst,
                           MaxIterations -> #MaxIterations],
                         eccfn = With[
                           {f = #Function},
                           With[{mn = f[0], mx = f[1]}, Function[(f[#] - mn) / (mx - mn)]]]},
                        With[
                          {data = <|
                           "Mesh" -> mesh,
                           "Options" -> #,
                           "PolarAngle" -> fields[[1]],
                           "Eccentricity" -> (#MaxValue * eccfn[fields[[2]]]),
                           "BoundaryMesh" -> b,
                           "V1BoundaryMesh" -> v1b,
                           "V2BoundaryMesh" -> v2b,
                           "V3BoundaryMesh" -> v3b,
                           "V4BoundaryMesh" -> v4b,
                           "V4VentralBoundaryMesh" -> v4vb,
                           "V4DorsalBoundaryMesh" -> v4db,
                           "PolarAngleConstraints" -> angConst,
                           "EccentricityContstraints" -> eccConst,
                           "VisualAreaFunction" -> Function@With[
                             {x0 = itx[#]},
                             With[
                               {inV1  = 1 - Unitize@RegionDistance[v1b,  x0],
                                inV2  = 1 - Unitize@RegionDistance[v2b,  x0],
                                inV3  = 1 - Unitize@RegionDistance[v3b,  x0],
                                inV4v = 1 - Unitize@RegionDistance[v4vb, x0],
                                inV4d = 1 - Unitize@RegionDistance[v4db, x0]},
                               1*inV1 + 2*inV2 + 3*inV3 + 4*inV4v + 5*inV4d]]|>},
                          V123MeshFunctions[data, tx]]]]]]]]]]]]]];
Protect[V123Model];

(* #V123ModelQ ************************************************************************************)
V123ModelQ[_] := False;
V123ModelQ[mdl_?AssociationQ] := And[
  KeyExistsQ[mdl, "VisualFieldToCorticalMap"],
  KeyExistsQ[mdl, "CorticalMapToVisualField"],
  KeyExistsQ[mdl, "Mesh"],
  KeyExistsQ[mdl, "PolarAngle"],
  KeyExistsQ[mdl, "Eccentricity"]];
Protect[V123ModelQ];

End[];
EndPackage[];

