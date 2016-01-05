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
  {"Neurotica`Global`", 
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
CorticalMapToRetinotopy::usage = "CorticalMapToRetinotopy[model, map] yields a list of predictions, one per vertex in map, of the {polar angle, eccentricity, visual area} for the vertex in the given SchiraModelObject model. For any vertex that lies outside of the model's bounds, the absolute value of the visual area will be greater than 4. A list of vertices or a single vertex may also be substituted for map.
CorticalMapToRetinotopy[model, X, Y] is equivalent to CorticalMapToRetinotopy[model, Thread[{X,Y}]].
CorticalMapToRetinotopy[model] yields a pure curried function that may be called with map directly.";
RetinotopyToCorticalMap::usage = "RetinotopyToCorticalMap[model, retinotopy] yields a list of 5 x 2 matrices, each row of which gives the {x, y} coordinate predictions for one of the visual areas. The result contains one such matrix for each retinotopic coordinate given. The retinotopy argument must be a list of {polarAngle, eccentricity}. The rows of each coordinate matrix returned represent, in order, the V1, V2, V3, HV4, and V3A predictions for the given retinotopic coordinate. A single retinotopy coordinate may be given, in which case, a single coordinate matrix is returned.
RetinotopyToCorticalMap[model, polarAngles, eccentricities] is equivalent to RetinotopyToCorticalMap[model, Thread[{polarAngles, eccentricitie}]].
RetinotopyToCorticalMap[model] yields a pure curried function that may be called with retinotopy directly.";

PolarAngleLegend::usage = "PolarAngleLegend[hemi] yields a graphic that is appropriate for a polar angle legend. All options that are valid for DensityPlot can be passed.";
EccentricityLegend::usage = "EccentricityLegend[hemi,max] yields a graphic that is appropriate for an eccentricity legend. All options that are valid for DensityPlot can be passed. The range is an optional argument that specifies the max eccentricity for this legend (default: 20)";

SchiraParametricPlot::usage = "SchiraParametricPlot[mdl, options...] is equivalent to ParametricPlot[f[th,r], {th, thMin, thMax}, {r, rMin, rMax}, options] where f is the Schira function for the given SchiraModelObject mdl, and thMin, thMax, rMin, and rMax can be controlled via the additional option Range. The visual areas specified by the VisualAreas argument are plotted. The Range argument may be of the following forms: {thRange, rRange} or rRange where thRange must be {thMin, thMax} and rRange may be rMax (with rMin 0) or {rMin, rMax}.";
SchiraParametricPlot::badarg = "Bad argument to SchiraParametricPlot: `1`";
VisualAreas::usage = "VisualAreas is an optional keyword argument to the SchiraParametricPlot function. VisualAreas may be a list of any of {-4,-3,-2,-1,1,2,3,4}, which represent areas hV4, V3ventral, V2ventral, V1Ventral, V1Dorsal, V2Dorsal, V3Dorsal, V3A, respectively. The default value or Automatic will yield {-3,-2,-1,1,2,3}.";

SchiraLinePlot::usage = "SchiraLinePlot[mdl,options...] is equivalent to calling ParametricPlot[] with the given options, but such that the first arguments are optimized to draw polar angle and eccentricity lines from the given SchiraModelObject mdl. The following additional parameters may be given:
PolarAngleLines gives the polar angle values at which to daw iso-eccentric lines.
EccentricityLines gives the eccentricity values at which to draw iso-angular lines.
PolarAngleStyleFunction gives a function that, given a polar angle value, yields the style instruction for that polar angle line.
EccentricityFunction gives a function that, given a polar angle value, yields the style instruction for that polar angle line.
VisualAreas is used as in SchiraParametricPlot.
Range is used as in SchiraParametricPlot.";
SchiraLinePlot::badarg = "Bad argument to SchiraLinePlot: `1`";
PolarAngleLines::usage = "PolarAngleLines is an option to SchiraLinePlot that specifies the polar angle values at which to draw the iso-eccentric lines.";
EccentricityLines::usage = "EccentricityLines is an option to SchiraLinePlot that specifies the eccentricity values at which to draw the iso-angular lines.";
PolarAngleStyleFunction::usage = "PolarAngleStyleFunction is an option to SchiraLinePlot that must accept a polar angle value and yield the style directive for plotting that angle's iso-eccentric line. The value Automatic results in ColorCortex[PolarAngle] beign used. The values Thick, Dashed, Dotted, etc. will result in the same coloring schele with the option directive appended.";
EccentricityStyleFunction::usage = "EccentricityStyleFunction is an option to SchiraLinePlot that must accept an eccentricity value and yield the style directive for plotting that eccentricity's iso-angular line. The value Automatic results in ColorCortex[Eccentricity] beign used. The values Thick, Dashed, Dotted, etc. will result in the same coloring schele with the option directive appended.";

ManipulateSchiraModel::usage = "ManipulateSchiraModel[map] yields a manipulate form that allows one to edit the parameters of the Schira model which is projected over a cortex plot of the given map.
ManipulateSchiraModel[map, model] uses the given model as the starting parameterization of the model.
ManipulateSchiraModel[map, model :> var] assigns the updated model value to the given var every time a change is entered, assuming that the \"Save on change?\" option is set to True.";

(* Registration potential functions ***************************************************************)
GaussianSchiraPotential::usage = "GaussianSchiraPotential[map, model] yields a potential function that describes the agreement between the given SchiraModel object, model, and the retinotopy data found in the properties \"PolarAngle\" and \"Eccentricity\" of map. If the property \"VertexWeight\" is also present, then all vertices with a value above 0 are included and they are weighted by the given weights. Otherwise, any vertex with missing polar angle or eccentricity values ($Failed, None, or Indeterminate) are excluded.";
GaussianSchiraPotential::badarg = "Bad argument given to GaussianSchiraPotential: `1`";

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
              Abs[a] < 1 || Abs[a] > 3, {None, None, None},
              NumberQ[pa] && pa > 90.0, {pa, e, -Abs[a]},
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
           CorticalMapToRetinotopy :> prf,
           RetinotopyToCorticalMap :> prc,
           All :> params,
           x_ :> Message[SchiraModelObject::badarg, x]}]]]},
    prf := With[
      {res = Check[CorticalMapToRetinotopy[mdl], $Failed]},
      If[res === $Failed, res, (prf = res)]];
    prc := With[
      {res = Check[RetinotopyToCorticalMap[mdl], $Failed]},
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

(* #CorticalMapToRetinotopy ***********************************************************************)
CorticalMapToRetinotopy[SchiraModelObject[disp_], map_?CorticalMapQ] := With[
  {inv = Replace[Inverse, disp],
   Z = VertexCoordinatesTr[map],
   r90 = Replace[\[CapitalRho]90, disp]},
  Quiet[
    Map[
      Append[ComplexToVisualAngle[#[[1]]], #[[2]]]&,
      inv[(Z[[1]] + I * Z[[2]])] * (90.0 / r90)],
    {FindRoot::cvmit, FindRoot::frmp}]];
CorticalMapToRetinotopy[SchiraModelObject[disp_], {x:Except[_List], y:Except[_List]}] := With[
  {inv = Replace[Inverse, disp],
   r90 = Replace[\[CapitalRho]90, disp]},
  With[
   {z = Quiet[inv[x + I*y] * (90.0 / r90), {FindRoot::cvmit, FindRoot::frmp}]},
   Append[ComplexToVisualAngle[z[[1]]], z[[2]]]]];
CorticalMapToRetinotopy[SchiraModelObject[disp_],
                        {x_List, y_List} /; Length[x] == Length[y] && Length[x] != 2] := With[
  {inv = Replace[Inverse, disp],
   r90 = Replace[\[CapitalRho]90, disp]},
  With[
    {res = Quiet[inv[X + I*Y] * (90.0 / r90), {FindRoot::cvmit, FindRoot::frmp}]},
    If[ListQ[First@res],
      Append[ComplexToVisualAngle[#[[1]]], #[[2]]]& /@ res,
      Append[ComplexToVisualAngle[res[[1]]], res[[2]]]]]];
CorticalMapToRetinotopy[SchiraModelObject[disp_], coords:{{_,_}..}] := With[
  {inv = Replace[Inverse, disp],
   Z = Transpose[coords],
   r90 = Replace[\[CapitalRho]90, disp]},
  Map[
    Append[ComplexToVisualAngle[#[[1]]], #[[2]]]&,
    Quiet[inv[Z[[1]] + I * Z[[2]]] * (90.0 / r90), {FindRoot::cvmit, FindRoot::frmp}]]];
CorticalMapToRetinotopy[SchiraModelObject[disp_], X_, Y_] := With[
  {inv = Replace[Inverse, disp],
   r90 = Replace[\[CapitalRho]90, disp]},
  With[
    {res = Quiet[inv[X + I*Y] * (90.0 / r90), {FindRoot::cvmit, FindRoot::frmp}]},
    If[ListQ[First@res],
      Append[ComplexToVisualAngle[#[[1]]], #[[2]]]& /@ res,
      Append[ComplexToVisualAngle[res[[1]]], res[[2]]]]]];
CorticalMapToRetinotopy[mdl_SchiraModelObject] := Function[CorticalMapToRetinotopy[mdl, ##]];

RetinotopyToCorticalMap[SchiraModelObject[disp_], retinotopy:{{_,_}..}] := With[
  {fun = Replace[Function, disp],
   tr = Transpose[retinotopy],
   r90 = Replace[\[CapitalRho]90, disp]},
  ComplexToCoordinate[fun[r90 / 90.0 * VisualAngleToComplex[tr[[1]], tr[[2]]]]]];
RetinotopyToCorticalMap[
  SchiraModelObject[disp_],
  {polarAngle:Except[_List], eccentricity:Except[_List]}
 ] := With[
  {fun = Replace[Function, disp],
   r90 = Replace[\[CapitalRho]90, disp]},
  ComplexToCoordinate[fun[r90 / 90.0 * VisualAngleToComplex[polarAngle, eccentricity]]]];
RetinotopyToCorticalMap[SchiraModelObject[disp_], polarAngles_, eccentricities_] := With[
  {fun = Replace[Function, disp],
   r90 = Replace[\[CapitalRho]90, disp]},
  ComplexToCoordinate[fun[r90 / 90.0 * VisualAngleToComplex[polarAngles, eccentricities]]]];
RetinotopyToCorticalMap[mdl_SchiraModelObject] := Function[RetinotopyToCorticalMap[mdl, ##]];

Protect[SchiraModel, SchiraModelObject, SchiraFunction, SchiraInverse,
        CorticalMapToRetinotopy, RetinotopyToCorticalMap];

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
     f = mdl[RetinotopyToCorticalMap],
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


(* #SchiraLinePlot ********************************************************************************)
Options[SchiraLinePlot] = Join[
  FilterRules[Options[ParametricPlot], Except[PlotStyle | ColorFunction | ColorFunctionScaling]],
  {EccentricityStyleFunction -> Automatic,
   PolarAngleStyleFunction -> Automatic,
   EccentricityLines -> Automatic,
   PolarAngleLines -> Automatic,
   VisualAreas -> Automatic,
   Range -> Full}];
SchiraLinePlot[mdl_SchiraModelObject, opts : OptionsPattern[]] := Catch[
  With[
    {epsilon = 0.000001,
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
     fn = mdl[RetinotopyToCorticalMap],
     f = (SetAttributes[#, Temporary]; #)& @ Unique["f"],
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
    f[t_?NumericQ, r_?NumericQ, k_Integer] := Part[fn[t, r], k];
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

Protect[PolarAngleLegend, EccentricityLegend, SchiraParametricPlot, VisualAreas,
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
       preds = Transpose[RetinotopyToCorticalMap[model, angle[[idcs]], eccen[[idcs]]], {3, 1, 2}]},
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
             

End[];
EndPackage[];

