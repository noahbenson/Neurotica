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
BeginPackage["Neurotica`Registration`", {"Neurotica`Global`","Neurotica`Util`","Neurotica`Mesh`"}];
Unprotect["Neurotica`Registration`*", "Neurotica`Registration`Private`*"];
ClearAll[ "Neurotica`Registration`*", "Neurotica`Registration`Private`*"];

CorticalPotentialFunction::usage = "CorticalPotentialFunction[{F, G}, X] yields a cortical potential function object with potential F and gradient G in terms of the coordinate matrix X, which is assumed to be a 2 or 3 by n matrix when the potential is evaluated. The following options may be given:
  * Print (default: Subscript[\"F\", \"Potential\"]) specifies what symbol should be used to display the potential function.
  * MetaInformation (default: {}) specifies any optional meta-information to attach to the potential function.
  * CorticalMesh (default: None) specifies the (optional) cortical mesh for which this potential function was defined.";
PotentialFunction::usage = "PotentialFunction[f] yields a pure functional form of the cortical potential function instance, f.";
GradientFunction::usage = "GradientFunction[f] yields a pure functional form of the gradient of the cortical potential function instance, f.";
HessianFunction::usage = "HessianFunction[f] yields a pure functional form of the Hessian of the cortical potential function instance, f.";

CalculateAngleIntermediateData::usage = "CalculateAngleIntermediateData[a,b,c] yields the numerical array equivalent to Flatten[{u1,u2,n1,n2,{d1,d2,cos}}, 1] where u1, u2, n1, and n2 are the vectors b - a, c - a, Normalize[b - a], and Normalize[c - a], respectively, and where d1 and d2 are the lengths of the b-a and c-a, respectively, and where cos is the cosine of the angle centered at a. Note that a, b, and c are expected to be 2 or 3 by n arrays where each column of the array corresponds to a single point.";
CalculateAngle::usage = "CalculateAngle[a,b,c] yields the angle between the vectors b-a and c-a; the three arguments must be 2D arrays with length 2 or 3 (for 2d or 3d points) and the return value is a vector of angles, one for each column of a,b, and c."
CalculateAngleGradient::usage = "CalculateAngleGradient[a,b,c] yields the gradient of the angle between the vectors b-a and c-a in terms of the vertices of a, b, and c; the three arguments must be 2D arrays with length 2 or 3 (for 2d or 3d points) and the return value is a 3x3xn or 3x2xn (for 3D or 2D points, respectively) matrix of gradients, one for each column of a,b, and c."
CalculateAngleHessian::usage = "CalculateAngleHessian[a,b,c] yields the Hessian of the angle between the vectors b-a and c-a in terms of the vertices of a, b, and c; the three arguments must be 2D arrays with length 2 or 3 (for 2d or 3d points) and the return value is a 3x3x3x3xn or 3x2x3x2xn (for 3D or 2D points, respectively) matrix of Hessians, one for each column of a,b, and c."

HarmonicEdgePotential::usage = "HarmonicEdgePotential[mesh] yields a function symbol f such that f[X] is the harmonic edge potential where X is a possible vertex coordinate list for the given cortical mesh. The potential is calculated as the total of (d - d0)^2 where d is the distance between two vertices in X and d0 is the distance in the coordinates of the given mesh. Note that Grad[f, X] yields the numerical gradient of the potential at the vertex configuration given in X.
The potential function of a HarmonicEdgePotential is U[d] = 0.5 / n * (d - d0)^2 where d is the distance between a pair of vertices connected by an edge, n is the number of edges in the system, and d0 is the distance, in the initial mesh, between the two vertices.";

HarmonicAnglePotential::usage = "HarmonicAnglePotential[mesh] yields a function symbol f such that f[X] is the harmonic angle potential where X is a possible vertex coordinate list for the given cortical mesh. The potential is calculated as the total of (a - a0)^2 where a is the angle of a face in X and a0 is the angle of the same face corner in the original mesh. Note that Grad[f, X] yields the numerical gradient of the potential at the vertex configuration given in X.
The potential function of a HarmonicAnglePotential is U[a] = 0.5 / n * (a - a0)^2 where a is the angle of a corner of a face, n is the number of faces in the system, and a0 is the angle of the same face in the initial mesh.";

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
  * \[Gradient]f[x] = w (d[x]/\[Sigma])^(\[Beta] - 1) / \[Sigma] Exp[-(d[x]/\[Sigma])^\[Beta] / \[Beta]]";
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
  * \[Gradient]f[x] = w/k (d[x]/k)^(\[Beta] - 1)";
HarmonicPotentialWell::badarg = "Bad argument given to HarmonicPotentialWell: `1`";

RegionDistancePotential::usage = "RegionDistancePotential[mesh, reg, {F, G}] yields a cortical potential function with potential function F and gradient G. Both F and G must be functions such that F[dists] and G[dists], for a vector dists of the distances of the relevant vertices to the region, yield a total potential and a vector of gradient magnitudes, respectively. Note that the direction of the gradient is calculated automatically. The following additional options may be given:
  * Print (default: Style[\"\[GothicCapitalG]\", FontWeight -> Bold]) specifies the default display name for the function.
  * MetaInformation (default: {}) specifies extra meta information attached to the function.
  * VertexWeight (default: Automatic) specifies the relative strength of each vertex in the field; the region will have a repulsive, neutral, or attractive effect on any vertex with a weight less than, equal to, or greater than 0, respectively. If a property is named by this argument, then its values are used. The default value, Automatic, applies the field to all vertices.
See also SignedRegionDistancePotential.";
RegionDistancePotential::badarg = "Bad argument given to RegionDistancePotential: `1`";

SignedRegionDistancePotential::usage = "SignedRegionDistancePotential[mesh, reg, {F, G}] is identical to RegionDistancePotential[mesh, reg, {F, G}] except that the functions F and G are given signed distances to the BoundaryMeshRegion reg for the relevant vertices instead of the absolute distances.";
SignedRegionDistancePotential::badarg = "Bad argument given to SignedRegionDistancePotential: `1`";

HarmonicPerimeterPotential::usage = "HarmonicPerimeterPotential[map] yields a cortical potential function that operates on the vertices on the perimeter of the given map to hold them in place using a harmonic potential well tied to their initial positions.";

CortexGradientPlot::usage = "CortexGradientPlot[mesh, functions] yields a plot of the edges in the given mesh with the arrows representing the gradient of the vertices in the mesh, according to the list of potential functions given in the list functions. In addition to all options that are valid for CortexPlot, the following options may be given:
  * Arrowheads (default: Small) indicates that the arrowheads should be the given size in the plot.
  * PlotStyle (default: Automatic) should be a list of style instructions for the arrows of the gradients; these are cycled across the potential functions as in ListPlot.
  * Scaled (default: Automatic) indicates the absolute plotting length of the largest single gradient vector for any of the vertices; effectively, all gradients are scaled such that the largest gradient is equal to this value. If Automatic is given, then the value used is 75% of the mean edge length.";

CurvaturePotential::usage = "CurvaturePotential[mesh, template, sigma] yields a cortical potential function that enforces an aligment of the curvature in the cortical mesh, mesh, to the curvature in the cortical mesh, template, using a Gaussian smoothing function with shape parameter sigma over the cortical surface.";

MapTangledQ::usage = "MapTangledQ[map] yields True if and only if the given cortical map is tangled (has faces that are inverted); otherwise yields False.
MapTangledQ[map, X] is identical to MapTangledQ[map] except that it uses the coordinates given in X.";
MapTangles::usage = "MapTangles[map] yields a list of the vertices in the given cortical map that are tangles (their neighbors are not in the correct counter-clockwise ordering).
MapTangles[map, X] uses the coordinates given in X for the map.";
MapUntangle::usage = "MapUntangle[map] attempts to untangle the map given by repeatedly moving tangled vertices to their centroids (with respect to their neighbors); this may not succeed, but will return the coordinates regardless after 50 such attempts.
MapUntangle[map, X] uses the coordinates X as the map coordinates.
MapUntangle[map, X, max] attempts at most max times to untangle the map.";

Begin["`Private`"];

(* #CorticalPotentialFunction *********************************************************************)
Options[CorticalPotentialFunction] = {
  Print -> Subscript["F", "Potential"],
  MetaInformation -> {},
  CorticalMesh -> None};
DefineImmutable[
  CorticalPotentialFunction[{F_, G_, H_}, X_Symbol, OptionsPattern[]] :> P,
  {(* Let's start with options! *)
   Options[P] = Map[# -> OptionValue[#]&, Options[CorticalPotentialFunction][[All,1]]],
   Options[P, opt_] := Replace[
     opt,
     Append[Options[P], x_ :> (Message[Options::optionf, x, CorticalPotentialFunction]; {})]],
   (* We need to define the potential functions and gradients... *)
   PotentialFunction[P] -> Block[{X},
     With[
       {pos = Position[Hold[F], X, Infinity],
        sym = Unique["arg"]},
       Function @@ Join[
         Hold[{sym}],
         ReplacePart[Hold[F], (# -> sym)& /@ pos]]]],
   GradientFunction[P]  -> Block[
     {X},
     With[
       {pos = Position[Hold[G], X, Infinity],
        sym = Unique["arg"]},
       Function @@ Join[
         Hold[{sym}],
         ReplacePart[Hold[G], (# -> sym)& /@ pos]]]],
   HessianFunction[P]   -> Block[
     {X},
     With[
       {pos = Position[Hold[H], X, Infinity],
        sym = Unique["arg"]},
       Function @@ Join[
         Hold[{sym}],
         ReplacePart[Hold[H], (# -> sym)& /@ pos]]]],
   (* And a call form for the gradient. *)
   Grad[P, M_ /; ArrayQ[M, 2, NumericQ]] := With[
     {f = GradientFunction[P]},
     Join @@ If[Length[M] <= Length[M[[1]]], f[M], Transpose[f[Transpose @ M]]]],
   (* And a call form for the Hessian. *)
   Hessian[P, M_ /; ArrayQ[M, 2, NumericQ]] := With[
     {f = HessianFunction[P]},
     If[Length[M] <= Length[M[[1]]], f[M], Transpose[f[Transpose @ M]]]],
   (* Finally, this one is private, but useful for combining potential functions *)
   HeldArguments[P] -> {Hold[F], Hold[G], Hold[H], Hold[X]}},
  SetSafe -> True,
  Symbol -> CorticalPotentialFunctionInstance];
SetAttributes[CorticalPotentialFunction, HoldAll];
Protect[PotentialFunction, GradientFunction];

(* We have a few more edits for the potential function though: *)
Unprotect[CorticalPotentialFunctionInstance];

(* The call form for the potential is here, where it can be defined: *)
(CPF:CorticalPotentialFunctionInstance[__])[X_ /; ArrayQ[X, 2, NumericQ]] := With[
   {f = PotentialFunction[CPF]},
   If[Length[X] <= Length[X[[1]]], f[X], Transpose @ f[Transpose @ X]]];
(CPF:CorticalPotentialFunctionInstance[__])[X_Symbol] := With[
   {tmp = TemporarySymbol["CPF"],
    f = PotentialFunction[CPF]},
   tmp[M_ /; ArrayQ[M, 2, NumericQ]] := f[M];
   tmp[X]];

(* Critically, we want a nice display form for these potential functions: *)
MakeBoxes[P_CorticalPotentialFunctionInstance, form_] := MakeBoxes[#]& @ Options[P, Print];

(* We need to define some simple combination operators... *)
SetAttributes[AutoExtendCorticalPotential, HoldRest];
Quiet[
  AutoExtendCorticalPotential[
    patt_ -> res_,
    P0_Symbol -> P_Symbol, 
    {F0_Symbol, G0_Symbol, H0_Symbol} -> {F_, G_, H_},
    X_Symbol
   ] := TagSetDelayed[
      CorticalPotentialFunctionInstance,
      Evaluate[patt /. HoldPattern[P0] :> P0_CorticalPotentialFunctionInstance],
      With[
        {opts = Options[P0],
         args = HeldArguments[P0]},
        With[
          {fns = ReplaceAll[
             Hold[{F, G, H}],
             MapThread[
               RuleDelayed @@ Join[Hold@@{HoldPattern @@ #1}, #2] &,
               {{Hold[F0], Hold[G0], Hold[H0]},
                ReplaceAll[args[[1;;3]], HoldPattern @@ args[[4]] :> X]}]]},
          CorticalPotentialFunction @@ Join[
            fns,
            Hold[X],
            Hold @@ opts]]]],
  {RuleDelayed::rhs}];
Protect[AutoExtendCorticalPotential];

AutoExtendCorticalPotential[
  Plus[x_?NumericQ, P0, rest___] -> Plus[P, rest], P0 -> P,
  {F0, G0, H0} -> {x + F0, G0, H0},
  X];
AutoExtendCorticalPotential[
  Times[x_?NumericQ, P0, rest___] -> Times[P, rest], P0 -> P,
  {F0, G0, H0} -> {x * F0, x * G0, x * H0},
  X];
(* One weird one... *)
CorticalPotentialFunctionInstance /: Plus[P0_CorticalPotentialFunctionInstance,
                                          P1_CorticalPotentialFunctionInstance,
                                          rest___] := Plus[
  With[
    {args0 = HeldArguments[P0], args1 = HeldArguments[P1],
     opts0 = Options[P0], opts1 = Options[P1]},
    With[
      {rule = RuleDelayed @@ Join[Hold @@ {HoldPattern @@ args1[[4]]}, args0[[4]]]},
      CorticalPotentialFunction @@ Join[
        Replace[
          Hold @@ {
            {Join[args0[[1]], ReplaceAll[args1[[1]], rule]],
             Join[args0[[2]], ReplaceAll[args1[[2]], rule]],
             Join[args0[[3]], ReplaceAll[args1[[3]], rule]]}},
          Hold[a__] :> Plus[a],
          {2}],
        args0[[4]],
        Hold @@ {
          Print -> Row[{Print /. opts0, "+", Print /. opts1}],
          CorticalMesh -> With[
            {msh = CorticalMesh /. opts0},
            If[msh == (CorticalMesh /. opts1), msh, None]],
          MetaInformation -> ((MetaInformation /. #)& /@ {opts0, opts1})}]]],
  rest];
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
Protect[CalculateAngleIntermediateData];
CalculateAngle = Compile[
  {{a, _Real, 2}, {b, _Real, 2}, {c, _Real, 2}},
  With[
    {data = CalculateAngleIntermediateData[a,b,c]},
    Last[data]],
  {{CalculateAngleIntermediateData[_,_,_], _Real, 2}},
  RuntimeOptions -> {"Speed", "EvaluateSymbolically" -> False}, 
  Parallelization -> True];
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
Protect[CalculateHarmonicAnglePotential, CalculateHarmonicAngleGradient,
        CalculateHarmonicAngleHessian];

CalculateVDWAnglePotential = Compile[
  {{a, _Real, 2}, {b, _Real, 2}, {c, _Real, 2}, {th0, _Real, 1}},
  With[
    {data = CalculateAngleIntermediateData[a,b,c]},
    With[
      {ratio = 2.0 * Last[data] / th0},
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
      {th = Last[data]
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
  {With[
     {CheckIndices = Replace @ {
        id_Integer :> {id},
        ids:{_Integer..} :> ids,
        _ :> Message[HarmonicPotentialWell::badarg, "index "<>ToString[idcsArg]<>" not found"]}},
     Which[
       IntegerQ[idcsArg], CheckIndices @ VertexIndex[mesh, idcsArg],
       ArrayQ[idcsArg, 1, IntegerQ], CheckIndices @ VertexIndex[mesh, idcsArg],
       All, Range[VertexCount[mesh]],
       True, Message[
         HarmonicPotentialWell::badarg,
         "Harmonic well spec's indices must be an integer or a list of integers"]]],
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
      {0.5 / m * Total[(EdgeLengthsTr[mesh, X] - D0)^2],
       (* Gradient is easy: along the axis between the points *)
       With[
         {x1 = X[[All, E[[1]]]], 
          x2 = X[[All, E[[2]]]]},
         With[
           {dX = CalculateDistanceGradient[x1, x2],
            r = CalculateDistance[x1, x2]},
           SumOverEdgesDirectedTr[
             mesh,
             ConstantArray[r, dims] * dX[[1]] / m]]],
       (* This code produces the correct grad, but has been depricated:
        With[
          {dX = X[[All, E[[2]]]] - X[[All, E[[1]]]]},
          With[
            {norms = Sqrt[Total[dX^2]]},
            With[
              {magnitude = (D0 - norms) / (m * norms)},
              SumOverEdgesDirectedTr[
                mesh,
                dX * ConstantArray[magnitude, Length[dX]]]]]]*)
       (* Hessian is trickier *)
       With[
         {x1 = X[[All, E[[1]]]], 
          x2 = X[[All, E[[2]]]]},
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
      X,
      Print -> Subscript[Style["\[GothicCapitalH]",Bold], Row[{"Edges",",",Length@X0}]],
      CorticalMesh -> mesh,
      MetaInformation -> OptionValue[MetaInformation]]]];
Protect[HarmonicEdgePotential];

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
      {With[
         {corners0 = X[[All, #]]& /@ Ft},
         const * Sum[
           pfun @@ Append[RotateLeft[corners0, i], A0[[i+1]]],
           {i, 0, 2}]],
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
           const * Table[SumOverFaceVerticesTr[mesh, facesGrad[[All, k]]], {k, 1, dims}]]],
       (* No hessian yet *)
       None},
      X,
      Print -> Subscript[Style["\[GothicCapitalH]",Bold], Row[{"Angles",",",Length@X0}]],
      CorticalMesh -> mesh,
      MetaInformation -> OptionValue[MetaInformation]]]];
Protect[HarmonicAnglePotential];

(* #GaussianPotentialWell *************************************************************************)
Options[GaussianPotentialWell] = {MetaInformation -> {}};
GaussianPotentialWell[mesh_?CorticalObjectQ, spec_] := Check[
  With[
    {wells = Transpose @ ParseGaussianPotentialWells[mesh, spec]},
    CorticalPotentialFunction[
      {CalculateGaussianPotential[wells, X],
       CalculateGaussianGradient[wells, X],
       CalculateGaussianHessian[wells, X]},
      X,
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
    {CalculateHarmonicPotential[wells, X], CalculateHarmonicGradient[wells, X], None},
    X,
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
    {idcs = Indices[weight, Except[0|0.0]],
     distFn = RegionDistance[reg],
     nearFn = RegionNearest[reg],
     W = Total[Abs@weight]},
    With[
      {w = weight[[idcs]], absw = Abs[weight[[idcs]]]},
      CorticalPotentialFunction[
        {Dot[F @ distFn @ Transpose @ X[[All, idcs]], absw] / W,
         With[
           {nears = nearFn @ Transpose @ X},
           With[
             {dX = X - Transpose[nears]},
             With[
               {dists = ColumnNorms[dX]},
               dX * ConstantArray[
                 w * G[dists] / (W * (# + (1 - Unitize[#]))&[Chop @ dists]),
                 Length[X]]]]],
         None},
        X,
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
     nearFn = RegionNearest[reg],
     W = Total[Abs@weight]},
    With[
      {w = weight[[idcs]], absw = Abs[weight[[idcs]]]},
      CorticalPotentialFunction[
        {Dot[F @ distFn @ Transpose @ X[[All, idcs]], absw] / W,
         With[
           {nears = nearFn @ Transpose @ X},
           With[
             {dX = X - Transpose[nears]},
             With[
               {dists = distFn[Transpose @ X[[All, idcs]]]},
               dX * ConstantArray[w * G[dists] / (W * dists), Length[X]]]]],
         None},
        X,
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
MapTangledQ[map_?CorticalMapQ, X_] := With[
  {nei = VertexIndex[map, NeighborhoodList[map]]},
  Catch[
    MapThread[
      Function[{x, n},
        If[Length[n] > 2,
          With[
            {xnei = X[[#]] & /@ n},
            With[
              {ord = Ordering[ArcTan[#[[1]] - x[[1]], #[[2]] - x[[2]]] & /@ xnei]},
              With[
                {rot = RotateLeft[ord, Position[ord, 1, {1}][[1, 1]] - 1]},
                If[Range[Length[rot]] != rot || !Graphics`Mesh`InPolygonQ[xnei, x],
                  Throw[True]]]]]]],
      {X, nei}];
    False]];
MapTangledQ[map_?CorticalMapQ] := MapTangledQ[map, VertexCoordinates[map]];
Protect[MapTangledQ];

(* #MapTangles ************************************************************************************)
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
(*
With[
  {P = HarmonicEdgePotential[map] + HarmonicAnglePotential[map],
   hold = Replace[OptionValue[Hold], None->{}],
   nei = VertexIndex[map, NeighborhoodList[map]]},
  NestWhile[
    Function @ With[
      {X0 = #,
       X0t = Transpose[#]},
      With[
        {tangles = With[
           {t0 = MapTangles[map, X0]},
           Complement[Union[Flatten[{t0, nei[[#]] & /@ t0}]], hold]]},
        With[
          {gradpart = Join[tangles, tangles + Length[X0]]},
          Block[
            {X, f, g},
            f[x_List] := P[MapThread[ReplacePart[#1, Thread[tangles -> #2]]&, {X0t, x}]];
            g[x_List] := Part[
              Grad[P, MapThread[ReplacePart[#1, Thread[tangles -> #2]]&, {X0t, x}]],
              gradpart];
            ReplacePart[
              X0,
              Thread[
                tangles -> Quiet[
                  First@FindArgMin[
                    f[X],
                    {X, X0t[[All, tangles]]},
                    Gradient :> g[X],
                    Method -> {"QuasiNewton",
                               "StepControl" -> {"LineSearch", "CurvatureFactor" -> 1.0}},
                    AccuracyGoal -> 6,
                    MaxIterations -> 100],
                  {FindArgMin::cvmit, FindArgMin::lstol}]]]]]]],
    Xtangled,
    MapTangledQ[map, #] &,
    1,
    max]];
*)


End[];
EndPackage[];

