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
  * Name (default: Subscript[\"F\", \"Potential\"]) specifies what symbol should be used to display the potential function.";
PotentialFunction::usage = "PotentialFunction[f] yields a pure functional form of the cortical potential function instance, f.";
GradientFunction::usage = "GradientFunction[f] yields a pure functional form of the gradient of the cortical potential function instance, f.";

HarmonicEdgePotential::usage = "HarmonicEdgePotential[mesh] yields a function symbol f such that f[X] is the harmonic edge potential where X is a possible vertex coordinate list for the given cortical mesh. The potential is calculated as the total of (d - d0)^2 where d is the distance between two vertices in X and d0 is the distance in the coordinates of the given mesh. Note that Grad[f, X] yields the numerical gradient of the potential at the vertex configuration given in X.
The potential function of a HarmonicEdgePotential is U[d] = 0.5 / n * (d - d0)^2 where d is the distance between a pair of vertices connected by an edge, n is the number of edges in the system, and d0 is the distance, in the initial mesh, between the two vertices.";

HarmonicAnglePotential::usage = "HarmonicAnglePotential[mesh] yields a function symbol f such that f[X] is the harmonic angle potential where X is a possible vertex coordinate list for the given cortical mesh. The potential is calculated as the total of (a - a0)^2 where a is the angle of a face in X and a0 is the angle of the same face corner in the original mesh. Note that Grad[f, X] yields the numerical gradient of the potential at the vertex configuration given in X.
The potential function of a HarmonicAnglePotential is U[a] = 0.5 / n * (a - a0)^2 where a is the angle of a corner of a face, n is the number of faces in the system, and a0 is the angle of the same face in the initial mesh.";

CosineAnglePotential::usage = "CosineAnglePotential[mesh] yields a function symbol f such that f[X] is the cosine angle potential where X is a possible vertex coordinate list for the given cortical mesh. The potential is calculated as the total of (Cos[a] - Cos[a0])^2 where a is the angle of a face in X and a0 is the angle of the same face corner in the original mesh. Note that Grad[f, X] yields the numerical gradient of the potential at the vertex configuration given in X.
The potential function of a HarmonicAnglePotential is U[a] = 0.5 / n * (Cos[a] - Cos[a0])^2 where a is the angle of a corner of a face, n is the number of faces in the system, and a0 is the angle of the same face in the initial mesh.";

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

MapToMeshPotential::usage = "MapToMeshPotential[map, f] yields a function equivalent to f projected onto the cortical mesh origin of the given map. A heuristic is used when necessary. Note that this may be slow for many maps, so the use of small orthographic maps is suggested.";

Begin["`Private`"];

(* #CorticalPotentialFunction *********************************************************************)
Options[CorticalPotentialFunction] = {
  Print -> Subscript["F", "Potential"],
  MetaInformation -> {},
  CorticalMesh -> None};
DefineImmutable[
  CorticalPotentialFunction[{F_, G_}, X_Symbol, OptionsPattern[]] :> P,
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
   (* And a call form for the gradient. *)
   Grad[P, M_ /; ArrayQ[M, 2, NumericQ]] := With[
     {f = GradientFunction[P]},
     Join @@ If[Length[M] <= Length[M[[1]]], f[M], Transpose[f[Transpose @ M]]]],
   (* Finally, this one is private, but useful for combining potential functions *)
   HeldArguments[P] -> {Hold[F], Hold[G], Hold[X]}},
  SetSafe -> True,
  Symbol -> CorticalPotentialFunctionInstance];
SetAttributes[CorticalPotentialFunction, HoldAll];
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

(* We need to define some simple combination operators... *)
SetAttributes[AutoExtendCorticalPotential, HoldRest];
Quiet[
  AutoExtendCorticalPotential[
    patt_ -> res_,
    P0_Symbol -> P_Symbol, 
    {F0_Symbol, G0_Symbol} -> {F_, G_}, X_Symbol
    ] := TagSetDelayed[
      CorticalPotentialFunctionInstance,
      Evaluate[patt /. HoldPattern[P0] :> P0_CorticalPotentialFunctionInstance],
      With[
        {opts = Options[P0],
         args = HeldArguments[P0]},
        With[
          {fns = ReplaceAll[
             Hold[{F, G}],
             MapThread[
               RuleDelayed @@ Join[Hold@@{HoldPattern @@ #1}, #2] &,
               {{Hold[F0], Hold[G0]},
                ReplaceAll[args[[1;;2]], HoldPattern @@ args[[3]] :> X]}]]},
          CorticalPotentialFunction @@ Join[
            fns,
            Hold[X],
            Hold @@ opts]]]],
  {RuleDelayed::rhs}];
Protect[AutoExtendCorticalPotential];

AutoExtendCorticalPotential[
  Plus[x_?NumericQ, P0, rest___] -> Plus[P, rest], P0 -> P,
  {F0, G0} -> {x + F0, G0},
  X];
AutoExtendCorticalPotential[
  Times[x_?NumericQ, P0, rest___] -> Times[P, rest], P0 -> P,
  {F0, G0} -> {x * F0, x * G0},
  X];
(* One weird one... *)
CorticalPotentialFunctionInstance /: Plus[P0_CorticalPotentialFunctionInstance,
                                          P1_CorticalPotentialFunctionInstance,
                                          rest___] := Plus[
  With[
    {args0 = HeldArguments[P0], args1 = HeldArguments[P1],
     opts0 = Options[P0], opts1 = Options[P1]},
    With[
      {rule = RuleDelayed @@ Join[Hold @@ {HoldPattern @@ args1[[3]]}, args0[[3]]]},
      CorticalPotentialFunction @@ Join[
        Replace[
          Hold @@ {
            {Join[args0[[1]], ReplaceAll[args1[[1]], rule]],
             Join[args0[[2]], ReplaceAll[args1[[2]], rule]]}},
          Hold[a__] :> Plus[a],
          {2}],
        args0[[3]],
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

(* Here we compile functions for calculating harmonic angle potentials *)
CalculateHarmonicAnglePotential3D = Compile[
  {{x0, _Real, 1}, {y0, _Real, 1}, {z0, _Real, 1},
   {x1, _Real, 1}, {y1, _Real, 1}, {z1, _Real, 1},
   {x2, _Real, 1}, {y2, _Real, 1}, {z2, _Real, 1},
   {th0, _Real, 1}},
  With[
    {u01 = {x1 - x0, y1 - y0, z1 - z0},
     u02 = {x2 - x0, y2 - y0, z2 - z0}},
    With[
      {l01 = Total[u01^2],
       l02 = Total[u02^2]},
      Total[(ArcCos[Total[u01 * u02] / Sqrt[l01 * l02]] - th0)^2]]],
  RuntimeOptions -> {"Speed", "EvaluateSymbolically" -> False},
  Parallelization -> True];
CalculateHarmonicAngleGradient3D = Compile[
    {{x0, _Real, 1}, {y0, _Real, 1}, {z0, _Real, 1},
     {x1, _Real, 1}, {y1, _Real, 1}, {z1, _Real, 1},
     {x2, _Real, 1}, {y2, _Real, 1}, {z2, _Real, 1},
     {th0, _Real, 1}},
    With[
      {u1 = {x1 - x0, y1 - y0, z1 - z0}, u2 = {x2 - x0, y2 - y0, z2 - z0}},
      With[
        {d1 = Sqrt[Total[u1^2]], d2 = Sqrt[Total[u2^2]]},
        With[
          {th = ArcCos[Total[u1*u2] / (d1*d2)],
           n1 = u1 / ConstantArray[d1, Length[u1]],
           n2 = u2 / ConstantArray[d2, Length[u2]]},
          With[
            {cos = ConstantArray[Cos[th], Length[u1]], sin = Sin[th]},
            With[
              {f1 = (cos*n1 - n2) * ConstantArray[(th - th0) / (d1 * sin), Length[u1]],
               f2 = (cos*n2 - n1) * ConstantArray[(th - th0) / (d2 * sin), Length[u1]]},
              {-(f1 + f2), f1, f2}]]]]],
    RuntimeOptions -> {"Speed", "EvaluateSymbolically" -> False},
    Parallelization -> True];
Protect[CalculateHarmonicAnglePotential3D, CalculateHarmonicAngleGradient3D];

CalculateHarmonicAnglePotential2D = Compile[
  {{x0, _Real, 1}, {y0, _Real, 1},
   {x1, _Real, 1}, {y1, _Real, 1},
   {x2, _Real, 1}, {y2, _Real, 1},
   {th0, _Real, 1}},
  With[
    {u01 = {x1 - x0, y1 - y0},
     u02 = {x2 - x0, y2 - y0}},
    With[
      {l01 = Total[u01^2],
       l02 = Total[u02^2]},
      Total[(ArcCos[Total[u01 * u02] / Sqrt[l01 * l02]] - th0)^2]]],
  RuntimeOptions -> {"Speed", "EvaluateSymbolically" -> False},
  Parallelization -> True];
CalculateHarmonicAngleGradient2D = Compile[
    {{x0, _Real, 1}, {y0, _Real, 1},
     {x1, _Real, 1}, {y1, _Real, 1},
     {x2, _Real, 1}, {y2, _Real, 1},
     {th0, _Real, 1}},
    With[
      {u1 = {x1 - x0, y1 - y0}, u2 = {x2 - x0, y2 - y0}},
      With[
        {d1 = Sqrt[Total[u1^2]], d2 = Sqrt[Total[u2^2]]},
        With[
          {th = ArcCos[Total[u1*u2] / (d1*d2)],
           n1 = u1 / ConstantArray[d1, Length[u1]],
           n2 = u2 / ConstantArray[d2, Length[u2]]},
          With[
            {cos = ConstantArray[Cos[th], Length[u1]], sin = Sin[th]},
            With[
              {f1 = (cos*n1 - n2) * ConstantArray[(th - th0) / (d1 * sin), Length[u1]],
               f2 = (cos*n2 - n1) * ConstantArray[(th - th0) / (d2 * sin), Length[u1]]},
              {-(f1 + f2), f1, f2}]]]]],
    RuntimeOptions -> {"Speed", "EvaluateSymbolically" -> False},
    Parallelization -> True];
(*
CalculateHarmonicAngleGradient2D = ReplacePart[
  Hold[
    {{x0, _Real, 1}, {y0, _Real, 1},
     {x1, _Real, 1}, {y1, _Real, 1},
     {x2, _Real, 1}, {y2, _Real, 1},
     {th0, _Real, 1}},
    Evaluate @ Block[
      {x0, y0, x1, y1, x2, y2, th0},
      With[
        {grad = Grad[
           Simplify[
             (ArcCos[
                Dot[
                  Normalize[{x1 - x0, y1 - y0}],
                  Normalize[{x2 - x0, y2 - y0}]]] - th0)^2,
             Assumptions -> Element[{x0, y0, x1, y1, x2, y2}, Reals]],
           {x0, y0, x1, y1, x2, y2}]},
        Hold[grad, 2]]],
    RuntimeOptions -> {"Speed", "EvaluateSymbolically" -> False},
    Parallelization -> True],
  {{2,0} -> Partition,
   {0}   -> Compile}];
*)
Protect[CalculateHarmonicAnglePotential2D, CalculateHarmonicAngleGradient2D];

(* And here we compile functions for calculating cosine angle potentials *)
CalculateCosineAnglePotential3D = Compile[
  {{x0, _Real, 1}, {y0, _Real, 1}, {z0, _Real, 1},
   {x1, _Real, 1}, {y1, _Real, 1}, {z1, _Real, 1},
   {x2, _Real, 1}, {y2, _Real, 1}, {z2, _Real, 1},
   {cos0, _Real, 1}},
  With[
    {u01 = {x1 - x0, y1 - y0, z1 - z0},
     u02 = {x2 - x0, y2 - y0, z2 - z0}},
    With[
      {l01 = Total[u01^2],
       l02 = Total[u02^2]},
      Total[(Total[u01 * u02] / Sqrt[l01 * l02] - cos0)^2]]],
  RuntimeOptions -> {"Speed", "EvaluateSymbolically" -> False},
  Parallelization -> True];
CalculateCosineAngleGradient3D = ReplacePart[
  Hold[
    {{x0, _Real, 1}, {y0, _Real, 1}, {z0, _Real, 1},
     {x1, _Real, 1}, {y1, _Real, 1}, {z1, _Real, 1},
     {x2, _Real, 1}, {y2, _Real, 1}, {z2, _Real, 1},
     {cos0, _Real, 1}},
    Evaluate @ Block[
      {x0, y0, z0, x1, y1, z1, x2, y2, z2, cos0},
      With[
        {grad = Grad[
           Simplify[
             (Dot[
                Normalize[{x1 - x0, y1 - y0, z1 - z0}],
                Normalize[{x2 - x0, y2 - y0, z2 - z0}]] - cos0)^2,
             Assumptions -> Element[{x0, y0, z0, x1, y1, z1, x2, y2, z2}, Reals]],
           {x0, y0, z0, x1, y1, z1, x2, y2, z2}]},
        Hold[grad, 2]]],
    RuntimeOptions -> {"Speed", "EvaluateSymbolically" -> False},
    Parallelization -> True],
  {{2,0} -> Partition,
   {0}   -> Compile}];
Protect[CalculateCosineAnglePotential3D, CalculateCosineAngleGradient3D];

CalculateCosineAnglePotential2D = Compile[
  {{x0, _Real, 1}, {y0, _Real, 1},
   {x1, _Real, 1}, {y1, _Real, 1},
   {x2, _Real, 1}, {y2, _Real, 1},
   {cos0, _Real, 1}},
  With[
    {u01 = {x1 - x0, y1 - y0},
     u02 = {x2 - x0, y2 - y0}},
    With[
      {l01 = Total[u01^2],
       l02 = Total[u02^2]},
      Total[(Total[u01 * u02] / Sqrt[l01 * l02] - cos0)^2]]],
  RuntimeOptions -> {"Speed", "EvaluateSymbolically" -> False},
  Parallelization -> True];
CalculateCosineAngleGradient2D = ReplacePart[
  Hold[
    {{x0, _Real, 1}, {y0, _Real, 1},
     {x1, _Real, 1}, {y1, _Real, 1},
     {x2, _Real, 1}, {y2, _Real, 1},
     {cos0, _Real, 1}},
    Evaluate @ Block[
      {x0, y0, x1, y1, x2, y2, cos0},
      With[
        {grad = Grad[
           Simplify[
             (Dot[
                Normalize[{x1 - x0, y1 - y0}],
                Normalize[{x2 - x0, y2 - y0}]] - cos0)^2,
             Assumptions -> Element[{x0, y0, x1, y1, x2, y2}, Reals]],
           {x0, y0, x1, y1, x2, y2}]},
        Hold[grad, 2]]],
    RuntimeOptions -> {"Speed", "EvaluateSymbolically" -> False},
    Parallelization -> True],
  {{2,0} -> Partition,
   {0}   -> Compile}];
Protect[CalculateCosineAnglePotential2D, CalculateCosineAngleGradient2D];

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
   weight = OptionValue["Weight"]},
  With[
    {sig = If[TrueQ[fwhm], 2.0*Sqrt[2.0*Log[2.0]] * sigma, sigma]},
    Block[
      {X, U, distance},
      {ReplacePart[ (* Potential function *)
         Hold[
           {{X, _Real, 2}},
           Evaluate @ Hold[
             {distance = Sqrt @ Total[MapThread[Subtract, {X, x0}]^2]},
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
             {U = MapThread[Subtract, {x0, X}]},
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
         {{2,2,2,0} -> Function[U * ConstantArray[#, Length[U]]],
          {2,2,0}   -> With,
          {2,0}     -> With,
          {0}       -> Compile}]}]]];
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

(* These functions are helpers for the Harmonic potential well functions *)
Options[HarmonicPotentialWellParseHarmonic] = {
  "Width" -> 1,
  "Shape" -> 2,
  "Weight" -> 1};
HarmonicPotentialWellParseHarmonic[x0_, OptionsPattern[]] := With[
  {width = OptionValue["Width"],
   shape = OptionValue["Shape"],
   weight = OptionValue["Weight"]},
  Block[
    {X, U, distance},
    {(* Potential function *)
     Compile[{{X, _Real, 2}},
       weight / shape * (Sqrt[Total[MapThread[Subtract, {X, x0}]^2]] / width)^shape,
       RuntimeOptions -> {"Speed", "EvaluateSymbolically" -> False},
       Parallelization -> True],
     (* Gradient function *)
     Compile[
       {{X, _Real, 2}},
       With[
         {U = MapThread[Subtract, {X, x0}]},
         With[
           {distances = Sqrt @ Total[U^2]},
           With[
             {magnitudes = weight * (distances / width)^(shape) / distances^2},
             U * ConstantArray[magnitudes, Length[U]]]]],
       RuntimeOptions -> {"Speed", "EvaluateSymbolically" -> False},
       Parallelization -> True]}]];
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
   m = EdgeCount[mesh]},
  CorticalPotentialFunction[
    {0.5 / m * Total[(EdgeLengthsTr[mesh, X] - D0)^2],
     With[
       {dX = X[[All, E[[2]]]] - X[[All, E[[1]]]]},
       With[
         {norms = Sqrt[Total[dX^2]]},
         With[
           {magnitude = (D0 - norms) / (m * norms)},
           SumOverEdgesDirectedTr[
             mesh,
             dX * ConstantArray[magnitude, Length[dX]]]]]]},
    X,
    Print -> Subscript["\[GothicCapitalH]", Row[{"Edges",",",Length@X0}]],
    CorticalMesh -> mesh,
    MetaInformation -> OptionValue[MetaInformation]]];
Protect[HarmonicEdgePotential];

(* #HarmonicAnglePotential ************************************************************************)
Options[HarmonicAnglePotential] = {MetaInformation -> {}};
HarmonicAnglePotential[mesh_?CorticalObjectQ, OptionsPattern[]] := With[
  {X0 = VertexCoordinatesTr[mesh],
   Ft = VertexIndex[mesh, FaceListTr[mesh]],
   A0 = FaceAnglesTr[mesh],
   n = 3 * FaceCount[mesh],
   dims = If[CorticalMeshQ[mesh], 3, 2],
   pfun = If[CorticalMeshQ[mesh],
     CalculateHarmonicAnglePotential3D,
     CalculateHarmonicAnglePotential2D],
   gfun = If[CorticalMeshQ[mesh],
     CalculateHarmonicAngleGradient3D,
     CalculateHarmonicAngleGradient2D]},
  CorticalPotentialFunction[
    {With[
       {corners0 = X[[All, #]]& /@ Ft},
       0.5 / n * Sum[
         pfun @@ Append[Join @@ RotateLeft[corners0, i], A0[[i+1]]],
         {i, 0, 2}]],
     With[
       {corners0 = X[[All, #]]& /@ Ft},
       (* corners0: {v1, v2, v3} (3xdxn); vi: {x, y, z} (3 x n) or {x, y} (2 x n) *)
       With[
         {facesGrad = Sum[
            RotateRight[
              gfun @@ Append[Join @@ RotateLeft[corners0, i], A0[[i + 1]]],
              i],
            {i, 0, 2}]},
         (* facesGrad: same format as corners0 *)
         Table[SumOverFaceVerticesTr[mesh, facesGrad[[All, k]]], {k, 1, dims}] / n]]},
    X,
    Print -> Subscript["\[GothicCapitalH]", Row[{"Angles",",",Length@X0}]],
    CorticalMesh -> mesh,
    MetaInformation -> OptionValue[MetaInformation]]];
Protect[HarmonicAnglePotential];

(* #CosineAnglePotential **************************************************************************)
Options[CosineAnglePotential] = {MetaInformation -> {}};
CosineAnglePotential[mesh_?CorticalObjectQ, OptionsPattern[]] := With[
  {X0 = VertexCoordinatesTr[mesh],
   Ft = VertexIndex[mesh, FaceListTr[mesh]],
   A0 = FaceAngleCosinesTr[mesh],
   n = 3 * FaceCount[mesh],
   dims = If[CorticalMeshQ[mesh], 3, 2],
   pfun = If[CorticalMeshQ[mesh],
     CalculateCosineAnglePotential3D,
     CalculateCosineAnglePotential2D],
   gfun = If[CorticalMeshQ[mesh],
     CalculateCosineAngleGradient3D,
     CalculateCosineAngleGradient2D]},
  CorticalPotentialFunction[
    {With[
       {corners0 = X[[All, #]]& /@ Ft},
       0.5 / n * Sum[
         pfun @@ Append[Join @@ RotateLeft[corners0, i], A0[[i+1]]],
         {i, 0, 2}]],
     With[
       {corners0 = X[[All, #]]& /@ Ft},
       (* corners0: {v1, v2, v3} (3x3xn); vi: {x, y, z} (3 x n) *)
       With[
         {facesGrad = Sum[
            RotateRight[
              gfun @@ Append[Join @@ RotateLeft[corners0, i], A0[[i+1]]],
              1],
            {i, 0, 2}]},
         (* facesGrad: same format as corners0 *)
         Table[SumOverFaceVerticesTr[mesh, facesGrad[[All, k]]], {k, 1, dims}] / n]]},
    X,
    Print -> Subscript["\[GothicCapitalH]", Row[{Subscript["Angles", "Cosine"],",",Length@X0}]],
    CorticalMesh -> mesh,
    MetaInformation -> OptionValue[MetaInformation]]];
Protect[CosineAnglePotential];

(* #GaussianPotentialWell *************************************************************************)
Options[GaussianPotentialWell] = {MetaInformation -> {}};
GaussianPotentialWell[mesh_?CorticalObjectQ, spec_] := Check[
  With[
    {wells = Transpose @ ParseGaussianPotentialWells[mesh, spec]},
    CorticalPotentialFunction[
      {CalculateGaussianGradient[wells, X], CalculateGaussianPotential[wells, X]},
      X,
      Print -> Subscript["\[GothicCapitalH]", Row[{Subscript["Angles", "Cosine"],",",Length@X0}]],
      CorticalMesh -> mesh,
      MetaInformation -> OptionValue[MetaInformation]]],
  $Failed];
Protect[GaussianPotentialWell];

(* #HarmonicPotentialWell *************************************************************************)
Options[HarmonicPotentialWell] = {MetaInformation -> {}};
HarmonicPotentialWell[mesh_?CorticalObjectQ, spec_] := With[
  {wells = Transpose @ ParseHarmonicPotentialWells[mesh, spec],
   dims = If[CorticalMeshQ[mesh], 3, 2]},
  CorticalPotentialFunction[
    {CalculateHarmonicPotential[wells, X], CalculateHarmonicGradient[wells, X]},
    X,
    Print -> Subscript[Style["\[GothicCapitalH]",Bold], Row[{Length[wells],",",Length[dims]}]],
    CorticalMesh -> mesh,
    MetaInformation -> OptionValue[MetaInformation]]];
Protect[HarmonicPotentialWell];

(* #MapToMeshPotential ****************************************************************************)
(* #here *)
MapToMeshPotential[map_?CoticalMapQ, f_] := None;

End[];
EndPackage[];

