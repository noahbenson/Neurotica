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

CorticalPotential::usage = "CorticalPotential[expression] ";

HarmonicEdgePotential::usage = "HarmonicEdgePotential[mesh] yields a function symbol f such that f[X] is the harmonic edge potential where X is a possible vertex coordinate list for the given cortical mesh. The potential is calculated as the total of (d - d0)^2 where d is the distance between two vertices in X and d0 is the distance in the coordinates of the given mesh. Note that Grad[f, X] yields the numerical gradient of the potential at the vertex configuration given in X.
The potential function of a HarmonicEdgePotential is U[d] = 0.5 / n * (d - d0)^2 where d is the distance between a pair of vertices connected by an edge, n is the number of edges in the system, and d0 is the distance, in the initial mesh, between the two vertices.";

HarmonicAnglePotential::usage = "HarmonicAnglePotential[mesh] yields a function symbol f such that f[X] is the harmonic angle potential where X is a possible vertex coordinate list for the given cortical mesh. The potential is calculated as the total of (a - a0)^2 where a is the angle of a face in X and a0 is the angle of the same face corner in the original mesh. Note that Grad[f, X] yields the numerical gradient of the potential at the vertex configuration given in X.
The potential function of a HarmonicAnglePotential is U[a] = 0.5 / n * (a - a0)^2 where a is the angle of a corner of a face, n is the number of faces in the system, and a0 is the angle of the same face in the initial mesh.";

CosineAnglePotential::usage = "CosineAnglePotential[mesh] yields a function symbol f such that f[X] is the cosine angle potential where X is a possible vertex coordinate list for the given cortical mesh. The potential is calculated as the total of (Cos[a] - Cos[a0])^2 where a is the angle of a face in X and a0 is the angle of the same face corner in the original mesh. Note that Grad[f, X] yields the numerical gradient of the potential at the vertex configuration given in X.
The potential function of a HarmonicAnglePotential is U[a] = 0.5 / n * (Cos[a] - Cos[a0])^2 where a is the angle of a corner of a face, n is the number of faces in the system, and a0 is the angle of the same face in the initial mesh.";

Begin["`Private`"];

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
      (ArcCos[Total[u01 * u02] / Sqrt[l01 * l02]] - th0)^2]],
  RuntimeOptions -> {"Speed", "EvaluateSymbolically" -> False},
  Parallelization -> True];
CalculateHarmonicAngleGradient3D = ReplacePart[
  Hold[
    {{x0, _Real, 1}, {y0, _Real, 1}, {z0, _Real, 1},
     {x1, _Real, 1}, {y1, _Real, 1}, {z1, _Real, 1},
     {x2, _Real, 1}, {y2, _Real, 1}, {z2, _Real, 1},
     {th0, _Real, 1}},
    Evaluate @ Block[
      {x0, y0, z0, x1, y1, z1, x2, y2, z2, th0},
      With[
        {grad = Grad[
           Simplify[
             (ArcCos[
                Dot[
                  Normalize[{x1 - x0, y1 - y0, z1 - z0}],
                  Normalize[{x2 - x0, y2 - y0, z2 - z0}]]] - th0)^2,
             Assumptions -> Element[{x0, y0, z0, x1, y1, z1, x2, y2, z2}, Reals]],
           {x0, y0, z0, x1, y1, z1, x2, y2, z2}]},
        Hold[grad, 3]]],
    RuntimeOptions -> {"Speed", "EvaluateSymbolically" -> False},
    Parallelization -> True],
  {{2,0} -> Partition,
   {0}   -> Compile}];
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
      (ArcCos[Total[u01 * u02] / Sqrt[l01 * l02]] - th0)^2]],
  RuntimeOptions -> {"Speed", "EvaluateSymbolically" -> False},
  Parallelization -> True];
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
        Hold[grad, 3]]],
    RuntimeOptions -> {"Speed", "EvaluateSymbolically" -> False},
    Parallelization -> True],
  {{2,0} -> Partition,
   {0}   -> Compile}];
Protect[CalculateHarmonicAnglePotential2D, CalculateHarmonicAngleGradient2D];



(* ============================================================================================== *)
(* ====================================== Public Functions ====================================== *)
(* ============================================================================================== *)

(* #HarmonicEdgePotential *************************************************************************)
HarmonicEdgePotential[mesh_?CorticalMeshQ] := With[
  {X0 = VertexCoordinatesTr[mesh],
   D0 = EdgeLengths[mesh],
   E = EdgePairsTr[mesh],
   m = EdgeCount[mesh],
   f = TemporarySymbol["edgePotential"]},
  With[
    {df = Function @ With[
       {deltaX = #[[E[[2]]]] - #[[E[[1]]]]},
       With[
         {norms = Sqrt[Total[deltaX^2]]},
         With[
           {magnitude = (norms - D0) / m},
           SumOverEdgesDirectedTr[
             dX / {norms, norms, norms} * {magnitude, magnitude, magnitude}]]]]},
    f /: Grad[f, Xarg_List] := Apply[
      Join,
      If[Length[Xarg] == 3, df[Xarg], Transpose @ df[Transpose @ Xarg]]];
    f[Xarg_List] := (0.5 / m * Total[(# - D0)^2])&[
      If[Length@Xarg == 3, EdgeLengthsTr[mesh, Xarg], EdgeLengths[mesh, Xarg]]];
    f]];
HarmonicEdgePotential[mesh_?CorticalMapQ] := With[
  {X0 = VertexCoordinatesTr[mesh],
   D0 = EdgeLengths[mesh],
   E = EdgePairsTr[mesh],
   m = EdgeCount[mesh],
   f = TemporarySymbol["edgePotential"]},
  With[
    {df = Function @ With[
       {deltaX = #[[E[[2]]]] - #[[E[[1]]]]},
       With[
         {norms = Sqrt[Total[deltaX^2]]},
         With[
           {magnitude = (norms - D0) / m},
           SumOverEdgesDirectedTr[dX / {norms, norms} * {magnitude, magnitude}]]]]},
    f /: Grad[f, Xarg_List] := Apply[
      Join,
      If[Length[Xarg] == 2, df[Xarg], Transpose @ df[Transpose @ Xarg]]];
    f[Xarg_List] := (0.5 / m * Total[(# - D0)^2])&[
      If[Length@Xarg == 2, EdgeLengthsTr[mesh, Xarg], EdgeLengths[mesh, Xarg]]];
    f]];
Protect[HarmonicEdgePotential];

(* #HarmonicEdgePotential *************************************************************************)
HarmonicAnglePotential[mesh_?CorticalMeshQ] := With[
  {X0 = VertexCoordinatesTr[mesh],
   Ft = FaceListTr[mesh],
   A0 = FaceAnglesTr[mesh],
   n = 3 * FaceCount[mesh],
   f = TemporarySymbol["anglePotential"]},
  f /: Grad[f, Xarg_List] := With[
    {Xt = If[Length[Xarg] == 3, Xarg, Transpose @ Xarg]},
    Join @@ With[
      {corners0 = Xt[[All, #]]& /@ Ft},
      With[
        {grad = Sum[
          Apply[
            CalculateHarmonicAngleGradient3D,
            Append[Join @@ RotateLeft[corners0, i], A0[[i+1]]]],
          {i, 0, 2}]},
        Join @@ If[Length[Xarg] == 3, Transpose @ grad, grad]]]];
  f[Xarg_List] := With[
    {Xt = If[Length[Xarg] == 3, Xarg, Transpose @ Xarg]},
    With[
      {corners0 = Xt[[All, #]]& /@ Ft},
      0.5 / m * Sum[
        Apply[
          CalculateHarmonicAnglePotential3D,
          Append[Join @@ RotateLeft[corners0, i], A0[[i+1]]]],
        {i, 0, 2}]]];
  f];
HarmonicAnglePotential[mesh_?CorticalMapQ] := With[
  {X0 = VertexCoordinatesTr[mesh],
   Ft = FaceListTr[mesh],
   A0 = FaceAnglesTr[mesh],
   n = 3 * FaceCount[mesh],
   f = TemporarySymbol["anglePotential"]},
  f /: Grad[f, Xarg_List] := With[
    {Xt = If[Length[Xarg] == 2, Xarg, Transpose @ Xarg]},
    Join @@ With[
      {corners0 = Xt[[All, #]]& /@ Ft},
      With[
        {grad = Sum[
          Apply[
            CalculateHarmonicAngleGradient2D,
            Append[Join @@ RotateLeft[corners0, i], A0[[i+1]]]],
          {i, 0, 2}]},
        Join @@ If[Length[Xarg] == 2, Transpose @ grad, grad]]]];
  f[Xarg_List] := With[
    {Xt = If[Length[Xarg] == 2, Xarg, Transpose @ Xarg]},
    With[
      {corners0 = Xt[[All, #]]& /@ Ft},
      0.5 / m * Sum[
        Apply[
          CalculateHarmonicAnglePotential2D,
          Append[Join @@ RotateLeft[corners0, i], A0[[i+1]]]],
        {i, 0, 2}]]];
  f];
Protect[HarmonicAnglePotential];


End[];
EndPackage[];

