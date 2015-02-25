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

ParabolicEdgePotential::usage = "ParabolicEdgePotential[mesh, order] yields a function symbol f such that f[X] is the edge potential of the given order where X is a possible vertex coordinate list for the given cortical mesh. A function of order q is calculated as the total of Abs[d - d0]^q where d is the distance between two vertices in X and d0 is the distance in the coordinates of the given mesh. Note that Grad[f, X] yields the numerical gradient of the potential at the vertex configuration given in X.
ParabolicEdgePotential[mesh] is equivalent to ParabolicEdgePotential[mesh, 2.0].";

ParabolicAnglePotential::usage = "ParabolicAnglePotential[mesh, order] yields a function symbol f such that f[X] is the angle potential of the given order where X is a possible vertex coordinate list for the given cortical mesh. A function of order q is calculated as the total of Abs[a - a0]^q where a is the angle of a face in X and a0 is the angle of the same face corner in the original mesh. Note that Grad[f, X] yields the numerical gradient of the potential at the vertex configuration given in X.
ParabolicAnglePotential[mesh] is equivalent to ParabolicEdgePotential[mesh, 2.0].";

Begin["`Private`"];

ParabolicEdgePotential[mesh_?CorticalMeshQ, order_] := With[
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
           {magnitude = order / m * Abs[norms - D0]^(order - 1)},
           SumOverEdgesDirectedTr[
             dX / {norms, norms, norms} * {magnitude, magnitude, magnitude}]]]]},
    f /: Grad[f, Xarg_List] := Apply[
      Join,
      If[Length[Xarg] == 3, df[Xarg], Transpose @ df[Transpose @ Xarg]]];
    f[Xarg_List] := (2.0 / m * Total[(# - D0)^2])&[
      If[Length@Xarg == 3, EdgeLengthsTr[mesh, Xarg], EdgeLengths[mesh, Xarg]]];
    f]];
ParabolicEdgePotential[mesh_?CorticalMapQ, order_] := With[
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
           {magnitude = order / m * Abs[norms - D0]^(order - 1)},
           SumOverEdgesDirectedTr[dX / {norms, norms} * {magnitude, magnitude}]]]]},
    f /: Grad[f, Xarg_List] := Apply[
      Join,
      If[Length[Xarg] == 2, df[Xarg], Transpose @ df[Transpose @ Xarg]]];
    f[Xarg_List] := (Total[(# - D0)^order] / m)&[
      If[Length@Xarg == 2, EdgeLengthsTr[mesh, Xarg], EdgeLengths[mesh, Xarg]]];
    f]];
ParabolicEdgePotential[mesh_?CorticalObjectQ] := ParabolicEdgePotential[mesh, 2.0];
Protect[ParabolicEdgePotential];

ParabolicAnglePotential[mesh_?CorticalMeshQ, order_] := With[
  {X0 = VertexCoordinatesTr[mesh],
   A0 = FaceAnglesTr[mesh],
   n = 3 * FaceCount[mesh],
   f = TemporarySymbol["anglePotential"]},
  With[
    {df = Function @ With[
       {A = FaceAnglesTr[mesh, #],
        dX = FaceBisectorsTr[mesh, #]},
       (* #here *)
       1.0 / n * SumOverFacesTr[
       (# / {#,#,#}&[Sqrt @ Total[#^2] - D0])& @ (#[[E[[1]]]] - #[[E[[2]]]])]]},
    f /: Grad[f, Xarg_List] := Apply[
      Join,
      If[Length[Xarg] == 3, df[Xarg], Transpose @ df[Transpose @ Xarg]]];
    f[Xarg_List] := (2.0 / n * Total[(# - D0)^2])&[
      If[Length@Xarg == 3, EdgeLengthsTr[mesh, Xarg], EdgeLengths[mesh, Xarg]]];
    f]];
ParabolicAnglePotential[mesh_?CorticalMapQ, order_] := With[
  {X0 = VertexCoordinatesTr[mesh],
   D0 = EdgeLengths[mesh],
   E = EdgePairsTr[mesh],
   m = EdgeCount[mesh],
   f = TemporarySymbol["edgePotential"]},
  With[
    {df = Function[2.0 / m * SumOverEdgesDirectedTr[
       (# / {#,#}&[Sqrt @ Total[#^2] - D0])& @ (#[[E[[1]]]] - #[[E[[2]]]])]]},
    f /: Grad[f, Xarg_List] := Apply[
      Join,
      If[Length[Xarg] == 2, df[Xarg], Transpose @ df[Transpose @ Xarg]]];
    f[Xarg_List] := (2.0 / m * Total[(# - D0)^2])&[
      If[Length@Xarg == 2, EdgeLengthsTr[mesh, Xarg], EdgeLengths[mesh, Xarg]]];
    f]];
ParabolicEdgePotential[mesh_?CorticalObjectQ] := ParabolicEdgePotential[mesh, 2.0];



End[];
EndPackage[];

