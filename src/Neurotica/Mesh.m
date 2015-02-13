(* Mesh.m
 *
 * The Neurotica`Mesh namespace contains functions related to meshes, which are either surface or
 * map objects representing topological geometries.
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
  "Neurotica`Mesh`",
  {"Neurotica`Global`", 
   "Neurotica`Util`",
   "Neurotica`Coordinates`",
   "ComputationalGeometry`"}];
Unprotect["Neurotica`Mesh`*", "Neurotica`Mesh`Private`*"];
ClearAll[ "Neurotica`Mesh`*", "Neurotica`Mesh`Private`*"];

CorticalMesh::usage = "CorticalMesh[vertexList, faceList] yields a CorticalMesh3D mesh form that can be used with the various CorticalMesh interface functions.
A cortical mesh resembles both a graph object and a boundary mesh region object. They are immutable data structures defined using the DefineImmutable function found in the Neurotica`Util` namespace; accordingly, they aren't designed to be edited, but to be copied (while reusing and sharing data between copies efficiently) via the Clone function. Cortical meshes support the following queries:
  * Properties, PropertyValue, SetProperty, RemoveProperty: All property functions supported by Graphs are also supported by cortical meshes. It is recommended that all data attached to nodes, faces, or edges be attached using properties. Note that the Properties option is accepted by CorticalMesh (see ?Properties) and Property wrappers in the vertex and face lists are parsed as well (see ?Property).
  * VertexList, VertexCount, EdgeList, EdgeCount, EdgePairs: Most graph functions can be used with a cortical mesh as if the cortical mesh were a graph. In addition, the function EdgePairs[mesh] yields the same result as EdgeList[mesh] but with lists if neighboring vertices rather than UndirectedEdge forms.
  * VertexCoordinates[mesh] additionally yields the list of vertex coordinates for the given mesh.
  * FaceList, FaceCount: In addition to edges and vertices, meshes also have faces that can be accessed with the FaceList and FaceCount functions.
  * EdgeLengths, EdgeWeights: EdgeLengths[mesh] and EdgeWeights[mesh] are identical; both yield the Euclidean distances between vertices in the mesh.
  * FaceAngles[mesh] yields the internal angles of each face in the same order as the vertices in FaceList[mesh].
  * FaceNormals[mesh] yields a list of one vector for each face in the mesh such that the vertex is orthogonal to the face.
  * FaceAxes[mesh] yields a list of orthogonal 2D axes that can be used to project each face into two dimensions, one pair for each face.
  * FaceCoordinates[mesh] yields a list of the faces but with the 3D mesh coordinates in the place of the vertex indices.
  * FaceRelativeCoordinates[mesh] yields a list of coordinates relative to the FaceAxes[mesh], such that each face has been flattened into the 2D plane it defines.
  * NeighborhoodList[mesh] yields a list of the neighbors of each vertex in the same order as VertexList[mesh].
  * NeighborhoodAngles[mesh] yields a list of the angles between the nodes in the NeighborhoodList; these angles are in the same order as the nodes in NeighborhoodList such that the first angle in a neighborhood is between the first vertex in the neighborhood, the central vertex, and second vertex in the neighborhood and the last angle in the neighborhood angles list is the angle between the first vertex in the neighborhood list, the center vertex, and the last vertex in the neighborhood list.
  * NeighborhoodBisectors[mesh] yeilds a list of the vectors that bisect each of the angles in NeighborhoodAngles[mesh].
  * NeighborhoodEdgeLengths[s] yields a list of the edge lengths for the neighborhood of each vertex in the surface or map s. The results are in the same order as NeighborhoodList[s] such that for the neighborhood list L, and the neighborhood edge length list G, G[[i,j]] is the length of the edge from vertex i to the vertex L[[i,j]].
  * EdgeIndex[mesh, e] yields the index in the EdgeList[mesh] of the edge e (e may be a list {u,v} or an UndirectedEdge[u,v]). The vertices of the edge may be in any order.
  * FaceIndex[mesh, f] yields the index in the FaceList[mesh] of the face f. The vertices in f may be in any order.
  * SourceImage[mesh] yields the source volume of the mesh (if specified).
  * MetaInformation may be passed to the cortical mesh function and may be accessed and modified via the MetaInformation[mesh] function.
  * Graph[mesh] yields a pure graph object for the mesh.
  * Most graph functions work natively with cortical meshes; e.g., FindShortestPath, BetweennessCentrality, and GraphRadius all work with cortical meshes in the first argument slot.
  * BoundaryMeshRegion[mesh] yields a pure boundary mesh object version of the mesh.
  * 3D Graphics: CorticalMesh[] accepts any option that can be passed to Graphics3D; these options will be used as default values when plotting the cortical mesh. The options may be accessed via Options[mesh] and may be changed (when using Clone) with Options -> {new-options}; Options -> Automatic will reset the options to the defaults accepted by Graphics3D with the following differences: Lighting -> \"Neutral\", ColorFunction -> None, ColorFunctionScaling -> False, Boxed -> False.";
CorticalMesh::badarg = "Bad argument given to CorticalMesh constructor: `1`";
CorticalMesh3D::usage = "CorticalMesh3D is a form used to store data for 3D surface mesh objects (see CorticalMesh).";
CorticalMeshQ::usage = "CorticalMeshQ[mesh] yields True if and only if mesh is a CorticalMesh object and False otherwise.";

CorticalMap::usage = "CorticalMap[mesh] yields a 2D flattened cortical projection of the given cortical mesh. The following options may be given:
  * Method (default: \"Equirectangular\") specifies the projection method to be used in projecting the cortical surface to a flat map. Possible values for Method include:
    * \"Mollenweide\", a Mollenweide projection, parameters: Center
    * \"Equirectangular\", a rectangular projection, parameters: Center
    * \"Mercator\", the Mercator projection, parameters: Center
    * \"Orthographic\", a projection as viewed from an infinite distance, parameters: Center
    * \"Graph\", embeds the cortex using a 2D graph embedding; note that Method -> {\"Graph\", options...} is allowed, and options are passed along to the graph embedding algorithm (see Graph, GraphPlot, GraphLayout); it is highly recommended that the Exclusions option be specified with graph embedding and that some cuts be made in order to produce a reasonable embedding.
  * Center (default: Automatic) specifies where the projection should be centered; this may be a 3D coordinate (in which case the closest point on the mesh to that coordinate is used) or a vertex identifier; if {center, orient} is given for the Center argument, then the orient point is placed on the positive x-axis in the map projection.
  * Exclusions (default: Automatic) specifies which vertices, edges, or faces are to be excluded from the projection; vertices should be specified as integer identifiers, edges as pairs (lists) or undirected edges, and faces as triples (lists) of vertex id's; in order for proper embedding, some cuts usually need to be made in the cortex; these are chosen heuristically if Automatic is specified.
  * Radius (default: Full) specifies that the radius of the projection should be limited to the given value; if this is not Full or All, then the projection algorithm excludes those faces, edges, and vertices that are farther (along the cortical surface) from the center of the projection than the given value.";
CorticalMap::badarg = "Bad argument given to CorticalMap: `1`";
CorticalMesh2D::usage = "CorticalMesh2D is a form used to store data for 2D projections of surface mesh objects (see CorticalMap).";
CorticalMapQ::usage = "CorticalMapQ[x] yields True if x is a 2D cortical projection and False otherwise.";
SourceMesh::usage = "SourceMesh[map] yields the mesh object from which the given cortical projection, map, was constructed.";

CorticalObjectQ::usage = "CorticalObjectQ[c] yields True if c is either a CorticalMesh object or a CorticalMap object and yields False otherwise.";

CortexSphericalQ::usage = "CortexSphericalQ[c] yields True if c is a CorticalMesh object that is approximately spherical.";
CortexToSphere::usage = "CortexToSphere[c] yields a CoticalMesh similar to c but inflated such that all vertices in c have a radius of approximately 100.";

VertexCoordinates::usage = "VertexCoordinates[mesh] yields the list of vertex coordinates for the given mesh in vertex order.";

FaceList::usage = "FaceList[s] yields a list of all faces in the surface or map s, if any, as lists of indices into VertexCoordinates[s].";
FaceCount::usage = "FaceCount[s] yields the count of all faces in the surface or map s.";
FaceAngles::usage = "FaceAngles[s] yields a list of the angles (in radians) of each edge in the FaceList[s] where s may be a surface or a map. FaceAngles[s, X] yields the angles for s if the vertices in s had coordinates equal to those in the list X.";
FaceIndex::usage = "FaceIndex[s] yields a list of the indices into FaceList[s] such that the i'th element in FaceIndex[s] is the list of indices at which the i'th vertex in VertexCoordinates[s] appears in FaceList[s].";
FaceNormals::usage = "FaceNormals[s] yields a list of normal vectors to each trianglular face in the surface mesh s. The vector yielded for each face is orthogonal to the plane of the face.";
FaceAxes::usage = "FaceAxe[s] yields a list of the orthogonal axes to each of the faces in the surface mesh s.";
FaceCoordinates::usage = "FaceCoordinates[s] yeilds a list of faces in which the coordinates insted of the vertex indices for each face are given.";
FaceRelativeCoordinates::usage = "FaceRelativeCoordinates[s] yeilds a list of coordinates, one for each face in the surface mesh s, that has been normalized to two dimensions.";

EdgePairs::usage = "EdgePairs[s] yields a list of the undirected edges between vertices in the surface mesh s; unlike the EdgeList function, this function yields the edges as lists instead of undirected edges.";
EdgeLengths::usage = "EdgeLengths[s] yields a list of the lengths (Euclidean norm) of each edge in the EdgeList[s] where s may be a surface mesh object or projection.";
EdgeCoordinates::usage = "EdgeCoordinates[s] yields a list identical to EdgePairs[s] except that the vertex ids in the list have been replaced with the coordinates of each vertex.";

NeighborhoodList::usage = "NeighborhoodList[s] yields a list of length N (where N is the number of vertices in s) of the neighboring vertices of each vertex; each entry k of the list is a list of the integer id's of the neighbors of the kth vertex. The neighbor id's are sorted such that they are listed in a counter-clockwise order around vertex k starting from the x-axis. The argument s may be either a map or a surface.";
NeighborhoodAngles::usage = "NeighborhoodAngles[mesh] yields a list of the angles between the nodes in the NeighborhoodList; these angles are in the same order as the nodes in NeighborhoodList such that the first angle in a neighborhood is between the first vertex in the neighborhood, the central vertex, and second vertex in the neighborhood and the last angle in the neighborhood angles list is the angle between the first vertex in the neighborhood list, the center vertex, and the last vertex in the neighborhood list.
NeighborhoodAngles[s, X] yields the neighborhood angles for s if the vertices of s were replaced with the vertices in the list X.";
NeighborhoodBisectors::usage = "NeighborhoodBisectors[mesh] yeilds a list of the vectors that bisect each of the angles in NeighborhoodAngles[mesh].
NeighborhoodBisectors[s, X] yields the neighborhood bisecting vectors for the points given by the coordinate matrix X.";
NeighborhoodEdgeLengths::usage = "NeighborhoodEdgeLengths[s] yields a list of the edge lengths for the neighborhood of each vertex in the surface or map s. The results are in the same order as NeighborhoodList[s] such that for the neighborhood list L, and the neighborhood edge length list G, G[[i,j]] is the length of the edge from vertex i to the vertex L[[i,j]].";

FaceIndex::usage = "FaceIndex[mesh, f] yields the index in the FaceList[mesh] of the face f. The vertices in f may be in any order.";
VertexFaceList::usage = "VertexFaceList[mesh] yields a list of the faces (in the same order as VertexList[mesh]) that each vertex is a member of.
VertexFaceList[mesh, vertex] yields a list of just the faces that the given vertex is a member of.";
VertexEdgeList::usage = "VertexEdgeList[mesh] yields a list of the edges (in the same order as VertexList[mesh]) that each vertex is a member of.
VertexEdgeList[mesh, vertex] yields a list of just the edges that the given vertex is a member of.";
EdgeFaceList::usage = "EdgeFaceList[mesh] yields a list of, for each edge (in the same order as EdgeList and EdgePairs), the faces to which that edge belongs.";

SourceImage::usage = "SourceImage[mesh] yields the source volume of the given mesh, if the given mesh has specified a volume; otherwise None is yielded.";

FaceRenderingFunction::usage = "FaceRenderingFunction is an option that can be given to CorticalMesh or CortexPlot3D, which specifies how to render the faces of a cortical mesh. The function is organized as with VertexRenderingFunction and EdgeRendering function; the arguments must be (1) the coordinates of the three corners of the face, (2) the vertex names for the three corners of the face, (3) the list of labels associated with the face, and (4) the list of vertex labels for the three corners of the face.";

CortexPlot3D::usage = "CortexPlot3D[mesh] yields a 3D Graphics form for the given CorticalMesh3D mesh. All options available to Graphics3D may be passed to CortexPlot3D. Note that some of the default options for Graphics3D have been altered in CortexPlot3D, and 3D graphics options that have been attached to the mesh will be used as well. See ?CorticalMesh for more details.";

CortexPlot::usage = "CortexPlot[mesh] yields a Graphics form for the given CorticalMesh2D or CorticalMesh3D mesh. If the given mesh is a 3D mesh, then the options accepted by CorticalMap are used to create the projection (Method, Center, etc.). All options available to Graphics may be passed to CortexPlot. Note that some of the default options for Graphics have been altered in CortexPlot, and 2D graphics options that have been attached to the mesh will be used as well. See ?CorticalMap for more details.";

CorticalCurvatureColor::usage = "CorticalCurvatureColor[c] yields the appropriate color for cortical curvature c in a CortexPlot or CortexPlot3D; c may be a list or a single value.";
CorticalCurvatureVertexColors::usage = "CorticalCurvatureVertexColors[m] yields the colors for each vertex in the given mesh m according to CorticalCurvatureColor[c] for the curvature c of each vertex; if there is no curvature proeprty defined, then Gray is used for each vertex.";

CalculateVertexNormals::usage = "CalculateVertexNormals[X, F] yields the normal vector to each vertex given in the 3 by n matrix of coordinates in X where the triangle faces are given by integer triples in F. These normal vectors are not normalized.";

CorticalMeshNeighborhood::usage = "CorticalMeshNeighborhood[mesh, u0, d] yields the list of vertices in the given cortical mesh such that the edge-weighted distance from u0 to each vertex returned is less than d.";

Begin["`Private`"];

(***************************************************************************************************
 * Private Helper functions
 **************************************************************************************************)

(* #NeighborhoodAnglesCompiled *)
NeighborhoodAnglesCompiled2D = Compile[{{x0, _Real, 1}, {xnei, _Real, 2}},
  If[Length[xnei] == 0,
    {},
    With[
      {x = Transpose[xnei]},
      With[
        {dx = {x[[1]] - x0[[1]], x[[2]] - x0[[2]]}},
        With[
          {norms = Sqrt[Total[dx^2]]},
          With[
            {normed = dx/{norms, norms}},
            With[
              {rot = RotateLeft /@ normed},
              ArcTan[rot[[1]], rot[[2]]] - ArcTan[normed[[1]], normed[[2]]]]]]]]],
  RuntimeOptions -> {"Speed", "EvaluateSymbolically" -> False},
  Parallelization -> True];
NeighborhoodAnglesCompiled3D = Compile[{{x0, _Real, 1}, {xnei, _Real, 2}},
  If[Length[xnei] == 0,
    {},
    With[
      {x = Transpose[xnei]},
      With[
        {zaxis = Normalize[x0],
         xaxis = Normalize[First@xnei]},
        With[
          {yaxis = Cross[zaxis, xaxis],
           dx = {x[[1]] - x0[[1]], x[[2]] - x0[[2]], x[[3]] - x0[[3]]}},
          NeighborhoodAnglesCompiled2D[{0.0, 0.0}, Dot[{xaxis, yaxis}, xnei]]]]]],
  {{NeighborhoodAnglesCompiled2D[__], _Real, 1}},
  RuntimeOptions -> {"Speed", "EvaluateSymbolically" -> False},
  Parallelization -> True];
Protect[NeighborhoodAnglesCompiled2D, NeighborhoodAnglesCompiled3D];

(* #NeighborhoodBisectorsCompiled *)
NeighborhoodBisectorsCompiled = Compile[{{x0, _Real, 1}, {xnei, _Real, 2}},
  With[
    {x = Transpose[xnei]},
    With[
      {dx = Table[x[[i]]-x0[[i]], {i,1,Length[x0]}]},
      With[
        {norms = Sqrt[Total[dx^2]]},
        With[
          {normed = dx/Table[norms, {Length[x0]}]},
          With[
            {means = 0.5*(normed + RotateLeft /@ normed)},
            With[
              {mnorms = Sqrt[Total[means^2]]},
              Transpose[means/Table[mnorms, {Length[x0]}]]]]]]]],
   RuntimeOptions -> {"Speed", "EvaluateSymbolically" -> False},
   Parallelization -> True];
Protect[NeighborhoodBisectorsCompiled];

(* #NeighborhoodEdgeLengthsCompiled *)
NeighborhoodEdgeLengthsCompiled = Compile[{{x0, _Real, 1}, {x, _Real, 2}},
  Sqrt[Total[MapThread[Subtract, {x0, Transpose[x]}]^2]],
  RuntimeOptions -> {"Speed", "EvaluateSymbolically" -> False},
  Parallelization -> True];
Protect[NeighborhoodEdgeLengthsCompiled];

(* #FaceNormalsCompiled *)
FaceNormalsCompiled = Compile[
  {{v, _Real, 2}, {cells, _Integer, 2}},
  With[
    {cellsTr = Transpose[cells]},
    With[
      {u1 = Transpose[v[[cellsTr[[2]]]] - v[[cellsTr[[1]]]]],
       u2 = Transpose[v[[cellsTr[[3]]]] - v[[cellsTr[[1]]]]]},
      With[
        {crosses = {
           u1[[2]]*u2[[3]] - u1[[3]]*u2[[2]],
           u1[[3]]*u2[[1]] - u1[[1]]*u2[[3]],
           u1[[1]]*u2[[2]] - u1[[2]]*u2[[1]]}},
        With[
          {norms = Sqrt@Total[crosses^2]},
          MapThread[
            If[#2 == 0.0, 0.0, #1/#2]&,
            {crosses, {norms, norms, norms}},
            2]]]]],
  RuntimeOptions -> {"Speed", "EvaluateSymbolically" -> False}];
Protect[FaceNormalsCompiled];

(* #FaceNormalsUnnormalizedCompiled *)
FaceNormalsUnnormalizedCompiled = Compile[
  {{v, _Real, 2}, {cells, _Integer, 2}},
  With[
    {cellsTr = Transpose[cells]},
    With[
      {u1 = Transpose[v[[cellsTr[[2]]]] - v[[cellsTr[[1]]]]],
       u2 = Transpose[v[[cellsTr[[3]]]] - v[[cellsTr[[1]]]]]},
      Transpose[{
        u1[[2]]*u2[[3]] - u1[[3]]*u2[[2]],
        u1[[3]]*u2[[1]] - u1[[1]]*u2[[3]],
        u1[[1]]*u2[[2]] - u1[[2]]*u2[[1]]}]]],
  RuntimeOptions -> {"Speed", "EvaluateSymbolically" -> False}];
Protect[FaceNormalsUnnormalizedCompiled];

(* #CalculateVertexNormals ************************************************************************)
CalculateVertexNormals[v_List, cells_List] := Check[
  With[
    {cellsTr = Transpose[cells],
     nx3 = 3 * Length[cells]},
    With[
      {u1 = Transpose[v[[cellsTr[[2]]]] - v[[cellsTr[[1]]]]],
       u2 = Transpose[v[[cellsTr[[3]]]] - v[[cellsTr[[1]]]]],
       P = SparseArray[
         Rule[
           Transpose@{Join @@ cellsTr, Range[nx3]}, 
           ConstantArray[1, nx3]]]},
      With[
        {q = Transpose[
           {u1[[2]]*u2[[3]] - u1[[3]]*u2[[2]],
            u1[[3]]*u2[[1]] - u1[[1]]*u2[[3]],
            u1[[1]]*u2[[2]] - u1[[2]]*u2[[1]]}]},
        Dot[P, Join[q, q, q]]]]],
  $Failed];

(* #FaceAxesCompiled *)
FaceAxesCompiled2D = Compile[{{x0, _Real, 1}, {xnei, _Real, 2}},
  With[
    {x = Transpose[xnei]},
    With[
      {u = (RotateLeft /@ x) - x},
      With[
        {n = Sqrt[Total[u^2]]},
        {{u[[1]]/n, u[[2]]/n}, {-u[[2]]/n, u[[1]]/n}}]]],
  RuntimeOptions -> {"Speed", "EvaluateSymbolically" -> False},
  Parallelization -> True];
FaceAxesCompiled3D = Compile[{{x0, _Real, 1}, {xnei, _Real, 2}, {tnorm, _Real, 2}},
  With[
    {x = Transpose[xnei]},
    With[
      {u = (RotateLeft /@ x) - x,
       dx = MapThread[Subtract, {x0, x}]},
      With[
        {n = Sqrt[Total[u^2]]},
        {{u[[1]]/n, u[[2]]/n, u[[3]]/n},
         {tnorm[[2]]*u[[3]]/n - tnorm[[3]]*u[[2]]/n, 
          tnorm[[3]]*u[[1]]/n - tnorm[[1]]*u[[3]]/n,
          tnorm[[1]]*u[[2]]/n - tnorm[[2]]*u[[1]]/n}}]]],
  RuntimeOptions -> {"Speed", "EvaluateSymbolically" -> False},
  Parallelization -> True];
Protect[FaceAxes2DCompiled, FaceAxes3DCompiled];

(* #TriangleCoordinatesCompiled *)
FaceCoordinatesCompiled = Compile[{{x0, _Real, 1}, {xnei, _Real, 2}, {axes, _Real, 3}},
  With[
    {x = Transpose[xnei]},
    With[
      {dx = MapThread[Subtract, {x0, x}]},
      {Total[axes[[1]] * dx], Total[axes[[2]] * dx]}]],
  RuntimeOptions -> {"Speed", "EvaluateSymbolically" -> False},
  Parallelization -> True];
Protect[FaceCoordinatesCompiled];

(* #ParseMeshProperties *)
ParseMeshProperties[X_, other_] := With[
  {others = If[$VersionNumber >= 10.0,
     Association[other],
     With[{disp = Dispatch[Append[other, k_ :> Missing["KeyAbsent", k]]]}, (# /. disp)&]],
   Xu = If[Head[X] === Property, X[[1]], X]},
  With[
    {outer = X /. {
       Property[x_List, name_ -> vals_List] /; Length[vals] == Length[x] :> {name -> vals},
       Property[x_List, rs:{Rule[_,_List]..}] /; Length[x] == Last@Dimensions@rs[[All, 2]] :> rs,
       Property[___] :> Message[CorticalMesh::badarg, "badly formatted outer property"],
       _ :> {}},
     inner = With[
       {pos = Flatten@Position[Xu, Property[__], {1}]},
       Last@Reap[
         MapThread[
           Function[{u, idx},
             Replace[others[u], (v:Except[_Missing]) :> Scan[Sow[idx -> #[[2]], #[[1]]]&, v]];
             Replace[
               u,
               {Property[_, name_ -> val_] :> Sow[idx -> val, name],
                Property[_, rs:{_Rule ..}] :> Scan[Sow[idx -> #[[2]], #[[1]]]&, rs],
                Property[___] :> Message[CorticalMesh::badarg, "badly formatted inner property"]}]],
           {Xu[[pos]], pos}],
         _,
         Rule]],
     blanks = ConstantArray[$Failed, Length[If[Head[X] === Property, X[[1]], X]]]},
    If[Length@Union[outer[[All, 1]]] != Length@outer,
      Message[CorticalMesh::badarg, "duplicate outer property lists"]];
    Fold[
      Function[{allProps, newProp},
        With[
          {name = newProp[[1]], list = newProp[[2]],
           idx = Position[allProps[[All,1]], newProp[[1]], {1}, Heads -> False]},
          If[idx == {},
            Append[allProps, name -> ReplacePart[blanks, list]],
            ReplacePart[
              allProps,
              idx[[1,1]] -> (name -> ReplacePart[allProps[[idx[[1,1]],2]], list])]]]],
      outer,
      inner]]];
Protect[ParseMeshProperties];  

(* #CorticalMeshNeighborhood *)
CorticalMeshNeighborhood[mesh_?CorticalObjectQ, u0_Integer, d_?NumericQ] := With[
  {nl = NeighborhoodEdgeLengths[mesh],
   nei = NeighborhoodList[mesh]},
  Module[
    {mark = Table[0, {VertexCount[mesh]}]},
    Last@NestWhile[
      Function[
        mark[[#[[1]]]] = 1;
        With[
          {newDisc = Join @@ (nei[[#[[1]]]]),
           newDist = Join @@ MapThread[(#2+nl[[#1]])&, #]},
          With[
            {new = Transpose@Map[
               {#[[1,1]], Min[#[[2]]]}&,
               Transpose/@Split[
                 SortBy[Transpose[{newDisc,newDist}], First],
                 (#1[[1]]==#2[[1]])&]]},
            With[
              {idcs = Pick[Range[Length@new[[1]]], mark[[new[[1]]]] + Sign[d - new[[2]]], 1]},
            {new[[1, idcs]], new[[2, idcs]]}]]]],
      {{u0}, {0.0}},
      Length[#[[1]]] > 0 &];
    Flatten@Position[mark, 1, {1}]]];
Protect[CorticalMeshNeighborhood, CorticalMeshNeighborhoodCompiled];

(* #CorticalMapTranslateCenter
 * Yields {center, orient-point}, both of which are 3D coordinates, when given a center argument;
 * The orient-point may be Automatic instead of a 3D coordinate, indicating that the first PC of the
 * resulting map should be used
 *)
CorticalMapTranslateCenter[mesh_, center_] := Check[
  Switch[
    center,
    (* center is a vertex *)
    _Integer, With[
      {idx = VertexIndex[mesh, center]},
      {If[idx === $Failed, $Failed, VertexCoordinates[mesh][[idx]]], Automatic}],
    (* center is a point *)
    {_?NumericQ, _?NumericQ, _?NumericQ}, {RegionNearest[MeshRegion[mesh], center], Automatic},
    (* center includes both a center and an orient point *)
    {_, _}, Map[
       Function[
         Switch[#,
           _Integer, With[
             {idx = VertexIndex[mesh, center]},
             If[idx === $Failed, $Failed, VertexCoordinates[mesh][[idx]]]],
           {_?NumericQ, _?NumericQ, _?NumericQ}, RegionNearest[MeshRegion[mesh], center],
           _, $Failed]],
       center],
    Automatic, {{1,0,0}, Automatic},
    (* center is not valid *)
    _, {$Failed, $Failed}],
  {$Failed, $Failed}];
Protect[CorticalMapTranslateCenter];

(* #CorticalMapAutomaticExclusions
 * Yields the lists of {vertices, edges, faces} that are to be excluded from a projection given
 * that Automatic exclusions were requested.
 *)
CorticalMapAutomaticExclusions[mesh_, method_, center_, excl_, prad_] := With[
  {X0 = Dot[
     VertexCoordinates[mesh],
     Transpose[RotationMatrix[{center[[1]], {1,0,0}}]],
     If[center[[2]] === Automatic, 
       IdentityMatrix[3], 
       RotationMatrix[-ArcTan @@ center[[2, 2;;3]], {1,0,0}]]]},
  Switch[
    If[StringQ[method], ToLowerCase[method], method],
    "mollenweide"|"equirectangular"|"mercator", Select[
      EdgePairs[mesh],
      And[X0[[#[[1]],1]] < 0, X0[[#[[2]],1]] < 0, X0[[#[[1]],2]] * X0[[#[[2]],2]] < 0] &],
    "orthographic", Part[
      VertexList[mesh],
      Select[Range[VertexCount[mesh]], X0[[#,1]] < 0&]],
    _, $Failed]];
Protect[CorticalMapAutomaticExclusions];

(* #CorticalMapTranslateExclusions
 * Yields the lists of {vertices, edges, faces} that are to be included in the projection; these
 * are only guaranteed that if you, for example, ask for VertexCoordinates[mesh][[vertices]] the
 * result will be correct; in this vein, vertices may be All.
 *)
CorticalMapTranslateExclusions[mesh_, method_, center_, excl_, prad_] := Check[
  Catch[
    With[
      {tr = Last@Reap[
         Sow[0, {1,2,3}];
         If[prad =!= Full && prad =!= All,
           (* we have a radius requirement; find everything with a small shortest path distance *)
           With[
             {u0 = First@Nearest[VertexCoordinates[mesh] -> VertexList[mesh], center[[1]]]},
             Map[
               Sow[#,1]&,
               Complement[Range[VertexCount[mesh]], CorticalMeshNeighborhood[mesh, u0, prad]]]]];
         Scan[
           Function@Which[
             IntegerQ[#],                Sow[VertexIndex[mesh, #], 1],
             ListQ[#] && Length[#] == 2, Sow[EdgeIndex[mesh, #], 2],
             Head[#] === UndirectedEdge, Sow[EdgeIndex[mesh, #], 2],
             ListQ[#] && Length[#] == 3, Sow[FaceIndex[mesh, #], 3],
             True, Throw[{$Failed, $Failed, $Failed}]],
           Which[
             ListQ[excl], Replace[
               excl, 
               Automatic :> Sequence@@CorticalMapAutomaticExclusions[
                 mesh, method, center, excl, prad],
               {1}],
             excl === Automatic, CorticalMapAutomaticExclusions[
               mesh, method, center, excl, prad],
             True, Throw[{$Failed,$Failed,$Failed}]]],
         {1,2,3},
         (Sequence @@ Union[Rest[#2]])&]},
      With[
        {Vs = Complement[Range[VertexCount[mesh]], tr[[1]]]},
        With[
          {Es = Complement[
             Range[EdgeCount[mesh]], 
             Join[
               tr[[2]],
               Flatten[VertexEdgeList[mesh, #]& /@ tr[[1]]]]]},
          With[
            {Fs = Complement[
               Range[FaceCount[mesh]],
               Join[
                 tr[[3]],
                 Flatten[VertexFaceList[mesh, #]& /@ tr[[1]]],
                 Flatten[EdgeFaceList[mesh][[tr[[2]]]]]]]},
            {If[Length[Vs] == VertexCount[mesh], All, Vs],
             If[Length[Es] == EdgeCount[mesh], All, Es],
             If[Length[Fs] == FaceCount[mesh], All, Fs]}]]]]],
  {$Failed, $Failed, $Failed}];
Protect[CorticalMapTranslateExclusions];

(* #CorticalMeshOrientForMap
 * Yields the coordinates of the subset of vertices subject to the map (via exclusions) rotated such
 * that the center is at {1,0,0} and the orient point (if there is one) is in the (+x, y)
 * half-plane.
 *)
CorticalMeshOrientForMap[mesh_?CorticalMeshQ, center_, excl_] := With[
  {RMain = If[center[[1]] === None || center[[1]] === Automatic, 
    IdentityMatrix[3],
    RotationMatrix[{center[[1]], {1,0,0}}]]},
  Dot[
    VertexCoordinates[mesh][[excl[[1]]]],
    Transpose[RMain],
    If[center[[2]] === None || center[[2]] === Automatic,
      IdentityMatrix[3], 
      With[
        {orientPt = Dot[RMain, center[[2]]]},
        Transpose[RotationMatrix[{1,0,0}, -ArcTan[orientPt[[1]], orientPt[[2]]]]]]]]];

(* #CorticalMapTranslateMethod
 * Yields a translation function given the Method argument and the translated parameters: Center,
 * vertex-exclusions, projection area, and projection radius.
 *)
CorticalMapTranslateMethod[method_, center_, excl_, prad_] := Check[
  With[
    {projFn = Replace[
       ToLowerCase[method],
       (* Here we actually define the functions for transformation;
          these should accept a list of vertices (already centered such that the center lies at 
          (1,0,0) and the orient point lies in the <positive X, Y> half-plane) and should return
          the 2D coordinates (while ignoring the orient point). *)
       {"mollenweide" :> Function[{mesh},
          With[
            {X = CorticalMeshOrientForMap[mesh, center, excl]},
            With[
              {S = Transpose@ConvertCoordinates[X, Cartesian -> {Longitude, Latitude}],
               meshRadius = Mean[Sqrt[Total[Transpose[X]^2]]]},
              With[
                {th = Block[{t}, 
                   Map[
                     Function[
                       FindRoot[
                         2.0*t + Sin[2.0*t] == Pi*Sin[#],
                         {t, #}
                         ][[1,2]]],
                     S[[2]]]]},
                Transpose[meshRadius * Sqrt[2.0] * {2.0 / Pi * S[[1]] * Cos[th], Sin[th]}]]]]],
        "equirectangular" :> Function[{mesh},
          With[
            {X = CorticalMeshOrientForMap[mesh, center, excl]},
            With[
              {S = Transpose@ConvertCoordinates[X, Cartesian -> {Longitude, Latitude}],
               meshRadius = Mean[Sqrt[Total[Transpose[X]^2]]]},
              Transpose[meshRadius * S]]]],
        "mercator" :> Function[{mesh},
          With[
            {X = CorticalMeshOrientForMap[mesh, center, excl]},
            With[
              {S = Transpose@ConvertCoordinates[X, Cartesian -> {Longitude, Latitude}],
               meshRadius = Mean[Sqrt[Total[Transpose[X]^2]]]},
              With[
                {sinPhi = Sin[S[[2]]]},
                Transpose[meshRadius * {S[[1]], 0.5*Log[(1 + sinPhi) / (1 - sinPhi)]}]]]]],
        "orthographic" :> Function[{mesh},
          With[
            {X = CorticalMeshOrientForMap[mesh, center, excl]},
            With[
              {S = Transpose@ConvertCoordinates[X, Cartesian -> {Longitude, Latitude}],
               meshRadius = Mean[Sqrt[Total[Transpose[X]^2]]]},
              Transpose[meshRadius * {Cos[S[[2]]]*Sin[S[[1]]], Sin[S[[2]]]}]]]],
        "graph" :> Function[{mesh},
          GraphEmbedding[
            Graph[
              Part[VertexList[mesh], excl[[1]]],
              Part[EdgeList[mesh], excl[[2]]],
              EdgeWeight -> Part[EdgeLengths[mesh], excl[[2]]],
              GraphLayout -> {"SpringElectricalEmbedding", "EdgeWeighted" -> True}]]],
        _ :> Message[CorticalMap::badarg, "Could not recognize projection type"]}],
     RMtxTr = Transpose@RotationMatrix[{center[[1]], {1,0,0}}]},
    With[
      {orientFn = If[center[[2]] === Automatic,
         Function[
           # . Transpose[RotationMatrix[{First@Eigenvectors[Covariance[#],1], {1,0}}]]],
         With[
           {orientRMtxTr = RotationMatrix[{projFn[{center[[2]]} . RMtxTr][[1]], {1,0}}]},
           Function[# . orientRMtxTr]]]},
      Function[orientFn @ projFn[#]]]],
  $Failed];
Protect[CorticalMapTranslateMethod];


(***************************************************************************************************
 * Exported Functions and Immutables
 **************************************************************************************************)

(* #Mesh ******************************************************************************************)

(* We define the mesh options in terms of the options to CortexPlot3D, but the default options for
   CortexPlot3D are actually all Automatic (so that they can defer to CorticalMesh). This makes it
   difficult to decide what order to define them in; to make this simpler, I've employed a private
   local variable: *)
$CortexPlot3DOptions = Join[
  FilterRules[
    Options[Graphics3D],
    Except[ColorFunction|ColorFunctionScaling|Lighting|Boxed]],
  {ColorFunction -> Automatic,
   ColorFunctionScaling -> False,
   Lighting -> "Neutral",
   Boxed -> False,
   VertexRenderingFunction -> None,
   EdgeRenderingFunction -> None,
   FaceRenderingFunction -> Automatic}];
Protect[$CortexPlot3DOptions];

(* #CorticalMesh *)
Options[CorticalMesh] = Join[
  $CortexPlot3DOptions,
  {Properties -> None,
   SourceImage -> None,
   MetaInformation -> {}}];
DefineImmutable[
  CorticalMesh[V_, F_, OptionsPattern[]] :> mesh,
  {(* First, we declare constants that depend on nothing else and thus can be used without fear of
    * forcing reevaluation of other values. These are all private.
    *)
   InitialVertexCoordinates[mesh] -> ReplaceRepeated[V, Property[x_, _] :> x],
   InitialVertexList[mesh]        -> Range[Length[InitialVertexCoordinates[mesh]]],
   InitialFaceList[mesh]          -> ReplaceRepeated[F, Property[f_, _] :> f],
   InitialEdgePairs[mesh]         -> Union@Flatten[
     {#[[{1,2}]], #[[{1,3}]], #[[{2,3}]]}& /@ Sort /@ FaceList[mesh],
     1],

   (* Options is settable, but depends on nothing but InitialOptions 
      (and should have nothing downstream but CorticalMeshQ) *)
   Options[mesh] = Map[(# -> OptionValue[#])&, $CortexPlot3DOptions[[All, 1]]],

   (* SourceImage also should not have dependencies *)
   SourceImage[mesh] = OptionValue[SourceImage],

   (* MetaInformation is a place-holder but also has no dependencies *)
   MetaInformation[mesh] = OptionValue[MetaInformation],

   (* Now we have the settable versions of the above *)
   VertexCoordinates[mesh] = InitialVertexCoordinates[mesh],
   FaceList[mesh] = InitialFaceList[mesh],
   EdgePairs[mesh] = InitialEdgePairs[mesh],
   EdgePairs[mesh, patt_] := Cases[EdgePairs[mesh], patt, {1}],
   FaceList[mesh, patt_] := Cases[FaceList[mesh], patt, {1}],

   (* VertexList is not settable, but it is public thus depends on the current vertex coordinates *)
   VertexList[mesh] -> Range[Length[VertexCoordinates[mesh]]],
   VertexList[mesh, patt_] := Cases[VertexList[mesh], patt, {1}],

   (* Here we have indices; the index arrays are private, but the index functions are public *)

   (* #FaceIndexArray [private] *)
   FaceIndexArray[mesh] :> SparseArray[
     Flatten[
       MapIndexed[
         Function[
           {#1 -> #2[[1]],
            {#1[[2]], #1[[3]], #1[[1]]} -> #2[[1]],
            {#1[[3]], #1[[1]], #1[[2]]} -> #2[[1]],
            {#1[[2]], #1[[1]], #1[[3]]} -> #2[[1]],
            {#1[[1]], #1[[3]], #1[[2]]} -> #2[[1]],
            {#1[[3]], #1[[2]], #1[[1]]} -> #2[[1]]}],
         FaceList[mesh]],
       1],
     {Max[VertexList[mesh]], Max[VertexList[mesh]], Max[VertexList[mesh]]},
     0],
   (* #EdgeIndexArray [private] *)
   EdgeIndexArray[mesh] :> SparseArray[
     Flatten[
       MapIndexed[
         Function[{#1 -> #2[[1]], Reverse[#1] -> #2[[1]]}],
         EdgePairs[mesh]],
       1],
     {Max[VertexList[mesh]], Max[VertexList[mesh]]},
     0],
   (* VertexIndexArray [private] *)
   VertexIndexArray[mesh] :> If[VertexList[mesh] == Range[Length@VertexList[mesh]],
     VertexList[mesh],
     Normal @ SparseArray[
       VertexList[mesh] -> Range[VertexCount[mesh]],
       Max[VertexList[mesh]],
       0]],

   (* #FaceIndex *)
   FaceIndex[mesh, {a_Integer, b_Integer, c_Integer}] := With[
     {id = FaceIndexArray[mesh][[a,b,c]]},
     If[id == 0, $Failed, id]],
   (* EdgeIndex *)
   EdgeIndex[mesh, (List|UndirectedEdge)[a_Integer, b_Integer]] := With[
     {id = EdgeIndexArray[mesh][[a,b]]},
     If[id == 0, $Failed, id]],     
   (* VertexIndex *)
   VertexIndex[mesh, i_Integer] := With[
     {id = If[i > Length[VertexIndexArray[mesh]], 0, VertexIndexArray[mesh][[i]]]},
     If[id == 0, $Failed, id]],

   (* These indices depend only on the face list and edge pairs and tell us to which faces a vertex
    * or an edge belongs.
    *)
   VertexEdgeList[mesh] :> Last@Reap[
     MapIndexed[
       Function[Sow[#2[[1]], #1]],
       EdgePairs[mesh]],
     VertexList[mesh],
     (Sequence@@#2)&],
   VertexFaceList[mesh] :> Last@Reap[
     MapIndexed[
       Function[Sow[#2[[1]], #1]],
       FaceList[mesh]],
     VertexList[mesh],
     (Sequence@@#2)&],
   VertexEdgeList[mesh, i_Integer] := Part[VertexEdgeList[mesh], VertexIndex[mesh, i]],
   VertexFaceList[mesh, i_Integer] := Part[VertexFaceList[mesh], VertexIndex[mesh, i]],
   EdgeFaceList[mesh] :> With[
     {Ft = Transpose @ FaceList[mesh]},
     Last @ Reap[
       MapThread[
         Sow[#1, Sort/@{{#2,#3}, {#2,#4}, {#3,#4}}]&,
         {Range[Length@First@Ft], Ft[[1]], Ft[[2]], Ft[[3]]}],
       EdgePairs[mesh],
       (Sequence @@ #2)&]],
   EdgeFaceList[mesh, e_] := Part[EdgeFaceList[mesh], EdgeIndex[mesh, e]],

   (* In this section, we deal with options.
    * Options must not depend on current values and instead must be changed via the SetOptions
    * interface, which is defined below, after the body of the immutable definition.
    *)

   (* #OptionalProperties [private] *)
   OptionalProperties[mesh] -> Last@Reap[
     Replace[
       OptionValue[Properties] /. None -> {},
       {r:Rule[_Integer, {_Rule..}] :> Sow[r, "Vertex"],
        r:Rule[{_Integer, _Integer}, {_Rule..}] :> Sow[r, "Edge"],
        Rule[UndirectedEdge[a_Integer, b_Integer], rs:{_Rule..}] :> Sow[{a,b} -> rs, "Edge"],
        r:Rule[{_Integer,_Integer,_Integer}, {_Rule..}] :> Sow[r, "Face"],
        x_ :> Message[CorticalMesh::badarg, "badly formatted property option"]},
       {1}],
     {"Vertex", "Edge", "Face"},
     (Sequence@@#2)&],

   (* #VertexProperties [private] *)
   VertexProperties[mesh] = ParseMeshProperties[V, OptionalProperties[mesh][[1]]],
   (* Note that VertexProperties does not depend on VertexCoordinates or FaceList so does not get
      recalculated when someone changes the coordinates, for example. This lets us intelligently
      carry over properties between clones. *)

   (* #EdgeProperties [private] *)
   EdgeProperties[mesh] = ParseMeshProperties[
     InitialEdgePairs[mesh],
     OptionalProperties[mesh][[2]]],

   (* #FaceProperties [private] *)
   FaceProperties[mesh] = ParseMeshProperties[F, OptionalProperties[mesh][[3]]],

   (* Now that those private mutables have been setup without dependence on the current state of
    * variables like VertexCoordinates, we can make the Properties value depend on them, so that
    * an update of any of them will result in an update of the properties and properties list, etc.
    *)

   (* #Properties *)
   Properties[mesh] :> With[
     {allPropNames = Union[
        VertexProperties[mesh][[All,1]],
        EdgeProperties[mesh][[All,1]],
        FaceProperties[mesh][[All,1]]]},
     Map[
       Function[{name},
         name -> MapThread[
           If[#1 === $Failed, #2 -> #1, ReplacePart[#1, 1 -> #2]]&,
           {{FirstCase[VertexProperties[mesh], (Rule|RuleDelayed)[name,val_],$Failed],
             FirstCase[EdgeProperties[mesh], (Rule|RuleDelayed)[name,val_], $Failed],
             FirstCase[FaceProperties[mesh], (Rule|RuleDelayed)[name,val_], $Failed]},
            {VertexList, EdgeList, FaceList}}]],
       allPropNames]],
   (* #PropertyList *)
   PropertyList[mesh] := Join[Properties[mesh][[All, 1]], {EdgeWeight, VertexCoordinates}],
   (* Note that property value and set property are below, outside the immutable definition *)

   (* #VertexCount *)
   VertexCount[mesh] -> Length[VertexCoordinates[mesh]],
   VertexCount[mesh, patt_] := Count[VertexList[mesh], patt, {1}],
   (* #EdgeCount *)
   EdgeCount[mesh] -> Length[EdgePairs[mesh]],
   EdgeCount[mesh, patt_] := Count[EdgePairs[mesh], patt, {1}],
   (* #FaceCount *)
   FaceCount[mesh] -> Length[FaceList[mesh]],
   FaceCount[mesh, patt_] := Count[FaceList[mesh], patt, {1}],
   (* #VertexDegree *)
   VertexDegree[mesh] :> Length /@ NeighborhoodList[mesh],
   VertexInDegree[mesh] := VertexDegree[mesh],
   VertexOutDegree[mesh] := VertexDegree[mesh],
   VertexDegree[mesh, i_Integer] := VertexDegree[mesh][[VertexIndex[mesh, i]]],
   VertexInDegree[mesh, i_Integer] := VertexDegree[mesh][[VertexIndex[mesh, i]]],
   VertexOutDegree[mesh, i_Integer] := VertexDegree[mesh][[VertexIndex[mesh, i]]],

   (* #FaceAngles *)
   FaceAngles[mesh] :> With[
     {FF = Transpose[FaceList[mesh]],
      XX = VertexCoordinates[mesh]},
     With[
       {Xf = XX[[#]]& /@ FF},
       With[
         {dX = Transpose /@ {
            (Xf[[2]] - Xf[[1]]),
            (Xf[[3]] - Xf[[2]]),
            (Xf[[1]] - Xf[[3]])}},
         With[
           {normed = Map[
              With[{lens = Sqrt[Total[#^2]]}, (#/lens)& /@ #]&,
              dX]},
           Transpose[
             ArcCos[
               {Total[normed[[1]] * -normed[[3]]],
                Total[normed[[2]] * -normed[[1]]],
                Total[normed[[3]] * -normed[[2]]]}]]]]]],
   (* #FaceNormals *)
   FaceNormals[mesh] :> FaceNormalsCompiled[VertexList[mesh], FaceList[mesh]],
   FaceNormals[mesh, X_] := FaceNormalsCompiled[X, FaceList[mesh]],
   (* #FaceAxes *)
   FaceAxes[mesh] :> With[
     {X = VertexCoordinates[mesh],
      nei = NeighborhoodList[mesh]},
     With[
       {Xnei = X[[#]]& /@ nei},
       MapThread[FaceAxesCompiled3D, {X, Xnei, MapThread[TriangleNormalsCompiled, {X, Xnei}]}]]],
   (* #FaceCoordinates *)
   FaceCoordinates[mesh] :> With[
     {X = VertexCoordinates[mesh],
      FF = FaceList[mesh]},
     X[[#]]& /@ FF],
   (* #FaceRelativeCoordinates *)
   FaceRelativeCoordinates[mesh] :> With[
     {X = VertexCoordinates[mesh],
      nei = NeighborhoodList[mesh],
      Fn = FaceNormals[mesh]},
     With[
       {Xnei = X[[#]]& /@ nei},
       MapThread[
         TriangleCoordinatesCompiled,
         {X, Xnei, MapThread[FaceAxesCompiled3D, {X, Xnei, Fn}]}]]],

   (* #EdgeList *)
   EdgeList[mesh] :> Map[Apply[UndirectedEdge, #]&, EdgePairs[mesh]],
   EdgeList[mesh, patt_] := Count[EdgeList[mesh], patt, {1}],
   (* #EdgeLengths *)
   EdgeLengths[mesh] :> With[
     {EL = Transpose[EdgePairs[mesh]],
      X = VertexCoordinates[mesh]},
     With[
       {X1 = Transpose[X[[EL[[1]]]]],
        X2 = Transpose[X[[EL[[2]]]]]},
       Sqrt[
         Plus[
           (X1[[1]] - X2[[1]])^2,
           (X1[[2]] - X2[[2]])^2,
           (X1[[3]] - X2[[3]])^2]]]],
   EdgeWeight[mesh] := EdgeLengths[mesh],
   EdgeWeight[mesh, e:(List|UndirectedEdge)[_Integer, _Integer]] := Part[
     EdgeLengths[mesh],
     EdgeIndex[e]],
   EdgeWeight[mesh, es:{(List|UndirectedEdge)[_Integer,_Integer]..}] := Part[
     EdgeLengths[mesh],
     EdgeIndex[mesh, #]& /@ es],
   (* #EdgeCoordinates *)
   EdgeCoordinates[mesh] :> With[{X = VertexCoordinates[mesh]}, X[[#]]& /@ EdgePairs[mesh]],
   
   (* #NeighborhoodList *)
   NeighborhoodList[mesh] :> With[
     {X = VertexCoordinates[mesh],
      E = EdgePairs[mesh]},
     Last[
       Reap[
         Scan[
           Function[{pair},
             Sow[pair[[1]], pair[[2]]];
             Sow[pair[[2]], pair[[1]]]],
           E],
         Range[Length[X]],
         Function[{id, neighbors},
           Sequence @@ With[
             {neis = Union[neighbors]},
             With[
               {U = Dot[
                  X[[neis]], 
                  Transpose[RotationMatrix[{X[[id]], {0,0,1}}]]
                  ][[All, 1;;2]]},
               SortBy[Thread[neis -> U], ArcTan[#[[2,1]], #[[2,2]]]&][[All,1]]]]]]]],
   (* #NeighborhoodAngles *)
   NeighborhoodAngles[mesh] :> With[
     {X = VertexCoordinates[mesh],
      nei = NeighborhoodList[mesh]},
     MapThread[NeighborhoodAnglesCompiled3D, {X, X[[#]]& /@ nei}]],
   NeighborhoodAngles[mesh, X_] := MapThread[
     NeighborhoodAnglesCompiled3D, 
     {X, X[[#]]& /@ NeighborhoodList[mesh]}],
   (* #NeighborhoodBisectors *)
   NeighborhoodBisectors[mesh] :> With[
     {X = VertexCoordinates[mesh],
      nei = NeighborhoodList[mesh]},
     MapThread[NeighborhoodBisectorsCompiled, {X, X[[#]] & /@ nei}]],
   NeighborhoodBisectors[mesh, X_] := MapThread[
     NeighborhoodBisectorsCompiled, 
     {X, X[[#]]& /@ NeighborhoodList[mesh]}],
   (* #NeighborhoodEdgeLengths *)
   NeighborhoodEdgeLengths[mesh] :> With[
     {X = VertexCoordinates[mesh],
      nei = NeighborhoodList[mesh]},
     MapThread[NeighborhoodEdgeLengthsCompiled, {X, X[[#]]& /@ nei}]],
   NeighborhoodEdgeLengths[mesh, X_] := MapThread[
     NeighborhoodEdgeLengthsCompiled, 
     {X, X[[#]]& /@ NeighborhoodList[mesh]}],


   (* #VertexNormals *)
   VertexNormals[mesh] :> CalculateVertexNormals[VertexCoordinates[mesh], FaceList[mesh]],
   VertexNormals[mesh, X_] := CalculateVertexNormals[X, FaceList[mesh]],
     
   (* #MeshRegion *)
   MeshRegion[mesh] :> MeshRegion[VertexCoordinates[mesh], Polygon[FaceList[mesh]]],

   (* #BoundaryMeshRegion *)
   BoundaryMeshRegion[mesh] :> BoundaryMeshRegion[VertexCoordinates[mesh], Polygon[FaceList[mesh]]],
   TriangulateMesh[mesh, opts___] := TriangulateMesh[BoundaryMeshRegion[mesh], opts],
   HighlightMesh[mesh, opts___] := HighlightMesh[BoundaryMeshRegion[mesh], opts],
   DimensionalMeshComponents[mesh, opts___] := DimensionalMeshComponents[
     BoundaryMeshRegion[mesh], 
     opts],
   ConnectedMeshComponents[mesh, opts___] := ConnectedMeshComponents[
     BoundaryMeshRegion[mesh], 
     opts],
   MeshCoordinates[mesh, opts___] := MeshCoordinates[BoundaryMeshRegion[mesh], opts],
   MeshCells[mesh, opts___] := MeshCells[BoundaryMeshRegion[mesh], opts],
   MeshPrimitives[mesh, opts___] := MeshPrimitives[BoundaryMeshRegion[mesh], opts],
   MeshCellIndex[mesh, opts___] := MeshCellIndex[BoundaryMeshRegion[mesh], opts],
   MeshCellCount[mesh, opts___] := MeshCellCount[BoundaryMeshRegion[mesh], opts],
   RegionDimension[mesh, opts___] := RegionDimension[BoundaryMeshRegion[mesh], opts],
   RegionEmbeddingDimension[mesh, opts___] := RegionEmbeddingDimension[
     BoundaryMeshRegion[mesh],
     opts],
   RegionMeasure[mesh, opts___] := RegionMeasure[BoundaryMeshRegion[mesh], opts],
   RegionNearest[mesh, opts___] := RegionNearest[BoundaryMeshRegion[mesh], opts],
   RegionDistance[mesh, opts___] := RegionDistance[BoundaryMeshRegion[mesh], opts],
   SignedRegionDistance[mesh, opts___] := SignedRegionDistance[BoundaryMeshRegion[mesh], opts],
   RegionMember[mesh, opts___] := RegionMember[BoundaryMeshRegion[mesh], opts],
   RegionBounds[mesh, opts___] := RegionBounds[BoundaryMeshRegion[mesh], opts],
   RegionBoundary[mesh, opts___] := RegionBoundary[BoundaryMeshRegion[mesh], opts],
   RegionCentroid[mesh, opts___] := RegionCentroid[BoundaryMeshRegion[mesh], opts],
   Volume[mesh, opts___] := Volume[BoundaryMeshRegion[mesh], opts],

   (* #Graph *)
   Graph[mesh] :> Graph[
     VertexList[mesh],
     EdgeList[mesh],
     VertexCoordinates -> VertexCoordinates[mesh],
     EdgeWeight -> EdgeLengths[mesh]],
   BetweennessCentrality[mesh, opts___]       := BetweennessCentrality[Graph[mesh], opts],
   ClosenessCentrality[mesh, opts___]         := ClosenessCentrality[Graph[mesh], opts],
   DegreeCentrality[mesh, opts___]            := DegreeCentrality[Graph[mesh], opts],
   EdgeBetweennessCentrality[mesh, opts___]   := EdgeBetweennessCentrality[Graph[mesh], opts],
   EdgeConnectivity[mesh, opts___]            := EdgeConnectivity[Graph[mesh], opts],
   EdgeCycleMatrix[mesh, opts___]             := EdgeCycleMatrix[Graph[mesh], opts],
   EigenvectorCentrality[mesh, opts___]       := EigenvectorCentrality[Graph[mesh], opts],
   EulerianGraphQ[mesh]                       := EulerianGraphQ[Graph[mesh]],
   FindCycle[mesh, opts___]                   := FindCycle[Graph[mesh], opts],
   FindEdgeIndendentPaths[mesh, opts___]      := FindEdgeIndendentPaths[Graph[mesh], opts],
   FindEulerianCycle[mesh, opts___]           := FindEulerianCycle[Graph[mesh], opts],
   FindFundamentalCycles[mesh, opts___]       := FindFundamentalCycles[Graph[mesh], opts],
   FindHamiltonianCycle[mesh, opts___]        := FindHamiltonianCycle[Graph[mesh], opts],
   FindMaximumFlow[mesh, opts___]             := FindMaximumFlow[Graph[mesh], opts],
   FindMinimumCostFlow[mesh, opts___]         := FindMinimumCostFlow[Graph[mesh], opts],
   FindMinimumFlow[mesh, opts___]             := FindMinimumFlow[Graph[mesh], opts],
   FindPath[mesh, opts___]                    := FindPath[Graph[mesh], opts],
   FindPostmanTour[mesh, opts___]             := FindPostmanTour[Graph[mesh], opts],
   FindShortestPath[mesh, opts___]            := FindShortestPath[Graph[mesh], opts],
   FindShortestTour[mesh, opts___]            := FindShortestTour[Graph[mesh], opts],
   FindVertexIndendentPaths[mesh, opts___]    := FindVertexIndendentPaths[Graph[mesh], opts],
   GraphAssortativity[mesh, opts___]          := GraphAssortativity[Graph[mesh], opts],
   GraphCenter[mesh, opts___]                 := GraphCenter[Graph[mesh], opts],
   GraphDiameter[mesh, opts___]               := GraphDiameter[Graph[mesh], opts],
   GraphDistance[mesh, opts___]               := GraphDistance[Graph[mesh], opts],
   GraphDistanceMatrix[mesh, opts___]         := GraphDistanceMatrix[Graph[mesh], opts],
   GraphPeriphery[mesh, opts___]              := GraphPeriphery[Graph[mesh], opts],
   GraphPower[mesh, opts___]                  := GraphPower[Graph[mesh], opts],
   GraphRadius[mesh, opts___]                 := GraphRadius[Graph[mesh], opts],
   GraphReciprocity[mesh, opts___]            := GraphReciprocity[Graph[mesh], opts],
   GlobalClusteringCoefficient[mesh, opts___] := GlobalClusteringCoefficient[mesh, opts],
   HamiltonianGraphQ[mesh]                    := HamiltonianGraphQ[Graph[mesh]],
   HITSCentrality[mesh, opts___]              := HITSCentrality[Graph[mesh], opts],
   KatzCentrality[mesh, opts___]              := KatzCentrality[Graph[mesh], opts],
   LocalClusteringCoefficient[mesh, opts___]  := LocalClusteringCoefficient[Graph[mesh], opts],
   MeanClusteringCoefficient[mesh, opts___]   := MeanClusteringCoefficient[Graph[mesh], opts],
   MeanDegreeConnectivity[mesh, opts___]      := MeanDegreeConnectivity[Graph[mesh], opts],
   MeanNeighborDegree[mesh, opts___]          := MeanNeighborDegree[Graph[mesh], opts],
   PageRankCentrality[mesh, opts___]          := PageRankCentrality[Graph[mesh], opts],
   RadialityCentrality[mesh, opts___]         := RadialityCentrality[Graph[mesh], opts],
   ShortestPathFunction[mesh, opts___]        := ShortestPathFunction[Graph[mesh], opts],
   StasusCentrality[mesh, opts___]            := StasusCentrality[Graph[mesh], opts],
   VertexConnectivity[mesh, opts___]          := VertexConnectivity[Graph[mesh], opts],
   VertexCorrelationSimilarity[mesh, opts___] := VertexCorrelationSimilarity[Graph[mesh], opts],
   VertexCosineSimilarity[mesh, opts___]      := VertexCosineSimilarity[Graph[mesh], opts],
   VertexDiceSimilarity[mesh, opts___]        := VertexDiceSimilarity[Graph[mesh], opts],
   VertexEccentricity[mesh, opts___]          := VertexEccentricity[Graph[mesh], opts],
   VertexInDegree[mesh, opts___]              := VertexDegree[mesh],
   VertexJaccardSimilarity[mesh, opts___]     := VertexJaccardSimilarity[Graph[mesh], opts],
   VertexOutDegree[mesh, opts___]             := VertexDegree[mesh],

   (* #CorticalMeshQ *)
   CorticalMeshQ[mesh] -> Which[
     !ArrayQ[VertexCoordinates[mesh], 2, NumericQ], Message[
       CorticalMesh::badarg,
       "VertexCoordinates must be a 2D numeric array"],
     !MatchQ[Dimensions[VertexCoordinates[mesh]], {_, 2|3}], Message[
       CorticalMesh::badarg,
       "VertexCoordinates must contain 2D or 3D points"],
     !ArrayQ[FaceList[mesh], 2, IntegerQ], Message[
       CorticalMesh::badarg, 
       "FaceList conatins non-integer values"],
     Min[Flatten@FaceList[mesh]] < 1, Message[
       CorticalMesh::badarg, 
       "FaceList conatins values less than 1"],
     Max[Flatten@FaceList[mesh]] > Length[VertexCoordinates[mesh]], Message[
       CorticalMesh::badarg, 
       "FaceList conatins values greater than the number of vertices"],
     Count[FaceList[mesh], f_List /; Length[Union[f]] != 3, {1}] > 0, Message[
       CorticalMesh::badarg,
       "FaceList must contain triangles only"],
     And[
       Options[mesh] =!= Automatic,
       Complement[Options[mesh][[All,1]], $CortexPlot3DOptions[[All, 1]]] != {}], Message[
         CorticalMesh::badarg,
         "Unrecognized option given to CorticalMesh"],
     True, True]},
  SetSafe -> True,
  Symbol -> CorticalMesh3D];

(* #CorticalMeshQ *)
CorticalMeshQ[_] := False;



(**************************************************************************************************)
(* #CorticalMap *)

(* This is essentially identical to the treatment of CortexPlot3D arguments above.
   We define the mesh options in terms of the options to CortexPlot, but the default options for
   CortexPlot are actually all Automatic (so that they can defer to CorticalMap). This makes
   it difficult to decide what order to define them in; to make this simpler, I've employed a
   private local variable: *)
$CortexPlotOptions = Join[
  FilterRules[
    Options[Graphics],
    Except[ColorFunction|ColorFunctionScaling]],
  {ColorFunction -> Automatic,
   ColorFunctionScaling -> False,
   VertexRenderingFunction -> None,
   EdgeRenderingFunction -> None,
   FaceRenderingFunction -> Automatic}];
Protect[$CortexPlot3DOptions];

Options[CorticalMap] = Join[
  $CortexPlotOptions,
  {MetaInformation -> {},
   Properties -> None,
   Method -> "Equirectangular",
   Center -> Automatic,
   Exclusions -> Automatic,
   Radius -> Full}];
DefineImmutable[
  CorticalMap[mesh_?CorticalMeshQ, OptionsPattern[]] :> map,
  {(* First we declare the simple constants for the projection *)
   SourceMesh[map] = mesh,

   (* Options is settable, but depends on nothing but the initial options
      (and should have nothing downstream but CorticalMapQ); 
      note that these are the graphics options only *)
   Options[map] = Map[(# -> OptionValue[#])&, $CortexPlotOptions[[All, 1]]],

   (* MetaInformation is a place-holder but also has no dependencies *)
   MetaInformation[map] = OptionValue[MetaInformation],

   (* Here we add options that allow the user to edit the projection details *)
   Method[map] = OptionValue[Method],
   Center[map] = OptionValue[Center],
   Exclusions[map] = OptionValue[Exclusions],
   Radius[map] = OptionValue[Radius],

   (* These ones are private... *)
   TranslatedCenter[map] -> CorticalMapTranslateCenter[SourceMesh[map], Center[map]],
   TranslatedExclusions[map] -> CorticalMapTranslateExclusions[
     SourceMesh[map],
     Method[map],
     TranslatedCenter[map],
     Exclusions[map],
     Radius[map]],
   TransformationFunction[map] -> CorticalMapTranslateMethod[
     Method[map],
     TranslatedCenter[map],
     TranslatedExclusions[map],
     Radius[map]],

   (* Now we have the settable versions of the above *)
   VertexCoordinates[map] -> With[
     {fn = TransformationFunction[map]},
     fn[SourceMesh[map]]],
   FaceList[map] -> FaceList[SourceMesh[map]][[TranslatedExclusions[map][[3]]]],
   EdgePairs[map] -> EdgePairs[SourceMesh[map]][[TranslatedExclusions[map][[2]]]],
   EdgePairs[map, patt_] := Cases[EdgePairs[map], patt, {1}],
   FaceList[map, patt_] := Cases[FaceList[map], patt, {1}],

   (* VertexList is not settable, but it is public thus depends on the current vertex coordinates *)
   VertexList[map] -> VertexList[SourceMesh[map]][[TranslatedExclusions[map][[1]]]],
   VertexList[map, patt_] := Cases[VertexList[map], patt, {1}],

   (* Here we have indices; the index arrays are private, but the index functions are public *)

   (* #FaceIndexArray [private] *)
   FaceIndexArray[map] :> SparseArray[
     Flatten[
       MapIndexed[
         Function[
           {#1 -> #2[[1]],
            {#1[[2]], #1[[3]], #1[[1]]} -> #2[[1]],
            {#1[[3]], #1[[1]], #1[[2]]} -> #2[[1]],
            {#1[[2]], #1[[1]], #1[[3]]} -> #2[[1]],
            {#1[[1]], #1[[3]], #1[[2]]} -> #2[[1]],
            {#1[[3]], #1[[2]], #1[[1]]} -> #2[[1]]}],
         FaceList[map]],
       1],
     {Max[VertexList[map]], Max[VertexList[map]], Max[VertexList[map]]},
     0],
   (* #EdgeIndexArray [private] *)
   EdgeIndexArray[map] :> SparseArray[
     Flatten[
       MapIndexed[
         Function[{#1 -> #2[[1]], Reverse[#1] -> #2[[1]]}],
         EdgePairs[map]],
       1],
     {Max[VertexList[map]], Max[VertexList[map]]},
     0],
   (* VertexIndexArray [private] *)
   VertexIndexArray[map] :> If[VertexList[map] == Range[Length@VertexList[map]],
     VertexList[map],
     Normal @ SparseArray[
       VertexList[map] -> Range[VertexCount[map]],
       Max[VertexList[map]],
       0]],

   (* #FaceIndex *)
   FaceIndex[map, {a_Integer, b_Integer, c_Integer}] := With[
     {id = FaceIndexArray[map][[a,b,c]]},
     If[id == 0, $Failed, id]],
   (* EdgeIndex *)
   EdgeIndex[map, (List|UndirectedEdge)[a_Integer, b_Integer]] := With[
     {id = EdgeIndexArray[map][[a,b]]},
     If[id == 0, $Failed, id]],     
   (* VertexIndex *)
   VertexIndex[map, i_Integer] := With[
     {id = If[i > Length[VertexIndexArray[map]], 0, VertexIndexArray[map][[i]]]},
     If[id == 0, $Failed, id]],

   (* These indices depend only on the face list and edge pairs and tell us to which faces a vertex
    * or an edge belongs.
    *)
   VertexEdgeList[map] :> Last@Reap[
     MapIndexed[
       Function[Sow[#2[[1]], #1]],
       EdgePairs[map]],
     VertexList[map],
     (Sequence@@#2)&],
   VertexFaceList[map] :> Last@Reap[
     MapIndexed[
       Function[Sow[#2[[1]], #1]],
       FaceList[map]],
     VertexList[map],
     (Sequence@@#2)&],
   VertexEdgeList[map, i_Integer] := Part[VertexEdgeList[map], VertexIndex[map, i]],
   VertexFaceList[map, i_Integer] := Part[VertexFaceList[map], VertexIndex[map, i]],

   (* In this section, we deal with options.
    * Options must not depend on current values and instead must be changed via the SetOptions
    * interface, which is defined below, after the body of the immutable definition.
    *)

   (* #OptionalProperties [private] identical to that of CorticalMesh *)
   OptionalProperties[map] -> Last@Reap[
     Replace[
       OptionValue[Properties] /. None -> {},
       {r:Rule[_Integer, {_Rule..}] :> Sow[r, "Vertex"],
        r:Rule[{_Integer, _Integer}, {_Rule..}] :> Sow[r, "Edge"],
        Rule[UndirectedEdge[a_Integer, b_Integer], rs:{_Rule..}] :> Sow[{a,b} -> rs, "Edge"],
        r:Rule[{_Integer,_Integer,_Integer}, {_Rule..}] :> Sow[r, "Face"],
        x_ :> Message[CorticalMap::badarg, "badly formatted property option"]},
       {1}],
     {"Vertex", "Edge", "Face"},
     (Sequence@@#2)&],

   (* Note that these Properties literally just steal the mesh's properties *)

   (* #VertexProperties [private] *)
   VertexProperties[map] = With[
     {idx = TranslatedExclusions[map][[1]]},
     Map[
       (#[[1]] -> #[[2]][[idx]]) &,
       VertexProperties[SourceMesh[map]]]],
   (* #EdgeProperties [private] *)
   EdgeProperties[map] = With[
     {idx = TranslatedExclusions[map][[2]]},
     Map[
       (#[[1]] -> #[[2]][[idx]]) &,
       EdgeProperties[SourceMesh[map]]]],
   (* #FaceProperties [private] *)
   FaceProperties[map] = With[
     {idx = TranslatedExclusions[map][[3]]},
     Map[
       (#[[1]] -> #[[2]][[idx]]) &,
       FaceProperties[SourceMesh[map]]]],

   (* Now that those private mutables have been setup without dependence on the current state of
    * variables like VertexCoordinates, we can make the Properties value depend on them, so that
    * an update of any of them will result in an update of the properties and properties list, etc.
    *)

   (* #Properties *)
   Properties[map] :> With[
     {allPropNames = Union[
        VertexProperties[map][[All,1]],
        EdgeProperties[map][[All,1]],
        FaceProperties[map][[All,1]]]},
     Map[
       Function[{name},
         name -> MapThread[
           If[#1 === $Failed, #2 -> #1, ReplacePart[#1, 1 -> #2]]&,
           {{FirstCase[VertexProperties[map], (Rule|RuleDelayed)[name,val_], $Failed],
             FirstCase[EdgeProperties[map], (Rule|RuleDelayed)[name,val_], $Failed],
             FirstCase[FaceProperties[map], (Rule|RuleDelayed)[name,val_], $Failed]},
            {VertexList, EdgeList, FaceList}}]],
       allPropNames]],
   (* #PropertyList *)
   PropertyList[map] := Join[Properties[map][[All, 1]], {EdgeWeight, VertexCoordinates}],
   (* Note that property value and set property are below, outside the immutable definition *)

   (* #VertexCount *)
   VertexCount[map] -> Length[VertexCoordinates[map]],
   VertexCount[map, patt_] := Count[VertexList[map], patt, {1}],
   (* #EdgeCount *)
   EdgeCount[map] -> Length[EdgePairs[map]],
   EdgeCount[map, patt_] := Count[EdgePairs[map], patt, {1}],
   (* #FaceCount *)
   FaceCount[map] -> Length[FaceList[map]],
   FaceCount[map, patt_] := Count[FaceList[map], patt, {1}],
   (* #VertexDegree *)
   VertexDegree[map] :> Length /@ NeighborhoodList[map],
   VertexInDegree[map] := VertexDegree[map],
   VertexOutDegree[map] := VertexDegree[map],
   VertexDegree[map, i_Integer] := VertexDegree[map][[VertexIndex[map, i]]],
   VertexInDegree[map, i_Integer] := VertexDegree[map][[VertexIndex[map, i]]],
   VertexOutDegree[map, i_Integer] := VertexDegree[map][[VertexIndex[map, i]]],

   (* #FaceAngles *)
   FaceAngles[map] :> With[
     {FF = Transpose[FaceList[map]],
      XX = VertexCoordinates[map]},
     With[
       {Xf = XX[[#]]& /@ FF},
       With[
         {dX = Transpose /@ {
            (Xf[[2]] - Xf[[1]]),
            (Xf[[3]] - Xf[[2]]),
            (Xf[[1]] - Xf[[3]])}},
         With[
           {normed = Map[
              With[{lens = Sqrt[Total[#^2]]}, (#/lens)& /@ #]&,
              dX]},
           Transpose[
             ArcCos[
               {Total[normed[[1]] * -normed[[3]]],
                Total[normed[[2]] * -normed[[1]]],
                Total[normed[[3]] * -normed[[2]]]}]]]]]],
   (* #FaceAxes *)
   FaceAxes[map] :> With[
     {X = VertexCoordinates[map],
      nei = NeighborhoodList[map]},
     With[
       {Xnei = X[[#]]& /@ nei},
       MapThread[FaceAxesCompiled2D, {X, Xnei, MapThread[TriangleNormalsCompiled, {X, Xnei}]}]]],
   (* #FaceCoordinates *)
   FaceCoordinates[map] :> With[
     {X = VertexCoordinates[map],
      FF = FaceList[map]},
     X[[#]]& /@ FF],
   (* #FaceRelativeCoordinates *)
   FaceRelativeCoordinates[map] :> With[
     {X = VertexCoordinates[map],
      nei = NeighborhoodList[map],
      Fn = FaceNormals[map]},
     With[
       {Xnei = X[[#]]& /@ nei},
       MapThread[
         TriangleCoordinatesCompiled,
         {X, Xnei, MapThread[FaceAxesCompiled2D, {X, Xnei, Fn}]}]]],

   (* #EdgeList *)
   EdgeList[map] :> Map[Apply[UndirectedEdge, #]&, EdgePairs[map]],
   EdgeList[map, patt_] := Count[EdgeList[map], patt, {1}],
   (* #EdgeLengths *)
   EdgeLengths[map] :> With[
     {EL = Transpose[EdgePairs[map]],
      X = VertexCoordinates[map]},
     With[
       {X1 = Transpose[X[[EL[[1]]]]],
        X2 = Transpose[X[[EL[[2]]]]]},
       Sqrt[Total[(X1 - X2)^2]]]],
   EdgeWeight[map] := EdgeLengths[map],
   EdgeWeight[map, e:(List|UndirectedEdge)[_Integer, _Integer]] := Part[
     EdgeLengths[map],
     EdgeIndex[map, e]],
   EdgeWeight[map, es:{(List|UndirectedEdge)[_Integer,_Integer]..}] := Part[
     EdgeLengths[map],
     EdgeIndex[map, #]& /@ es],
   (* #EdgeCoordinates *)
   EdgeCoordinates[map] :> With[{X = VertexCoordinates[map]}, X[[#]]& /@ EdgePairs[map]],
   
   (* #NeighborhoodList *)
   NeighborhoodList[map] :> With[
     {X = VertexCoordinates[map],
      E = EdgePairs[map]},
     Last[
       Reap[
         Scan[
           Function[{pair},
             Sow[pair[[1]], pair[[2]]];
             Sow[pair[[2]], pair[[1]]]],
           E],
         Range[Length[X]],
         Function[{id, neighbors},
           Sequence @@ With[
             {neis = Union[neighbors]},
             With[
               {U = X[[neis]]},
               SortBy[Thread[neis -> U], ArcTan[#[[2,1]], #[[2,2]]]&][[All,1]]]]]]]],
   (* #NeighborhoodAngles *)
   NeighborhoodAngles[map] :> With[
     {X = VertexCoordinates[map],
      nei = NeighborhoodList[map]},
     MapThread[NeighborhoodAnglesCompiled2D, {X, X[[#]]& /@ nei}]],
   NeighborhoodAngles[map, X_] := MapThread[
     NeighborhoodAnglesCompiled2D, 
     {X, X[[#]]& /@ NeighborhoodList[map]}],
   (* #NeighborhoodBisectors *)
   NeighborhoodBisectors[map] :> With[
     {X = VertexCoordinates[map],
      nei = NeighborhoodList[map]},
     MapThread[NeighborhoodBisectorsCompiled, {X, X[[#]] & /@ nei}]],
   NeighborhoodBisectors[map, X_] := MapThread[
     NeighborhoodBisectorsCompiled, 
     {X, X[[#]]& /@ NeighborhoodList[map]}],
   (* #NeighborhoodEdgeLengths *)
   NeighborhoodEdgeLengths[map] :> With[
     {X = VertexCoordinates[map],
      nei = NeighborhoodList[map]},
     MapThread[NeighborhoodEdgeLengthsCompiled, {X, X[[#]]& /@ nei}]],
   NeighborhoodEdgeLengths[map, X_] := MapThread[
     NeighborhoodEdgeLengthsCompiled, 
     {X, X[[#]]& /@ NeighborhoodList[map]}],

   (* #MeshRegion *)
   MeshRegion[map] :> MeshRegion[VertexCoordinates[map], Polygon[FaceList[map]]],

   (* #BoundaryMeshRegion *)
   BoundaryMeshRegion[map] :> BoundaryMeshRegion[VertexCoordinates[map], Polygon[FaceList[map]]],
   TriangulateMesh[map, opts___] := TriangulateMesh[BoundaryMeshRegion[map], opts],
   HighlightMesh[map, opts___] := HighlightMesh[BoundaryMeshRegion[map], opts],
   DimensionalMeshComponents[map, opts___] := DimensionalMeshComponents[
     BoundaryMeshRegion[map], 
     opts],
   ConnectedMeshComponents[map, opts___] := ConnectedMeshComponents[
     BoundaryMeshRegion[map], 
     opts],
   MeshCoordinates[map, opts___] := MeshCoordinates[BoundaryMeshRegion[map], opts],
   MeshCells[map, opts___] := MeshCells[BoundaryMeshRegion[map], opts],
   MeshPrimitives[map, opts___] := MeshPrimitives[BoundaryMeshRegion[map], opts],
   MeshCellIndex[map, opts___] := MeshCellIndex[BoundaryMeshRegion[map], opts],
   MeshCellCount[map, opts___] := MeshCellCount[BoundaryMeshRegion[map], opts],
   RegionDimension[map, opts___] := RegionDimension[BoundaryMeshRegion[map], opts],
   RegionEmbeddingDimension[map, opts___] := RegionEmbeddingDimension[
     BoundaryMeshRegion[map],
     opts],
   RegionMeasure[map, opts___] := RegionMeasure[BoundaryMeshRegion[map], opts],
   RegionNearest[map, opts___] := RegionNearest[BoundaryMeshRegion[map], opts],
   RegionDistance[map, opts___] := RegionDistance[BoundaryMeshRegion[map], opts],
   SignedRegionDistance[map, opts___] := SignedRegionDistance[BoundaryMeshRegion[map], opts],
   RegionMember[map, opts___] := RegionMember[BoundaryMeshRegion[map], opts],
   RegionBounds[map, opts___] := RegionBounds[BoundaryMeshRegion[map], opts],
   RegionBoundary[map, opts___] := RegionBoundary[BoundaryMeshRegion[map], opts],
   RegionCentroid[map, opts___] := RegionCentroid[BoundaryMeshRegion[map], opts],
   Volume[map, opts___] := Volume[BoundaryMeshRegion[map], opts],

   (* #Graph *)
   Graph[map] :> Graph[
     VertexList[map],
     EdgeList[map],
     VertexCoordinates -> VertexCoordinates[map],
     EdgeWeight -> EdgeLengths[map]],
   BetweennessCentrality[map, opts___]       := BetweennessCentrality[Graph[map], opts],
   ClosenessCentrality[map, opts___]         := ClosenessCentrality[Graph[map], opts],
   DegreeCentrality[map, opts___]            := DegreeCentrality[Graph[map], opts],
   EdgeBetweennessCentrality[map, opts___]   := EdgeBetweennessCentrality[Graph[map], opts],
   EdgeConnectivity[map, opts___]            := EdgeConnectivity[Graph[map], opts],
   EdgeCycleMatrix[map, opts___]             := EdgeCycleMatrix[Graph[map], opts],
   EigenvectorCentrality[map, opts___]       := EigenvectorCentrality[Graph[map], opts],
   EulerianGraphQ[map]                       := EulerianGraphQ[Graph[map]],
   FindCycle[map, opts___]                   := FindCycle[Graph[map], opts],
   FindEdgeIndendentPaths[map, opts___]      := FindEdgeIndendentPaths[Graph[map], opts],
   FindEulerianCycle[map, opts___]           := FindEulerianCycle[Graph[map], opts],
   FindFundamentalCycles[map, opts___]       := FindFundamentalCycles[Graph[map], opts],
   FindHamiltonianCycle[map, opts___]        := FindHamiltonianCycle[Graph[map], opts],
   FindMaximumFlow[map, opts___]             := FindMaximumFlow[Graph[map], opts],
   FindMinimumCostFlow[map, opts___]         := FindMinimumCostFlow[Graph[map], opts],
   FindMinimumFlow[map, opts___]             := FindMinimumFlow[Graph[map], opts],
   FindPath[map, opts___]                    := FindPath[Graph[map], opts],
   FindPostmanTour[map, opts___]             := FindPostmanTour[Graph[map], opts],
   FindShortestPath[map, opts___]            := FindShortestPath[Graph[map], opts],
   FindShortestTour[map, opts___]            := FindShortestTour[Graph[map], opts],
   FindVertexIndendentPaths[map, opts___]    := FindVertexIndendentPaths[Graph[map], opts],
   GraphAssortativity[map, opts___]          := GraphAssortativity[Graph[map], opts],
   GraphCenter[map, opts___]                 := GraphCenter[Graph[map], opts],
   GraphDiameter[map, opts___]               := GraphDiameter[Graph[map], opts],
   GraphDistance[map, opts___]               := GraphDistance[Graph[map], opts],
   GraphDistanceMatrix[map, opts___]         := GraphDistanceMatrix[Graph[map], opts],
   GraphPeriphery[map, opts___]              := GraphPeriphery[Graph[map], opts],
   GraphPower[map, opts___]                  := GraphPower[Graph[map], opts],
   GraphRadius[map, opts___]                 := GraphRadius[Graph[map], opts],
   GraphReciprocity[map, opts___]            := GraphReciprocity[Graph[map], opts],
   GlobalClusteringCoefficient[map, opts___] := GlobalClusteringCoefficient[map, opts],
   HamiltonianGraphQ[map]                    := HamiltonianGraphQ[Graph[map]],
   HITSCentrality[map, opts___]              := HITSCentrality[Graph[map], opts],
   KatzCentrality[map, opts___]              := KatzCentrality[Graph[map], opts],
   LocalClusteringCoefficient[map, opts___]  := LocalClusteringCoefficient[Graph[map], opts],
   MeanClusteringCoefficient[map, opts___]   := MeanClusteringCoefficient[Graph[map], opts],
   MeanDegreeConnectivity[map, opts___]      := MeanDegreeConnectivity[Graph[map], opts],
   MeanNeighborDegree[map, opts___]          := MeanNeighborDegree[Graph[map], opts],
   PageRankCentrality[map, opts___]          := PageRankCentrality[Graph[map], opts],
   RadialityCentrality[map, opts___]         := RadialityCentrality[Graph[map], opts],
   ShortestPathFunction[map, opts___]        := ShortestPathFunction[Graph[map], opts],
   StasusCentrality[map, opts___]            := StasusCentrality[Graph[map], opts],
   VertexConnectivity[map, opts___]          := VertexConnectivity[Graph[map], opts],
   VertexCorrelationSimilarity[map, opts___] := VertexCorrelationSimilarity[Graph[map], opts],
   VertexCosineSimilarity[map, opts___]      := VertexCosineSimilarity[Graph[map], opts],
   VertexDiceSimilarity[map, opts___]        := VertexDiceSimilarity[Graph[map], opts],
   VertexEccentricity[map, opts___]          := VertexEccentricity[Graph[map], opts],
   VertexInDegree[map, opts___]              := VertexDegree[map],
   VertexJaccardSimilarity[map, opts___]     := VertexJaccardSimilarity[Graph[map], opts],
   VertexOutDegree[map, opts___]             := VertexDegree[map],

   (* #CorticalMapQ *)
   CorticalMapQ[map] -> Which[
     And[
       Options[map] =!= Automatic,
       Complement[Options[map][[All,1]], $CortexPlotOptions[[All, 1]]] != {}], Message[
         CorticalMap::badarg,
         "Unrecognized option given to CorticalMap"],
     True, True]},
  SetSafe -> True,
  Symbol -> CorticalMesh2D];

(* #CorticalMapQ *)
CorticalMapQ[_] := False;

(* #CorticalObjectQ *)
CorticalObjectQ[c_] := Or[CorticalMeshQ[c], CorticalMapQ[c]];

(* We want to make sure to have a nicely formatted output for a mesh or projection *)

MakeBoxes[mesh_CorticalMesh3D, form_] := RowBox[
  {"CorticalMesh3D","[",
   "<"<>ToString[VertexCount[mesh]]<>" vertices>", ",",
   "<"<>ToString[EdgeCount[mesh]]<>" edges>", ",",
   "<"<>ToString[FaceCount[mesh]]<>" faces>","]"}];
MakeBoxes[mesh_CorticalMesh2D, form_] := RowBox[
  {"CorticalMesh2D","[",
   "<"<>ToString[VertexCount[mesh]]<>" vertices>", ",",
   "<"<>ToString[EdgeCount[mesh]]<>" edges>", ",",
   "<"<>ToString[FaceCount[mesh]]<>" faces>","]"}];

(* Protect these functions... *)
Protect[CorticalMesh, CorticalMeshQ, InitialVertexCoordinates, InitialVertexList, InitialFaceList,
        InitialEdgePairs, VertexCoordinates, EdgePairs, FaceIndexArray, EdgeIndexArray, FaceIndex,
        EdgeIndex, EdgeCoordinates, OptionalProperties, VertexProperties, EdgeProperties,
        FaceProperties, FaceAxes, FaceAngles, FaceCoordinates, FaceNormals, FaceRelativeCoordinates,
        CornerList, EdgeLengths, NeighborhoodList, NeighborhoodAngles, NeighborhoodBisectors,
        NeighborhoodEdgeLengths, SourceImage, VertexEdgeList, VertexFaceList, EdgeFaceList];


(**************************************************************************************************)
(* Properties for Cortical Objects *)

(* Here we do the property stuff related to meshes; these need to be overrided *)
Unprotect[PropertyValue, SetProperty, RemoveProperty, PropertyList];

(* PropertyValue extractors *)
PropertyValue[mesh_?CorticalObjectQ, prop:Except[_List]] := With[
  {list = Replace[
     prop,
     Join[
       Properties[mesh],
       {EdgeWeight :> {
          VertexList -> $Failed,
          EdgeList -> EdgeLengths[mesh],
          FaceList -> $Failed},
        VertexCoordinates :> {
          VertexList -> VertexCoordinates[mesh],
          EdgeList -> $Failed,
          FaceList -> $Failed}}]]},
  If[list === prop || list === $Failed,
    $Failed,
    list]];
PropertyValue[{mesh_?CorticalObjectQ, VertexList}, prop:Except[_List]] := With[
  {list = Replace[prop, Properties[mesh]]},
  If[list === prop || list === $Failed,
    If[prop === VertexCoordinates, VertexCoordinates[mesh], $Failed],
    list[[1,2]]]];
PropertyValue[{mesh_?CorticalObjectQ, EdgeList}, prop:Except[_List]] := With[
  {list = Replace[prop, Properties[mesh]]},
  If[list === prop || list === $Failed,
    If[prop === EdgeWeight, EdgeLengths[mesh], $Failed],
    list[[2,2]]]];
PropertyValue[{mesh_?CorticalObjectQ, FaceList}, prop:Except[_List]] := With[
  {list = Replace[prop, Properties[mesh]]},
  If[list === prop || list === $Failed,
    $Failed,
    list[[3,2]]]];
PropertyValue[{mesh_?CorticalObjectQ, i_Integer}, prop:Except[_List]] := With[
  {list = Replace[prop, Properties[mesh]]},
  If[list === prop || list === $Failed || list[[1]] === $Failed,
    If[prop === VertexCoordinates, VertexCoordinates[mesh][[i]], $Failed],
    list[[1, 2]][[i]]]];
PropertyValue[{mesh_?CorticalObjectQ, e:(List|UndirectedEdge)[_Integer, _Integer]}, 
              prop:Except[_List]] := With[
  {list = Replace[prop, Properties[mesh]]},
  If[list === prop || list === $Failed || list[[2]] === $Failed,
    If[prop === EdgeWeight, EdgeLengths[mesh][[EdgeIndex[mesh, e]]], $Failed],
    list[[2, 2]][[EdgeIndex[mesh, e]]]]];
PropertyValue[{mesh_?CorticalObjectQ, f:{_Integer, _Integer, _Integer}}, prop:Except[_List]] := With[
  {list = Replace[prop, Properties[mesh]]},
  If[list === prop || list === $Failed || list[[3]] === $Failed,
    $Failed,
    list[[3, 2]][[FaceIndex[mesh, f]]]]];
PropertyValue[mesh_?CorticalObjectQ, prop_List] := Map[PropertyValue[mesh, #]&, prop];
PropertyValue[{mesh_?CorticalObjectQ, x_}, prop_List] := Map[PropertyValue[{mesh, x}, #]&, prop];

(* PropertyList stuff *)
PropertyList[{mesh_?CorticalObjectQ, type:(VertexList|EdgeList|FaceList)}] := With[
  {props = Properties[mesh],
   idx = Replace[type, {VertexList -> 1, EdgeList -> 2, FaceList ->3}]},
  Join[
    Select[props, (Hold @@ #[[2,idx]])[[{2}]] =!= Hold[$Failed]&][[All, 1]],
    Which[type === VertexList, {VertexCoordinates}, type === EdgeList, {EdgeWeight}, True, {}]]];
PropertyList[{mesh_?CorticalObjectQ, i_Integer}] := With[
  {props = Properties[mesh]},
  Append[
    Select[props, And[ListQ[#[[2,1,2]]], #[[2,1,2,i]] =!= $Failed]&][[All, 1]],
    VertexCoordinates]];
PropertyList[{mesh_?CorticalObjectQ, e:(List|UndirectedEdge)[_Integer, _Integer]}] := With[
  {props = Properties[mesh],
   idx = EdgeIndex[mesh, e]},
  Append[
    Select[props, And[ListQ[#[[2,2,2]]], #[[2,2,2,idx]] =!= $Failed]&][[All, 1]],
    EdgeWeight]];
PropertyList[{mesh_?CorticalObjectQ, f:{_Integer, _Integer, _Integer}}] := With[
  {props = Properties[mesh],
   idx = FaceIndex[mesh, f]},
  Select[props, And[ListQ[#[[2,3,2]]], #[[2,3,2,idx]] =!= $Failed]&][[All, 1]]];

(* SetProperty stuff... *)

SetProperty[{mesh_?CorticalObjectQ, type:(VertexList|EdgeList|FaceList)}, prop_ -> vals_List] := With[
  {allList = Switch[type, 
     VertexList, VertexProperties[mesh],
     EdgeList, EdgeProperties[mesh],
     FaceList, FaceProperties[mesh]],
   propType = Switch[type, 
     VertexList, VertexProperties, EdgeList, EdgeProperties, FaceList, FaceProperties]},
  With[
    {list = Replace[prop, allList]},
    Which[
      Length[vals] != Length[type[mesh]], $Failed,
      (* If the property doesn't yet exist, add it. *)
      list === prop || list === $Failed, Which[
        prop === VertexCoordinates, If[type === VertexList,
          Clone[mesh, VertexCoordinates -> vals],
          $Failed],
        prop === EdgeWeight, $Failed,
        True, Clone[mesh, propType -> Append[allList, prop -> vals]]],
      (* If the property already exists, we overwrite it *)
      True, Clone[
        mesh,
        propType -> Replace[list, (Rule|RuleDelayed)[prop, _] -> (prop -> vals), {1}]]]]];
(* for delayed rules, we handle things slightly differently *)
SetProperty[{mesh_?CorticalObjectQ, type:(VertexList|EdgeList|FaceList)}, prop_ :> vals_] := With[
  {allList = Switch[type, 
     VertexList, VertexProperties[mesh],
     EdgeList, EdgeProperties[mesh],
     FaceList, FaceProperties[mesh]],
   propType = Switch[type, 
     VertexList, VertexProperties, EdgeList, EdgeProperties, FaceList, FaceProperties]},
  With[
    {list = Replace[prop, allList],
     sym = Unique["delayedProperty"]},
    sym := With[
      {res = vals},
      If[!ListQ[res] || Length[res] != Length[type[mesh]],
        Table[$Failed, {Length[type[mesh]]}],
        (sym = res)]];
    If[list === prop || list === $Failed, 
      (* If the property doesn't yet exist, add it. *)
      Which[
        prop === VertexCoordinates, If[type === VertexList,
          Clone[mesh, VertexCoordinates -> vals],
          $Failed],
        prop === EdgeWeight, $Failed,
        True, Clone[mesh, propType -> Append[allList, (prop :> sym)]]],
      (* If the property already exists, we overwrite it *)
      Clone[mesh, propType -> Replace[list, (Rule|RuleDelayed)[prop, _] -> (prop :> sym), {1}]]]]];

(* bulk property setting *)
SetProperty[{mesh_?CorticalObjectQ,
             type:(VertexList|EdgeList|FaceList)}, rs:{(_Rule|_RuleDelayed)..}] := Fold[
  SetProperty[{#1, type}, #2]&,
  mesh,
  rs];
(* confused property setting (note that :> is not valid here) *)
SetProperty[mesh_?CorticalObjectQ, rs:{_Rule..}] := Fold[SetProperty, mesh, rs];
SetProperty[mesh_?CorticalObjectQ, r:(_ -> vals_List)] := Switch[
  Length[vals],
  VertexCount[mesh], SetProperty[{mesh, VertexList}, r],
  EdgeCount[mesh],   SetProperty[{mesh, EdgeList}, r],
  FaceCount[mesh],   SetProperty[{mesh, FaceList}, r],
  _, $Failed];

(* Set individual vertices/edges/faces *)
SetProperty[{mesh_?CorticalObjectQ, obj_}, prop_ -> val_] := With[
  {type = Switch[obj,
     _Integer,                                  VertexProperties -> obj,
     (List|UndirectedEdge)[_Integer, _Integer], EdgeProperties -> EdgeIndex[mesh, obj],
     {_Integer, _Integer, _Integer},            FaceProperties -> FaceIndex[mesh, obj],
     _, $Failed]},
  If[type === $Failed,
    $Failed,
    Which[
      prop === VertexCoordinates, If[type[[1]] =!= VertexProperties, 
        $Failed,
        Clone[mesh, VertexCoordinates -> ReplacePart[VertexCoordinates[mesh], obj -> val]]],
      prop === EdgeWeight, $Failed,
      True, With[
        {list = Replace[prop, type[[1]][mesh]]},
        Clone[
          mesh,
          type[[1]] -> If[list === prop,
            With[
              {n = Switch[type[[1]], 
                 VertexProperties, VertexCount[mesh],
                 EdgeProperties, EdgeCount[mesh],
                 FaceProperties, FaceCount[mesh]]},
              Append[
                type[[1]][mesh],
                prop -> ReplacePart[ConstantArray[$Failed, n], type[[2]] -> val]]],
            Replace[
              type[[1]][mesh],
              (prop -> v_) :> (prop -> ReplacePart[v, type[[2]] -> val]),
              {1}]]]]]]];
SetProperty[{mesh_?CorticalObjectQ, obj_}, prop_ :> val_] := With[
  {type = Switch[obj,
     _Integer,                                  VertexProperties -> obj,
     (List|UndirectedEdge)[_Integer, _Integer], EdgeProperties -> EdgeIndex[mesh, obj],
     {_Integer, _Integer, _Integer},            FaceProperties -> FaceIndex[mesh, obj],
     _, $Failed]},
  If[type === $Failed,
    $Failed,
    Which[
      prop === VertexCoordinates, If[type[[1]] =!= VertexProperties, 
        $Failed,
        Clone[mesh, VertexCoordinates -> ReplacePart[VertexCoordinates[mesh], obj :> val]]],
      prop === EdgeWeight, $Failed,
      True, With[
        {list = Replace[prop, type[[1]][mesh]]},
        Clone[
          mesh,
          type[[1]] -> If[list === prop,
            With[
              {n = Switch[type[[1]], 
                 VertexProperties, VertexCount[mesh],
                 EdgeProperties, EdgeCount[mesh],
                 FaceProperties, FaceCount[mesh]]},
              Append[
                type[[1]][mesh],
                prop -> ReplacePart[ConstantArray[$Failed, n], type[[2]] :> val]]],
            Replace[
              type[[1]][mesh],
              (prop -> v_) :> (prop -> ReplacePart[v, type[[2]] :> val]),
              {1}]]]]]]];
SetProperty[{mesh_?CorticalObjectQ, obj_}, rs:{(_Rule|_RuleDelayed)..}] := Fold[
  SetProperty[{#1, obj}, #2]&,
  mesh,
  rs];

(* Remove Property stuff *)
RemoveProperty[mesh_?CorticalObjectQ, prop:Except[_List]] := With[
  {list = Replace[prop, Properties[mesh]]},
  If[list === prop, 
    mesh,
    Clone[
      mesh,
      VertexProperties -> DeleteCases[VertexProperties[mesh], Rule[prop, _]],
      EdgeProperties -> DeleteCases[EdgeProperties[mesh], Rule[prop, _]],
      FaceProperties -> DeleteCases[FaceProperties[mesh], Rule[prop, _]]]]];
RemoveProperty[mesh_?CorticalObjectQ, prop:{_Rule..}] := Fold[RemoveProperty, mesh, prop];
RemoveProperty[{mesh_?CorticalObjectQ, t:(VertexList|EdgeList|FaceList)}, prop:Except[_List]] := With[
  {type = Switch[t,
     VertexList, VertexProperties,
     EdgeList,   EdgeProperties,
     FaceList,   FaceProperties,
     _, $Failed]},
  With[
    {list = Replace[prop, type[mesh]]},
    If[list === $Failed || list === prop,
      mesh,
      Clone[mesh, type -> DeleteCases[type[mesh], prop -> _]]]]];
RemoveProperty[{mesh_?CorticalObjectQ, obj_}, prop:Except[_List]] := With[
  {type = Switch[
     _Integer,                                  VertexProperties -> obj,
     (List|UndirectedEdge)[_Integer, _Integer], EdgeProperties -> EdgeIndex[mesh, obj],
     {_Integer, _Integer, _Integer},            FaceProperties -> FaceIndex[mesh, obj],
     _, $Failed]},
  With[
    {list = Replace[prop, type[[1]][mesh]]},
    If[list === prop || list === $Failed,
      mesh,
      Clone[
        mesh,
        type[[1]] -> Replace[
          type[[1]][mesh],
          (prop -> list) :> (prop -> ReplacePart[list, type[[2]] -> $Failed]),
          {1}]]]]];
RemoveProperty[{mesh_?CorticalObjectQ, obj_}, prop_List] := Fold[
  RemoveProperty[{#1, obj}, #2]&,
  mesh,
  prop];

Protect[PropertyValue, SetProperty, RemoveProperty, PropertyList];


(* #CorticalCurvatureColor ************************************************************************)
CorticalCurvatureColor[c_] := If[c > 0, GrayLevel[0.2], Gray];
SetAttributes[CorticalCurvatureColor, Listable];
Protect[CorticalCurvatureColor];

(* #CorticalCurvatureVertexColors *****************************************************************)
CorticalCurvatureVertexColors[mesh_?CorticalObjectQ] := With[
  {curv = Replace[
     PropertyValue[{mesh, VertexList}, Curvature],
     $Failed :> PropertyValue[{mesh, VertexList}, "Curvature"]]},
  If[curv === $Failed,
    ConstantArray[Gray, VertexCount[mesh]],
    CorticalCurvatureColor /@ curv]];
Protect[CorticalCurvatureVertexColors];


(* #CortexSphericalQ ******************************************************************************)
CortexSphericalQ[c_?CorticalMeshQ] := With[
  {centroid = Mean[VertexCoordinates[c]]},
  With[
    {distances = Sqrt[Total[MapThread[Subtract, {Transpose[VertexCoordinates[c]], centroid}]^2]]},
    With[
      {mu = Mean[distances], sigma = StandardDeviation[distances]},
      sigma / mu < 0.01]]];
(* #CortexToSphere ********************************************************************************)
CortexToSphere[c_?CorticalMeshQ /; CortexSphericalQ[c]] := c;
CortexToSphere[c_?CorticalMeshQ /; Not[CortexSphericalQ[c]]] := With[
  {X0 = With[{tmp = VertexCoordinates[c]}, MapThread[Subtract, {Transpose[tmp], Mean[tmp]}]],
   EL0 = NeighborhoodEdgeLengths[c],
   nei = NeighborhoodList[c]},
  With[
    {mu = Mean[Sqrt@Total[X0^2]]},
    Block[
      {f, df, X},
      f[X_List] := Plus[
        Total[(mu - Sqrt@Total[Transpose[X]^2])^2],
        5.0 * Total[Flatten[EL0 - NeighborhoodEdgeLengths[c, X]]^2]];
      df[X_List] := With[
        {norms = Sqrt@Total[Transpose[X]^2]},
        Flatten@Plus[
          -2.0 * (mu - norms) / norms * X,
          MapThread[
            Function[{x0, el0, xnei},
              With[
                {dx = MapThread[Subtract, {xnei, x0}]},
                With[
                  {neiDists = Sqrt@Total[dx^2]},
                  -4.0 * Dot[(el0 - neiDists)/neiDists, Transpose[dx]]]]],
            {X, EL0, Transpose[X[[#]]]& /@ nei}]]];
      With[
        {res = First@Quiet[
           FindArgMin[
             f[X],
             {X, Transpose@X0},
             Gradient :> df[X]],
           {FindArgMin::cvmit, FindArgMin::lstol, FindArgMin::sdprec}]},
        Clone[c, VertexCoordinates -> res]]]]];


(* #CortexPlot3D **********************************************************************************)

(* We have a few default options for color schema; they get picked in this order *)
$CortexPlotDefaultColorSchemas = {
  Curvature -> CorticalCurvatureVertexColors,
  "Curvature" -> CorticalCurvatureVertexColors,
  "Data" -> Function[{mesh},
    With[
      {data = PropertyValue[{mesh, VertexList}, "Data"]},
      With[
        {range = (#[[2]] + 3 * {-(#[[2]] - #[[1]]), (#[[2]] - #[[1]])})& @ Quartiles[data]},
        Map[
          Blend[{Blue, Darker[Cyan, 1/6], Darker[Green, 1/3], Darker[Yellow, 1/6], Red}, #]&,
          (data - range[[1]]) / (range[[2]] - range[[1]])]]]]};
Protect[$CortexPlotDefaultColorSchemas];

Options[CortexPlot3D] = Map[(#[[1]] -> Automatic)&, $CortexPlot3DOptions];
CortexPlot3D[mesh_?CorticalMeshQ, opts:OptionsPattern[]] := With[
  {Opt = Function[{name},
     Fold[
       Replace,
       OptionValue[name],
       {Automatic :> Replace[name, Options[mesh]],
        (name|Automatic) :> Replace[name, $CortexPlot3DOptions]}]],
   GetProperties = Function[{type},
     Replace[
       Thread[Rule[#, PropertyValue[{mesh, type}, #]]]& /@ PropertyList[{mesh, type}],
       {{} :> Table[{}, {Length[type[mesh]]}],
        l_ :> Transpose[l]}]],
   U = VertexList[mesh],
   X = VertexCoordinates[mesh],
   F = FaceList[mesh],
   vnorms = Replace[
     PropertyValue[{mesh, VertexList}, VertexNormals],
     $Failed :> Replace[
       PropertyValue[{mesh, VertexList}, "VertexNormals"],
       $Failed :> VertexNormals[mesh]]]},
  With[
    {vcolors = Replace[
       Opt[ColorFunction],
       {Automatic :> With[
          {sel = SelectFirst[
             $CortexPlotDefaultColorSchemas,
             ArrayQ[PropertyValue[{mesh, VertexList}, #[[1]]], 1, NumericQ]&,
             None]},
          If[Head[sel] === None, 
            ConstantArray[Gray, VertexCount[mesh]],
            (sel[[2]])[mesh]]],
        None -> None,
        f_ :> With[
          {known = Replace[f, $CortexPlotDefaultColorSchemas]},
          If[f === known,
            MapThread[f, {X, U, GetProperties[VertexList]}],
            known[mesh]]]}]},
    With[
      {vfn = Opt[VertexRenderingFunction],
       efn = Opt[EdgeRenderingFunction],
       ffn = Opt[FaceRenderingFunction],
       vprop = If[Or[(vfn =!= Automatic && vfn =!= None),
                     (efn =!= Automatic && efn =!= None),
                     (ffn =!= Automatic && ffn =!= None)],
         GetProperties[VertexList],
         None]},
      Graphics3D[
        GraphicsComplex[
          VertexCoordinates[mesh],
          {If[ffn =!= Automatic,
             {EdgeForm[], Gray, 
              MapThread[
                ffn,
                {FaceCoordinates[mesh], F, GetProperties[FaceList], vprop[[#]]& /@ F}]},
             {EdgeForm[], Gray, Polygon[F]}],
           If[efn === None || efn === Automatic, 
             {},
             MapThread[
               efn,
               {EdgeCoordinates[mesh], EdgePairs[mesh],
                GetProperties[EdgeList], vprop[[#]]& /@ EdgePairs[mesh]}]],
           If[vfn === None || vfn === Automatic,
             {},
             MapThread[vfn, {X, U, vprop}]]},
          VertexColors -> If[(ffn === Automatic || ffn === None) && vcolors =!= None, 
            vcolors,
            None],
          VertexNormals -> vnorms],
        Sequence@@FilterRules[
          Join[{opts}, Options[mesh], $CortexPlot3DOptions],
          Options[Graphics3D][[All,1]]]]]]];
Protect[CortexPlot3D];


(* #CortexPlot ************************************************************************************)
Options[CortexPlot] = Map[(#[[1]] -> Automatic)&, Options[CorticalMap]];
CortexPlot[mesh_?CorticalMapQ, opts:OptionsPattern[]] := With[
  {Opt = Function[{name},
     Fold[
       Replace,
       OptionValue[name],
       {Automatic :> Replace[name, Options[mesh]],
        (name|Automatic) :> Replace[name, $CortexPlotOptions]}]],
   GetProperties = Function[{type},
     Replace[
       Thread[Rule[#, PropertyValue[{mesh, type}, #]]]& /@ PropertyList[{mesh, type}],
       {{} :> Table[{}, {Length[type[mesh]]}],
        l_ :> Transpose[l]}]],
   U = VertexList[mesh],
   X = VertexCoordinates[mesh],
   F = Partition[Normal[VertexIndexArray[mesh][[Flatten@FaceList[mesh]]]], 3]},
  With[
    {vcolors = Replace[
       Opt[ColorFunction],
       {Automatic :> With[
          {sel = SelectFirst[
             $CortexPlotDefaultColorSchemas,
             ArrayQ[PropertyValue[{mesh, VertexList}, #[[1]]], 1, NumericQ]&,
             None]},
          If[Head[sel] === None, 
            ConstantArray[Gray, VertexCount[mesh]],
            (sel[[2]])[mesh]]],
        None -> None,
        f_ :> With[
          {known = Replace[f, $CortexPlotDefaultColorSchemas]},
          If[f === known,
            MapThread[f, {X, U, GetProperties[VertexList]}],
            known[mesh]]]}]},
    With[
      {vfn = Opt[VertexRenderingFunction],
       efn = Opt[EdgeRenderingFunction],
       ffn = Opt[FaceRenderingFunction],
       vprop = If[Or[(vfn =!= Automatic && vfn =!= None),
                     (efn =!= Automatic && efn =!= None),
                     (ffn =!= Automatic && ffn =!= None)],
         GetProperties[VertexList],
         None]},
      Graphics[
        GraphicsComplex[
          VertexCoordinates[mesh],
          {If[ffn =!= Automatic,
             {EdgeForm[], Gray, 
              MapThread[
                ffn,
                {FaceCoordinates[mesh], F, GetProperties[FaceList], vprop[[#]]& /@ F}]},
             {EdgeForm[], Gray, Polygon[F]}],
           If[efn === None || efn === Automatic, 
             {},
             MapThread[
               efn,
               {EdgeCoordinates[mesh], EdgePairs[mesh],
                GetProperties[EdgeList], vprop[[#]]& /@ EdgePairs[mesh]}]],
           If[vfn === None || vfn === Automatic,
             {},
             MapThread[vfn, {X, U, vprop}]]},
          VertexColors -> If[(ffn === Automatic || ffn === None) && vcolors =!= None, 
            vcolors,
            None]],
        Sequence@@FilterRules[
          Join[{opts}, Options[mesh], $CortexPlotOptions],
          Options[Graphics][[All,1]]]]]]];
Protect[CortexPlot];


End[];
EndPackage[];
