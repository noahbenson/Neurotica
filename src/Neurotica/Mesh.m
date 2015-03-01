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
A cortical mesh resembles both a graph object and a boundary mesh region object. They are immutable data structures defined using the DefineImmutable function found in the Neurotica`Util` namespace; accordingly, they aren't designed to be edited, but to be copied (while reusing and sharing data between copies efficiently). To edit the options or data in a CorticalMesh m, call CorticalMesh[m, <option> -> <value>...]. Valid options for such a call are any option that may be passed to CorticalMesh normally, FaceList, and VertexCoordinates. Cortical meshes support the following queries:
  * Properties, PropertyValue, SetProperty, RemoveProperty: All property functions supported by Graphs are also supported by cortical meshes. It is recommended that all data attached to nodes, faces, or edges be attached using properties. Note that the Properties option is accepted by CorticalMesh (see ?Properties) and Property wrappers in the vertex and face lists are parsed as well (see ?Property).
  * VertexList, VertexCount, EdgeList, EdgeCount, EdgePairs: Most graph functions can be used with a cortical mesh as if the cortical mesh were a graph. In addition, the function EdgePairs[mesh] yields the same result as EdgeList[mesh] but with lists if neighboring vertices rather than UndirectedEdge forms.
  * VertexCoordinates[mesh] additionally yields the list of vertex coordinates for the given mesh.
  * VertexNormals[mesh], FaceNormals[mesh] yields a list of the normal vector to each vertex or face in the mesh.
  * FaceList, FaceCount: In addition to edges and vertices, meshes also have faces that can be accessed with the FaceList and FaceCount functions.
  * EdgeLengths, EdgeWeights: EdgeLengths[mesh] and EdgeWeights[mesh] are identical; both yield the Euclidean distances between vertices in the mesh.
  * FaceAngles[mesh] yields the internal angles of each face in the same order as the vertices in FaceList[mesh].
  * FaceNormals[mesh] yields a list of one vector for each face in the mesh such that the vertex is orthogonal to the face.
  * FaceBisectors[mesh] yields a list of the bisecting vectors of each of the corners of each face in the mesh.
  * FaceAxes[mesh] yields a list of orthogonal 2D axes that can be used to project each face into two dimensions, one pair for each face.
  * FaceCoordinates[mesh] yields a list of the faces but with the 3D mesh coordinates in the place of the vertex indices.
  * FaceRelativeCoordinates[mesh] yields a list of coordinates relative to the FaceAxes[mesh], such that each face has been flattened into the 2D plane it defines.
  * NeighborhoodList[mesh] yields a list of the neighborhood of each vertex in the same order as VertexList[mesh]; the neighborhood of a vertex u is a list of the vertices that are adjacent to u in counterclockwise order around the vector normal to u.
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
  * 3D Graphics: CorticalMesh[] accepts any option that can be passed to Graphics3D; these options will be used as default values when plotting the cortical mesh. The options may be accessed via Options[mesh] and may be changed (when cloning a mesh m via CorticalMesh[m, options...]) with Options -> {new-options}; Options -> Automatic will reset the options to the defaults accepted by Graphics3D with the following differences: Lighting -> \"Neutral\", ColorFunction -> None, ColorFunctionScaling -> False, Boxed -> False.";
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

CalculateFaceAngleCosinesTr::usage = "CalculateFaceAngleCosinesTr[Ft, Xt] is a compiled function that yields a list of the cosines of the angles of each face in given in the transposed face list Ft using the given transposed coordinate list Xt. This is a low-level function used by CorticalMesh and CorticalMap that should generally not be called by the user.";
CalculateFaceAnglesTr::usage = "CalculateFaceAnglesTr[Ft, Xt] is a compiled function that yields a list of the angles of each face in given in the transposed face list Ft using the given transposed coordinate list Xt. This is a low-level function used by CorticalMesh and CorticalMap that should generally not be called by the user.";
CalculateFaceAngleCosines::usage = "CalculateFaceAngleCosines[F, X] is a compiled function that yields a list of the angles of each face in given in the face list F using the given coordinate list X. This is a low-level function used by CorticalMesh and CorticalMap that should generally not be called by the user.";
CalculateFaceAngles::usage = "CalculateFaceAngles[F, X] is a compiled function that yields a list of the angles of each face in given in the face list F using the given coordinate list X. This is a low-level function used by CorticalMesh and CorticalMap that should generally not be called by the user.";
CalculateFaceNormalsTr::usage = "CalculateFaceNormalsTr[Ft, Xt] is a compiled function that yields a list of the {x, y, z} coordinates of the normal vectors to each face in the transposed face list Ft using the given transposed coordinate list Xt. This is a low-level function used by CorticalMesh and CorticalMap that should generally not be called by the user.";
CalculateFaceNormals::usage = "CalculateFaceNormals[F, X] is a compiled function that yields a list of the normal vectors of each face given in the face list F using the given transposed coordinate list X. This is a low-level function used by CorticalMesh and CorticalMap that should generally not be called by the user.";
CalculateFaceAxesTr::usage = "CalculateFaceAxesTr[Ft, Xt] is a compiled function that yields a list of the {{x1, x2}, {y1, y2}, {z1, z2}} coordinate lists of the two axes ({x1, y1, z1} and {x2, y2, z2}) of each face in the given transposed face list Ft and the given transposed coordinate list Xt. This is a low-level function used by CorticalMesh and CorticalMap that should generally not be called by the user.";
CalculateFaceAxes::usage = "CalculateFaceAxes[F, X] is a compiled function that yields a list of the axes for each face given in the face list F using the given coordinate list X. This is a low-level function used by CorticalMesh and CorticalMap that should generally not be called by the user.";
CalculateFacePlaneCoordinatesTr::usage = "CalculateFacePlaneCoordinatesTr[Ft, Xt, At] is a compiled function that yields a list of the {x,y} 2D coordinates for each point of each face given in the transposed face list Ft using the given transposed coordinate list Xt. The result of this calculation is equivalent to Transpose[q, {3,2,1}] where q is the result of the non-transposed version of the function. This is a low-level function used by CorticalMesh and CorticalMap that should generally not be called by the user.";
CalculateFacePlaneCoordinates::usage = "CalculateFacePlaneCoordinates[F, X, A] is a compiled function that yields a list of the {x,y} 2D coordinates for each point of each face given in the face list F using the given coordinate list X. This is a low-level function used by CorticalMesh and CorticalMap that should generally not be called by the user.";
CalculateFaceBisectors::usage = "CalculateFaceBisectors[F, X] yields the bisector unit vectors for each corner of each face in the face list F given the coordinate matrix X.";
CalculateFaceBisectorsTr::usage = "CalculateFaceBisectorsTr[Ft, Xt] yields the transposed bisector unit vectors for each corner of each face; this is the equivalent of Transpose[CalculateFaceBisectors[Transpose @ Ft, Transpose @ Xt], {3,2,1}].";
CalculateVertexNormalsTr::usage = "CalculateVertexNormalsTr[Ft, Xt] yields the transposed normal vectors to each vertex given in the 3 by n matrix of coordinates in Xt where the triangle faces are given by columns of integer triples in Ft.";
CalculateVertexNormals::usage = "CalculateVertexNormals[F, X] yields the normal vector to each vertex given in the 3 by n matrix of coordinates in X where the triangle faces are given by integer triples in F.";
CorticalMeshNeighborhood::usage = "CorticalMeshNeighborhood[mesh, u0, d] yields the list of vertices in the given cortical mesh such that the edge-weighted distance from u0 to each vertex returned is less than d.";

VertexCoordinatesTr::usage = "VertexCoordinatesTr[mesh] yields the transposed list of vertex coordinates for the given mesh: Transpose[VertexCoordinates[mesh]].";
VertexNormalsTr::usage = "VertexNormalsTr[mesh] yields the transposed list of vertex normals for each vertex in the given mesh: Transpose[VertexNormals[mesh]].";
FaceListTr::usage = "FaceListTr[s] yields a transposition of the FaceList[s] for the given cortical mesh or map s.";
FaceAngleCosinesTr::usage = "FaceAngleCosinesTr[s] yields a transposition of the FaceAngleCosines[s] for the given cortical mesh or map s.";
FaceAnglesTr::usage = "FaceAnglesTr[s] yields a transposition of the FaceAngles[s] for the given cortical mesh or map s.";
FaceNormalsTr::usage = "FaceNormalsTr[s] yields a transposition of the FaceNormals[s] for the given cortical mesh or map s such that FaceNormalsTr[s] == Transpose[FaceNormals[s], {3,2,1}].";
FaceBisectorsTr::usage = "FaceBisectorsTr[s] yields the transpose of the list of the vectors which bisect each of the corners of the faces in s.";
FaceAxesTr::usage = "FaceAxesTr[s] yields a transposition of the FaceAxes[s] for the given cortical mesh or map s such that FaceAxesTr[s] == Transpose[FaceAxes[s], {3,2,1}].";
FaceCoordinatesTr::usage = "FaceCoordinatesTr[s] yields a transposition of the FaceCoordinates[s] for the given cortical mesh or map s such that FaceCoordinatesTr[s] == Transpose[FaceCoordinates[s], {3,2,1}].";
FacePlaneCoordinatesTr::usage = "FacePlaneCoordinatesTr[s] yields a transposition of the FacePlaneCoordinates[s] for the given cortical mesh or map s such that FacePlaneCoordinatesTr[s] == Transpose[FacePlaneCoordinates[s], {3,2,1}].";

FaceList::usage = "FaceList[s] yields a list of all faces in the surface or map s, if any, as lists of indices into VertexCoordinates[s].";
FaceCount::usage = "FaceCount[s] yields the count of all faces in the surface or map s.";
FaceAngleCosines::usage = "FaceAngles[s] yields a list of the cosines of the angles of each edge in the FaceList[s] where s may be a surface or a map. FaceAngles[s, X] yields the angles for s if the vertices in s had coordinates equal to those in the list X.";
FaceAngles::usage = "FaceAngles[s] yields a list of the angles (in radians) of each edge in the FaceList[s] where s may be a surface or a map. FaceAngles[s, X] yields the angles for s if the vertices in s had coordinates equal to those in the list X.";
FaceIndex::usage = "FaceIndex[s] yields a list of the indices into FaceList[s] such that the i'th element in FaceIndex[s] is the list of indices at which the i'th vertex in VertexCoordinates[s] appears in FaceList[s].";
FaceNormals::usage = "FaceNormals[s] yields a list of normal vectors to each trianglular face in the surface mesh s. The vector yielded for each face is orthogonal to the plane of the face.";
FaceAxes::usage = "FaceAxe[s] yields a list of the orthogonal axes to each of the faces in the surface mesh s.";
FaceBisectors::usage = "FaceBisectors[s] yields a list of the vectors which bisect each of the corners of the faces in s.";
FaceCoordinates::usage = "FaceCoordinates[s] yeilds a list of faces in which the coordinates insted of the vertex indices for each face are given.";
FacePlaneCoordinates::usage = "FacePlaneCoordinates[s] yeilds a list of coordinates, one for each face in the surface mesh s, that has been normalized to two dimensions.";

EdgePairsTr::usage = "EdgePairsTr[s] yields a transposition of EdgePairs[s] for the given cortical mesh or map s.";
EdgeCoordinatesTr::usage = "EdgeCoordinatesTr[s] yields a transposition of EdgeCoordinates[s] for the given cortical mesh or map s.";

EdgePairs::usage = "EdgePairs[s] yields a list of the undirected edges between vertices in the surface mesh s; unlike the EdgeList function, this function yields the edges as lists instead of undirected edges.";
EdgeLengths::usage = "EdgeLengths[s] yields a list of the lengths (Euclidean norm) of each edge in the EdgeList[s] where s may be a surface mesh object or projection.
EdgeLengths[s, X] yields a list of the lengths of the edges if the given mesh s had the coordinates given in the list X.";
EdgeLengthsTr::usage = "EdgeLengthsTr[s, Xt] is equivalent to EdgeLengths[s, Transpose[Xt]].";
EdgeCoordinates::usage = "EdgeCoordinates[s] yields a list identical to EdgePairs[s] except that the vertex ids in the list have been replaced with the coordinates of each vertex.";

VertexFaceMatrix::usage = "VertexFaceMatrix[M] yields a SparseArray matrix S such that each element S[[m*k + i,j]] of S is equal to 0 if vertex j is the k'th part of face i and 0 if not (where m is the number of faces in the mesh M).";
VertexEdgeMatrix::usage = "VertexEdgeMatrix[M] yields a SparseArray matrix S such that each element S[[m*k + i,j]] of S is equal to 0 if vertex j is the k'th part of edge i and 0 if not (where m is the number of edges in the mesh M).";
EdgeFaceMatrix::usage = "VertexFaceMatrix[M] yields a SparseArray matrix S such that each element S[[m*k + i,j]] of S is equal to 0 if edge j is the kth part of face i and 0 if not (where m is the number of faces in the mesh M).";

SumOverFacesTr::usage = "SumOverFacesTr[M, Q] yields the 3 x n result of summing over the faces given in the matrix Q whose first dimension must be equal to 3 and whose last dimension must be equal to the number of faces in cortical mesh or map M.";
SumOverEdgesTr::usage = "SumOverEdgesTr[M, Q] yields the 3 x n result of summing over the edges given in the matrix Q whose first dimension must be equal to 2 and whose last dimension must be equal to the number of edges in cortical mesh or map M.";
SumOverFaces::usage = "SumOverFaces[M, Q] yields the result of summing over the faces given in the matrix Q whose last dimension must be equal to 3 and whose first dimension must be equal to the number of faces in cortical mesh or map M.";
SumOverEdges::usage = "SumOverEdges[M, Q] yields the result of summing over the edges given in the matrix Q whose last dimension must be equal to 2 and whose first dimension must be equal to the number of faces in cortical mesh or map M.";
SumOverEdgesDirected::usage = "SumOverEdgesDirected[M, Q] yields the equivalent of SumOverEdges[M, Q] except that the values in Q are subtracted from the second element of each edge instead of added to it.";
SumOverEdgesDirectedTr::usage = "SumOverEdgesDirectedTr[M, Q] yields the equivalent of SumOverEdgesTr[M, Q] except that the values in Q are subtracted from the second element of each edge instead of added to it..";

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

Inclusions::usage = "Inclusions[map] yields a list of the {vertexIndives, edgeIndices, faceIndices} of the vertices, edges, and faces from the source mesh of the given map that are included in map.";

CortexPlot3D::usage = "CortexPlot3D[mesh] yields a 3D Graphics form for the given CorticalMesh3D mesh. All options available to Graphics3D may be passed to CortexPlot3D. Note that some of the default options for Graphics3D have been altered in CortexPlot3D, and 3D graphics options that have been attached to the mesh will be used as well. See ?CorticalMesh for more details.";

CortexPlot::usage = "CortexPlot[mesh] yields a Graphics form for the given CorticalMesh2D or CorticalMesh3D mesh. If the given mesh is a 3D mesh, then the options accepted by CorticalMap are used to create the projection (Method, Center, etc.). All options available to Graphics may be passed to CortexPlot. Note that some of the default options for Graphics have been altered in CortexPlot, and 2D graphics options that have been attached to the mesh will be used as well. See ?CorticalMap for more details.";

CorticalCurvatureColor::usage = "CorticalCurvatureColor[c] yields the appropriate color for cortical curvature c in a CortexPlot or CortexPlot3D; c may be a list or a single value.";
CorticalCurvatureVertexColors::usage = "CorticalCurvatureVertexColors[m] yields the colors for each vertex in the given mesh m according to CorticalCurvatureColor[c] for the curvature c of each vertex; if there is no curvature proeprty defined, then Gray is used for each vertex.";

Begin["`Private`"];

(***************************************************************************************************
 * Private Helper functions
 **************************************************************************************************)

(* #NeighborhoodAnglesCompiled *)
NeighborhoodAnglesCompiled2D = Compile[{{x0, _Real, 1}, {xnei, _Real, 2}},
  If[Length[xnei] == 0,
    {},
    With[
      {x = {#[[1]] - x0[[1]], #[[2]] - x0[[2]]}& @ Transpose[xnei]},
      With[
        {xu = (x / {#,#})& @ Sqrt @ Total[x^2]},
        ArcCos @ Total[xu * (RotateLeft /@ xu)]]]],
  RuntimeOptions -> {"Speed", "EvaluateSymbolically" -> False},
  Parallelization -> True];
NeighborhoodAnglesCompiled3D = Compile[{{x0, _Real, 1}, {xnei, _Real, 2}},
  If[Length[xnei] == 0,
    {},
    With[
      {x = {#[[1]] - x0[[1]], #[[2]] - x0[[2]], #[[3]] - x0[[3]]}& @ Transpose[xnei]},
      With[
        {xu = (x / {#,#,#})& @ Sqrt @ Total[x^2]},
        ArcCos @ Total[xu * (RotateLeft /@ xu)]]]],
  RuntimeOptions -> {"Speed", "EvaluateSymbolically" -> False},
  Parallelization -> True];
Protect[NeighborhoodAnglesCompiled2D, NeighborhoodAnglesCompiled3D];

(* #NeighborhoodBisectorsCompiled *)
NeighborhoodBisectorsCompiled = Compile[{{x0, _Real, 1}, {xnei, _Real, 2}},
  If[Length[xnei] == 0,
    {{}},
    With[
      {x = {#[[1]] - x0[[1]], #[[2]] - x0[[2]], #[[3]] - x0[[3]]}& @ Transpose[xnei]},
      With[
        {xu = (x / {#,#,#})& @ Sqrt @ Total[x^2]},
        With[
          {means = 0.5*(xu + RotateLeft /@ xu)},
          With[
            {mnorms = Sqrt @ Total[means^2]},
            Transpose[means / If[Length[mnorms] == 2, {#,#}, {#,#,#}]&[mnorms]]]]]]],
   RuntimeOptions -> {"Speed", "EvaluateSymbolically" -> False},
   Parallelization -> True];
Protect[NeighborhoodBisectorsCompiled];

(* #NeighborhoodEdgeLengthsCompiled *)
NeighborhoodEdgeLengthsCompiled = Compile[{{x0, _Real, 1}, {x, _Real, 2}},
  Sqrt[Total[MapThread[Subtract, {x0, Transpose[x]}]^2]],
  RuntimeOptions -> {"Speed", "EvaluateSymbolically" -> False},
  Parallelization -> True];
Protect[NeighborhoodEdgeLengthsCompiled];

(* #NeighborhoodSort3DCompiled ********************************************************************)
NeighborhoodSort3DCompiled = Compile[
  {{x0, _Real, 1}, {neivec, _Complex, 1}, {nei, _Integer, 1}, {xnei, _Real, 2}},
  Part[
    nei,
    Ordering[
      Arg@Dot[
        neivec,
        MapThread[Subtract, {Transpose[xnei], x0}]]]],
  RuntimeOptions -> "Speed"];

(* #NeighborhoodSort3DCompiled ********************************************************************)
NeighborhoodSort2DCompiled = Compile[
  {{x0, _Real, 1}, {nei, _Integer, 1}, {xnei, _Real, 2}},
  Part[
    nei,
    Ordering[
      Arg@Dot[{1, I}, {#[[1]] - x0[[1]], #[[2]] - x0[[2]]}&@Transpose[xnei]]]],
  RuntimeOptions -> "Speed"];

(* #CalculateVertexNormals ************************************************************************)
CalculateVertexNormalsTr[Ft_List, Xt_List] := With[
  {nx3 = 3 * Length[First@Ft]},
  With[
    {u1 = Xt[[All, Ft[[2]]]] - Xt[[All, Ft[[1]]]],
     u2 = Xt[[All, Ft[[3]]]] - Xt[[All, Ft[[1]]]],
     P = SparseArray[
       Rule[
         Transpose@{Join[Ft[[1]], Ft[[2]], Ft[[3]]], Range[nx3]}, 
         ConstantArray[1, nx3]]]},
    With[
      {q = Transpose[
         {u1[[2]]*u2[[3]] - u1[[3]]*u2[[2]],
          u1[[3]]*u2[[1]] - u1[[1]]*u2[[3]],
          u1[[1]]*u2[[2]] - u1[[2]]*u2[[1]]}]},
      NormalizeRows @ Dot[P, Join[q, q, q]]]]];
CalculateVertexNormals[F_List, X_List] := Transpose @ CalculateVertexNormalsTr[
  Transpose @ F,
  Transpose @ X];
Protect[CalculateVertexNormalsTr, CalculateVertexNormals];

(* #CalculateFaceAngleCosinesTr *******************************************************************)
CalculateFaceAngleCosinesTr = Compile[
  {{Ft, _Integer, 2}, {Xtr, _Real, 2}},
  With[
    {sides = Transpose @ Map[
       Partition[#, Length@First@Ft]&,
       NormalizeColumns @ Subtract[
         Xtr[[All, Join[Ft[[2]], Ft[[3]], Ft[[1]]]]],
         Xtr[[All, Join[Ft[[1]], Ft[[2]], Ft[[3]]]]]]]},
    {Total[sides[[1]] * -sides[[3]]],
     Total[sides[[2]] * -sides[[1]]],
     Total[sides[[3]] * -sides[[2]]]}],
  RuntimeOptions -> {"Speed", "EvaluateSymbolically" -> False},
  Parallelization -> True];
(* #CalculateFaceAnglesTr *************************************************************************)
CalculateFaceAnglesTr[Ft_List, Xtr_List] := ArcCos @ CalculateFaceAngleCosinesTr[Ft, Xtr];
(* #CalculateFaceAngleCosines *********************************************************************)
CalculateFaceAngleCosines[F_List, X_List] := Transpose @ CalculateFaceAngleCosinesTr[
  Transpose @ F,
  Transpose @ X];
(* #CalculateFaceAngles ***************************************************************************)
CalculateFaceAngles[Ft_List, Xtr_List] := Transpose @ ArcCos @ CalculateFaceAngleCosinesTr[
  Transpose @ F,
  Transpose @ X];
Protect[CalculateFaceAngleCosinesTr, CalculateFaceAngleCosines, 
        CalculateFaceAnglesTr, CalculateFaceAngles];

(* #CalculateFaceNormals *)
CalculateFaceNormalsTr = Compile[
  {{Ft, _Integer, 2}, {Xt, _Real, 2}},
  With[
      {u1 = Xt[[All, Ft[[2]]]] - Xt[[All, Ft[[1]]]],
       u2 = Xt[[All, Ft[[3]]]] - Xt[[All, Ft[[1]]]]},
      With[
        {crosses = {
           u1[[2]]*u2[[3]] - u1[[3]]*u2[[2]],
           u1[[3]]*u2[[1]] - u1[[1]]*u2[[3]],
           u1[[1]]*u2[[2]] - u1[[2]]*u2[[1]]}},
        With[
          {norms = (Sqrt[#] + 1.0 - Sign[Chop[#]])& @ Total[crosses^2]},
          (* the number with the sign prevents divide by 0 for 0-norm cross-products *)
          crosses / {norms, norms, norms}]]],
  RuntimeOptions -> {"Speed", "EvaluateSymbolically" -> False},
  Parallelization -> True];
CalculateFaceNormals[F_List, X_List] := Transpose @ CalculateFaceNormalsTr[
  Transpose @ F,
  Transpose @ X];
Protect[CalculateFaceNormalsTr, CalculateFaceNormals];

(* #CalculateFaceAxes *)
CalculateFaceAxes3DTr = Compile[
  {{Ft, _Integer, 2}, {Xt, _Real, 2}, {normals, _Real, 2}},
  With[
    {side1s = Sqrt @ Total[(Xt[[All, Ft[[2]]]] - Xt[[All, Ft[[1]]]])^2]},
    Transpose @ {
      side1s,
      {normals[[2]] * side1s[[3]] - side1s[[2]] * normals[[3]],
       side1s[[1]] * normals[[3]] - normals[[1]] * side1s[[3]],
       side1s[[2]] * normals[[1]] - normals[[2]] * side1s[[1]]}}],
  RuntimeOptions -> {"Speed", "EvaluateSymbolically" -> False},
  Parallelization -> True];
CalculateFaceAxes2DTr = Compile[
  {{Ft, _Integer, 2}, {Xt, _Real, 2}},
  With[
    {side1s = Sqrt @ Total[(Xt[[All, Ft[[2]]]] - Xt[[All, Ft[[1]]]])^2]},
    Transpose @ {side1s, {-side1s[[2]], side1s[[1]]}}],
  RuntimeOptions -> {"Speed", "EvaluateSymbolically" -> False},
  Parallelization -> True];
CalculateFaceAxesTr[Ft_List, Xt_List, normt_List] := CalculateFaceAxes3DTr[Ft, Xt, normt];
CalculateFaceAxes[F_List, X_List, norm_List] := CalculateFaceAxes3DTr[
  Transpose @ F,
  Transpose @ X,
  Transpose @ norm];
CalculateFaceAxesTr[Ft_List, Xt_List] := If[Length[Xt] == 2,
  CalculateFaceAxes2DTr[Ft, Xt],
  CalculateFaceAxes3DTr[Ft, Xt, CalculateFaceNormalsTr[Ft, Xt]]];
CalculateFaceAxes[F_List, X_List] := If[Length[First@X] == 2,
  CalculateFaceAxes2DTr[Transpose @ F, Transpose @ X],
  With[
    {Ft = Transpose@F, Xt = Transpose@X}, 
    CalculateFaceAxes3DTr[Ft, Xt, CalculateFaceNormalsTr[Ft, Xt]]]];
Protect[CalculateFaceAxes, CalculateFaceAxesTr, CalculateFaceAxes3DTr, CalculateFaceAxes2DTr];

(* #CalculateFacePlaneCoordinates *)
CalculateFacePlaneCoordinatesTr = Compile[
  {{Ft, _Integer, 2}, {Xt, _Real, 2}, {At, _Real, 3}},
  Outer[
    Total[#1 * #2]&,
    Transpose[At],
    {Xt[[All, Ft[[1]]]], Xt[[All, Ft[[2]]]], Xt[[All, Ft[[3]]]]},
    1, 1],
  RuntimeOptions -> {"Speed", "EvaluateSymbolically" -> False},
  Parallelization -> True];
CalculateFacePlaneCoordinates[F_List, X_List, A_List] := CalculateFacePlaneCoordinatesTr[
  Transpose @ F,
  Transpose @ X,
  Transpose[A, {3,2,1}]];
Protect[CalculateFacePlaneCoordinates, CalculateFacePlaneCoordinatesTr];

(* #CalculateFaceBisectorsTr *)
CalculateFaceBisectorsTr = Compile[
  {{Ft, _Integer, 2}, {Xt, _Real, 2}},
  With[
    {sides = {
       Xt[[Ft[[2]]]] - Xt[[Ft[[1]]]],
       Xt[[Ft[[3]]]] - Xt[[Ft[[2]]]],
       Xt[[Ft[[1]]]] - Xt[[Ft[[3]]]]}},
    With[
      {norms = Sqrt @ Map[Total[#^2]&, sides]},
      With[
        {U = MapThread[(#1 / Table[#2, {Length@#1}])&, {sides, norms}]},
        With[
          {sums = {U[[1]] - U[[3]], U[[2]] - U[[1]], U[[3]] - U[[2]]}},
          sums / Table[#, {Length@Xt}]& /@ Sqrt[Total /@ sums^2]]]]],
  RuntimeOptions -> {"Speed", "EvaluateSymbolically" -> False},
  Parallelization -> True];
CalculateFaceBisectors[F_List, X_List] := Transpose @ CalculateFaceBisectorsTr[
  Transpose @ F,
  Transpose @ X];
Protect[CalculateFaceBisectorsTr, CalculateFaceBisectors];

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

(* #CorticalMeshOrientForMap
 * Yields the coordinates of the subset of vertices subject to the map (via exclusions) rotated such
 * that the center is at {1,0,0} and the orient point (if there is one) is in the (+x, y)
 * half-plane.
 *)
CorticalMeshOrientForMap[Xt_, center_] := With[
  {RMain = If[center[[1]] === None || center[[1]] === Automatic, 
    IdentityMatrix[3],
    RotationMatrix[{center[[1]], {1,0,0}}]]},
  Dot[
    If[center[[2]] === None || center[[2]] === Automatic,
      IdentityMatrix[3], 
      With[
        {orientPt = Dot[RMain, center[[2]]]},
        RotationMatrix[-ArcTan[orientPt[[1]], orientPt[[2]]], {1,0,0}]]],
    RMain,
    Xt]];

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
             {idx = VertexIndex[mesh, #]},
             If[idx === $Failed, $Failed, VertexCoordinates[mesh][[idx]]]],
           {_?NumericQ, _?NumericQ, _?NumericQ}, RegionNearest[MeshRegion[mesh], #],
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
  {X0t = CorticalMeshOrientForMap[VertexCoordinatesTr[mesh], center],
   EPt = EdgePairsTr[mesh],
   Ft = FacePairsTr[mesh]},
  Switch[
    If[StringQ[method], ToLowerCase[method], method],
    "mollenweide"|"equirectangular"|"mercator", Pick[
      Transpose @ EPt,
      Plus[
        Sign[X0t[[1, EPt[[1]]]]],
        Sign[X0t[[1, EPt[[2]]]]],
        Sign[X0t[[2, EPt[[1]]]] * X0t[[2, EPt[[2]]]]]],
      -3],
    "orthographic", Pick[VertexList[mesh], Sign[X0t[[1]]], -1],
    "polarstretching", With[
      {id = First @ Nearest[VertexCoordinates[mesh] -> Automatic, -center[[1]]]},
      With[
        {faces = VertexFaceList[mesh][[id]],
         FX = FaceCoordinates[mesh]},
        {Part[
           FaceList[mesh],
           First @ SortBy[
             faces,
             RegionDistance[Polygon[FX[[#]]], -center[[1]]]&]]}]],
    _, {}]];
Protect[CorticalMapAutomaticExclusions];

(* #CorticalMapTranslateExclusions
 * Yields the lists of {vertices, edges, faces} that are to be included in the projection; these
 * are only guaranteed that if you, for example, ask for VertexCoordinates[mesh][[vertices]] the
 * result will be correct; in this vein, vertices may be All.
 *)
CorticalMapTranslateExclusions[mesh_, method_, center_, excl_, prad_] := Check[
  Catch @ With[
    {excl0 = Apply[
       Function[{vertices, edges, faces},
         {VertexIndex[mesh, vertices],
          EdgeIndex[mesh, ReplaceAll[edges, UndirectedEdge -> List]],
          FaceIndex[mesh, faces]}],
       Union[Rest[#]]& /@ GatherBy[
         Join[
           (* these enforce the ordering [vertices, edges, faces] in the result *)
           {0, {0,0}, {0,0,0}},
           Which[
             ListQ[excl], Replace[
               excl, 
               Automatic :> Sequence@@CorticalMapAutomaticExclusions[
                 mesh, method, center, excl, prad],
               {1}],
             excl === Automatic, CorticalMapAutomaticExclusions[
               mesh, method, center, excl, prad],
             True, Throw[{$Failed,$Failed,$Failed}]]],
         Function@Which[
           IntegerQ[#], 1,
           Head[#] === UndirectedEdge, 2,
           ListQ[#] && Length[#] == 2, 2,
           ListQ[#] && Length[#] == 3, 3,
           True, Throw[{$Failed, $Failed, $Failed}]]]]},
    With[
      {incl0 = MapThread[
         Complement,
         {{Range[VertexCount[mesh]], Range[EdgeCount[mesh]], Range[FaceCount[mesh]]},
          excl0}]},
      With[
        {incl = If[prad === Full || prad === All,
           incl0,
           ReplacePart[
             incl0,
             1 -> If[QuantityQ[prad] || ListQ[prad],
               (* we have an angular radius requirement; just use distances *)
               Pick[
                 incl0[[1]],
                 Sign @ Subtract[
                   Which[
                     prad[[2]] == "AngularDegrees", Cos[prad[[1]] * Pi / 180.0],
                     prad[[2]] == "Radians", Cos[prad[[1]]],
                     True, Message[
                       CorticalMap::badarg,
                       "Radius quantities must be in degrees or radians"]],
                   Total @ MapThread[
                     Times,
                     {NormalizeColumns @ Part[VertexCoordinatesTr[mesh], All, incl0[[1]]],
                      Normalize @ center[[1]]}]],
                  -1],
               (* we have a radius requirement; find everything with a small shortest path dist *)
               True, Intersection[
                 incl0[[1]],
                 CorticalMeshNeighborhood[
                   mesh,
                   First @ Ordering[
                     Total[
                       MapThread[
                         Subtract,
                         {Transpose@VertexCoordinates[mesh], center[[1]]}]^2],
                     1],
                   prad]]]]]},
        With[
          {Vs = incl[[1]],
           Es = Intersection[
             incl[[2]],
             Pick @@ Append[
               Transpose@Tally[Join @@ Part[VertexEdgeList[mesh], incl[[1]]]],
               2]]},
          With[
            {Fs = Intersection[
               incl[[3]],
               Pick @@ Append[Transpose@Tally[Join @@ Part[VertexFaceList[mesh], Vs]], 3],
               Pick @@ Append[Transpose@Tally[Join @@ Part[EdgeFaceList[mesh], Es]], 3]]},
            {Vs, Es, Fs}]]]]],
  {$Failed, $Failed, $Failed}];
Protect[CorticalMapTranslateExclusions];

(* #CorticalMapTranslateMethod
 * Yields a translation function given the Method argument and the translated parameters: Center,
 * vertex-exclusions, projection area, and projection radius.
 *)
CorticalMapTranslateMethod[mesh_, method_, center_, incl_, prad_] := Check[
  With[
    {projFn = Replace[
       ToLowerCase[method],
       (* Here we actually define the functions for transformation;
          these should accept a list of vertices (already centered such that the center lies at 
          (1,0,0) and the orient point lies in the <positive X, Y> half-plane) and should return
          the 2D coordinates (while ignoring the orient point). *)
       {"mollenweide" :> Function[{X},
          With[
            {S = Transpose @ ConvertCoordinates[Transpose@X, Cartesian -> {Longitude, Latitude}],
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
              meshRadius * Sqrt[2.0] * {2.0 / Pi * S[[1]] * Cos[th], Sin[th]}]]],
        "equirectangular" :> Function[{X},
          With[
            {S = ConvertCoordinates[Transpose@X, Cartesian -> {Longitude, Latitude}],
             meshRadius = Mean[Sqrt[Total[Transpose[X]^2]]]},
            Transpose[meshRadius * S]]],
        "mercator" :> Function[{X},
          With[
            {S = Transpose @ ConvertCoordinates[Transpose@X, Cartesian -> {Longitude, Latitude}],
             meshRadius = Mean[Sqrt[Total[Transpose[X]^2]]]},
            With[
              {sinPhi = Sin[S[[2]]]},
              meshRadius * {S[[1]], 0.5*Log[(1 + sinPhi) / (1 - sinPhi)]}]]],
        "orthographic" :> Function[{X}, X[[2;;3]]],
        "polarstretching" :> Function[{X},
          With[
            {polar = Transpose @ ConvertCoordinates[
               Transpose @ X[[{2,3,1}]],
               Cartesian -> {Longitude, SphericalPolarAngle}]},
            Global`dbg = {X, polar};
            {# * Cos[polar[[1]]], # * Sin[polar[[1]]]}& @ Tan[polar[[2]] / 2]]],
        "graph" :> With[
          {em = Transpose @ GraphEmbedding[
             Graph[
               Part[VertexList[mesh], incl[[1]]],
               Part[EdgeList[mesh], incl[[2]]],
               EdgeWeight -> Part[EdgeLengths[mesh], incl[[2]]],
               GraphLayout -> {"SpringElectricalEmbedding", "EdgeWeighted" -> True}]],
           idx = VertexList[mesh][[incl[[1]]]],
           coords = VertexCoordinatesTr[mesh][[All, incl[[1]]]],
           near = Nearest[VertexCoordinates[mesh][[incl[[1]]]] -> Automatic]},
          Function[{Xt},
            If[Dimensions[Xt] === {3, Length[incl[[1]]]},
              em,
              Message[CorticalMap::badarg, "Graph embeddings do not support transforms"]]]],
        _ :> Message[CorticalMap::badarg, "Could not recognize projection type"]}],
     RMtx = RotationMatrix[{center[[1]], {1,0,0}}]},
    With[
      {orientFn = If[center[[2]] === Automatic,
         Function[RotationMatrix[{First@Eigenvectors[Covariance[Transpose@#],1], {1,0}}] . #],
         With[
           {orientRMtx = RotationMatrix[
              {projFn[CorticalMeshOrientForMap[List /@ center[[2]], center]][[All, 1]],
               {1,0}}]},
           Function[orientRMtx . #]]]},
      Function[
        If[Length[#] == 3,
          orientFn @ projFn[CorticalMeshOrientForMap[#, center]],
          Transpose @ orientFn[projFn[CorticalMeshOrientForMap[Transpose @ #, center]]]]]]],
  $Failed];
Protect[CorticalMapTranslateMethod];



(**************************************************************************************************)
(**************************************************************************************************)
(**************** CorticalMesh, CorticalMap, and the CortexPlots are below ************************)
(**************************************************************************************************)
(**************************************************************************************************)



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
  CorticalMesh[V_List, F_List, OptionsPattern[]] :> mesh,
  (* We do want to declare a few locals here *)
  With[
    {optionalProperties = Last@Reap[
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
     cleaned = With[
       {step1 = ReplaceRepeated[{V, F}, Property[x_, _] :> x]},
       {step1[[1]], 
        With[
          {Ftr = Transpose @ step1[[2]]},
          Transpose @ Union[
            Sort /@ Transpose @ {
              Join[Ftr[[1]],Ftr[[1]],Ftr[[2]]],
              Join[Ftr[[2]],Ftr[[3]],Ftr[[3]]]}]],
        step1[[2]]}]},
    {(* ======================================= Privates ======================================== *)

     (* #VertexProperties, #EdgeProperties, #FaceProperties *)
     VertexProperties[mesh] = ParseMeshProperties[V, optionalProperties[[1]]],
     EdgeProperties[mesh] = ParseMeshProperties[Transpose @ cleaned[[2]], optionalProperties[[2]]],
     FaceProperties[mesh] = ParseMeshProperties[F, optionalProperties[[3]]],
     
     (* #FaceIndexArray [private] *)
     FaceIndexArray[mesh] :> With[
       {Ft = FaceListTr[mesh],
        RF = Range[FaceCount[mesh]]},
       SparseArray @ Rule[
         Transpose @ MapThread[
           Join,
           {Ft,
            {Ft[[2]], Ft[[3]], Ft[[1]]},
            {Ft[[3]], Ft[[1]], Ft[[2]]},
            {Ft[[3]], Ft[[2]], Ft[[1]]},
            {Ft[[2]], Ft[[1]], Ft[[3]]},
            {Ft[[1]], Ft[[3]], Ft[[2]]}}],
         Join[RF, RF, RF, RF, RF, RF]]],
     (* #EdgeIndexArray [private] *)
     EdgeIndexArray[mesh] :> With[
       {Et = EdgePairsTr[mesh],
        RE = Range[EdgeCount[mesh]]},
       SparseArray[Transpose[MapThread[Join, {Et, Reverse[Et]}]] -> Join[RE, RE]]],
     (* VertexIndexArray [private] *)
     VertexIndexArray[mesh] :> If[VertexList[mesh] == Range[VertexCount[mesh]],
       VertexList[mesh],
       Normal @ SparseArray[VertexList[mesh] -> Range[VertexCount[mesh]]]],
     
     (* #SumOverFacesMatrix and #SumOverEdgesMatrix *)
     SumOverFacesMatrix[mesh] :> SparseArray[
       Transpose[{Range[3*FaceCount[mesh]], Join @@ FaceListTr[mesh]}] -> 1],
     SumOverEdgesMatrix[mesh] :> SparseArray[
       Transpose[{Range[2*EdgeCount[mesh]], Join @@ EdgePairsTr[mesh]}] -> 1],


     (* ======================================= Settables ======================================= *)

     (* All data is actually stored transposed for efficiency *)
     VertexCoordinatesTr[mesh] = Transpose @ cleaned[[1]],
     FaceListTr[mesh]          = Transpose @ cleaned[[3]],
     Options[mesh]             = Cases[
       Options[CorticalMesh][[All,1]],
       opt:Except[Properties] :> (opt -> OptionValue[opt]),
       {1}],


     (* ====================================== Immediates ======================================= *)

     (* #EdgePairs is deducible strictly from the FaceList *)
     EdgePairsTr[mesh] -> cleaned[[2]],

     (* #VertexCount *)
     VertexCount[mesh] -> Length @ First @ VertexCoordinatesTr[mesh],
     (* #EdgeCount *)
     EdgeCount[mesh] -> Length @ First @ EdgePairsTr[mesh],
     (* #FaceCount *)
     FaceCount[mesh] -> Length @ First @ FaceListTr[mesh],

     (* VertexList is not settable but is public thus depends on the current vertex coordinates *)
     VertexList[mesh] -> Range[VertexCount[mesh]],

     (* #CorticalMeshQ tests whether the mesh is valid and yields True if so *)
     CorticalMeshQ[mesh] -> Which[
       !ArrayQ[VertexCoordinatesTr[mesh], 2, NumericQ], Message[
         CorticalMesh::badarg,
         "VertexCoordinates must be a 2D numeric array"],
       !MatchQ[Dimensions[VertexCoordinatesTr[mesh]], {3, _}], Message[
         CorticalMesh::badarg,
         "VertexCoordinates must contain 3D points"],
       !ArrayQ[FaceListTr[mesh], 2, IntegerQ], Message[
         CorticalMesh::badarg, 
         "FaceList conatins non-integer values"],
       Length[FaceListTr[mesh]] != 3, Message[
         CorticalMesh::badarg,
         "FaceList must contain triangles only"],
       Min[Join@@FaceListTr[mesh]] < 1, Message[
         CorticalMesh::badarg, 
         "FaceList conatins values less than 1"],
       Max[Join@@FaceListTr[mesh]] > VertexCount[mesh], Message[
         CorticalMesh::badarg, 
         "FaceList conatins values greater than the number of vertices"],
       And[
         Options[mesh] =!= Automatic,
         Complement[
           Options[mesh][[All,1]], 
           Cases[Options[CorticalMesh][[All, 1]], Except[Properties], {1}]] != {}], Message[
           CorticalMesh::badarg,
           "Unrecognized option given to CorticalMesh"],
       True, True],



     (* ======================================== Delays ========================================= *)

     (* These indices depend only on the face list and edge pairs and tell us to which faces a
      * vertex or an edge belongs.
      *)
     VertexEdgeList[mesh] :> With[
       {RE = Range[EdgeCount[mesh]],
        EP = VertexIndex[mesh, EdgePairsTr[mesh]]},
       Part[
         SplitBy[
           SortBy[
             Transpose[{Join[EP[[1]], EP[[2]]], Join[RE, RE]}],
             First],
           First],
         All, All, 2]],
     VertexFaceList[mesh] :> With[
       {T = VertexIndex[mesh, FaceListTr[mesh]],
        RT = Range[FaceCount[mesh]]},
       Part[
         SplitBy[
           SortBy[
             Transpose[{Join @@ T, Join[RT, RT, RT]}],
             First],
           First],
         All, All, 2]],
     EdgeFaceList[mesh] :> With[
       {Ft = VertexIndex[mesh, Transpose @ FaceList[mesh]],
        RT = Range[FaceCount[mesh]],
        eidx = EdgeIndexArray[mesh]},
       With[
         {edges = Extract[
            eidx,
            Transpose @ {
              Join[Ft[[1]], Ft[[2]], Ft[[3]]],
              Join[Ft[[2]], Ft[[3]], Ft[[1]]]}]},
         Part[
           SplitBy[
             SortBy[
               Transpose[{edges, Join[RT, RT, RT]}],
               First],
             First],
           All, All, 2]]],

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
             {{FirstCase[VertexProperties[mesh], (Rule|RuleDelayed)[name,val_], $Failed],
               FirstCase[EdgeProperties[mesh], (Rule|RuleDelayed)[name,val_],   $Failed],
               FirstCase[FaceProperties[mesh], (Rule|RuleDelayed)[name,val_],   $Failed]},
              {VertexList, EdgeList, FaceList}}]],
         allPropNames]],



     (* ======================================== Methods ======================================== *)

     (* A few simple extensions of the immediates that count things *)
     VertexCount[mesh, patt_] := Count[VertexList[mesh], patt, {1}],
     FaceCount[mesh, patt_] := Count[FaceList[mesh], patt, {1}],
     EdgeCount[mesh, patt_] := Count[EdgePairs[mesh], patt, {1}],
     
     (* EdgeList is here just for the compatibility with graph-like operations (use EdgePairs) *)
     EdgeList[mesh] := ReplacePart[EdgePairs[mesh], {_,0} -> UndirectedEdge],
     EdgeList[mesh, patt_] := Cases[EdgeList[mesh], patt, {1}],
     
     (* Untransposed accessors to basic data *)
     VertexCoordinates[mesh] := Transpose @ VertexCoordinatesTr[mesh],
     FaceList[mesh] := Transpose @ FaceListTr[mesh],
     EdgePairs[mesh] := Transpose @ EdgePairsTr[mesh],
     EdgePairs[mesh, patt_] := Cases[Transpose @ EdgePairsTr[mesh], patt, {1}],
     FaceList[mesh, patt_] := Cases[Transpose @ FaceListTr[mesh], patt, {1}],
     VertexList[mesh, patt_] := Cases[VertexList[mesh], patt, {1}],
     
     (* #PropertyList method extension *)
     PropertyList[mesh] := Join[Properties[mesh][[All, 1]], {EdgeWeight, VertexCoordinates}],
     
     (* A few simple extensions of the private indices above *)
     (* #FaceIndex *)
     FaceIndex[mesh, {a_Integer, b_Integer, c_Integer}] := With[
       {id = FaceIndexArray[mesh][[a,b,c]]},
       If[id == 0, $Failed, id]],
     FaceIndex[mesh, list_List] := With[
       {id = FaceIndexArray[mesh]},
       ReplaceAll[
         Map[Extract[id, #]&, list, {-2}],
         0 -> $Failed]],
     (* EdgeIndex *)
     EdgeIndex[mesh, (List|UndirectedEdge)[a_Integer, b_Integer]] := With[
       {id = EdgeIndexArray[mesh][[a,b]]},
       If[id == 0, $Failed, id]],     
     EdgeIndex[mesh, list_List] := With[
       {id = EdgeIndexArray[mesh]},
       ReplaceAll[
         Map[Extract[id, #]&, list, {-2}],
         0 -> $Failed]],
     (* VertexIndex *)
     VertexIndex[mesh, i_Integer] := With[
       {id = If[i > Length[VertexIndexArray[mesh]], 0, VertexIndexArray[mesh][[i]]]},
       If[id == 0, $Failed, id]],
     VertexIndex[mesh, is_List] := With[
       {idx = VertexIndexArray[mesh]},
       Map[Part[idx, #]&, is, {-2}]],

     (* extensions of the opposite indices; edge/face lists *)
     VertexEdgeList[mesh, i_Integer] := Part[VertexEdgeList[mesh], VertexIndex[mesh, i]],
     VertexEdgeList[mesh, l_List] := With[
       {idx = VertexEdgeList[mesh]},
       Map[Part[idx, #]&, VertexIndex[mesh, i], {-2}]],
     VertexFaceList[mesh, i_Integer] := Part[VertexFaceList[mesh], VertexIndex[mesh, i]],
     VertexFaceList[mesh, l_List] := With[
       {idx = VertexFaceList[mesh]},
       Map[Part[idx, #]&, VertexIndex[mesh, i], {-2}]],
     EdgeFaceList[mesh, e_] := With[
       {idx = EdgeFaceList[mesh]},
       Which[
         Head[e] === UndirectedEdge, Part[idx, EdgeIndex[mesh, e]],
         MatchQ[e, {_Integer, _Integer}], Part[idx, EdgeIndex[mesh, e]],
         ListQ[e], Map[Part[idx, #]&, e, {-2}],
         True, $Failed]],
     
     (* #VertexDegree *)
     VertexDegree[mesh] := Length /@ NeighborhoodList[mesh],
     VertexInDegree[mesh] := VertexDegree[mesh],
     VertexOutDegree[mesh] := VertexDegree[mesh],
     VertexDegree[mesh, i_Integer] := VertexDegree[mesh][[VertexIndex[mesh, i]]],
     VertexInDegree[mesh, i_Integer] := VertexDegree[mesh][[VertexIndex[mesh, i]]],
     VertexOutDegree[mesh, i_Integer] := VertexDegree[mesh][[VertexIndex[mesh, i]]],

     (* #SumOverFaces *)
     SumOverFacesTr[mesh, datat_] := Dot[
       MapThread[Join, {datat, datat, datat}],
       SumOverFacesMatrix[mesh]],
     SumOverFaces[mesh, data_] := Dot[Transpose @ Join[data, data, data], SumOverFacesMatrix[mesh]],
     SumOverFaceVerticesTr[mesh, datat_] := Dot[Join @@ datat, SumOverFacesMatrix[mesh]],
     SumOverFaceVertices[mesh, data_] := Dot[Join @@ Transpose[data], SumOverFacesMatrix[mesh]],
     (* #SumOverEdges *)
     SumOverEdgesTr[mesh, datat_] := Dot[
       MapThread[Join, {datat, datat}],
       SumOverEdgesMatrix[mesh]],
     SumOverEdges[mesh, data_] := Dot[
       Transpose @ Join[data, data],
       SumOverEdgesMatrix[mesh]],
     SumOverEdgesDirectedTr[mesh, datat_] := Dot[
       MapThread[Join, {datat, -datat}],
       SumOverEdgesMatrix[mesh]],
     SumOverEdgesDirected[mesh, data_] := Dot[
       Transpose @ Join[data, -data],
       SumOverEdgesMatrix[mesh]],
     SumOverEdgeVerticesTr[mesh, datat_] := Dot[Join @@ datat, SumOverEdgesMatrix[mesh]],
     SumOverEdgeVertices[mesh, data_] := Dot[Join @@ Transpose[data], SumOverEdgesMatrix[mesh]],


     (* ================================ Geometrical Properties ================================= *)

     (* ------------------------------------------ Faces ---------------------------------------- *)

     (* #FaceAngles *)
     FaceAngleCosinesTr[mesh] :> CalculateFaceAngleCosinesTr[
       FaceListTr[mesh],
       VertexCoordinatesTr[mesh]],
     FaceAngleCosinesTr[mesh, Xtr_] := CalculateFaceAngleCosinesTr[FaceListTr[mesh], Xtr],
     FaceAnglesTr[mesh, Xtr_]       := CalculateFaceAnglesTr[FaceListTr[mesh], Xtr],
     FaceAngleCosines[mesh, Xx_]    := CalculateFaceAngleCosines[FaceList[mesh], Xx],
     FaceAngles[mesh, Xx_]          := CalculateFaceAngles[FaceList[mesh], Xx],
     FaceAnglesTr[mesh]       := ArcCos @ FaceAngleCosinesTr[mesh],
     FaceAngleCosines[mesh]   := Transpose @ FaceAngleCosinesTr[mesh],
     FaceAngles[mesh]         := Transpose @ ArcCos @ FaceAngleCosinesTr[mesh],

     (* #FaceNormals *)
     FaceNormalsTr[mesh] :> CalculateFaceNormalsTr[FaceListTr[mesh], VertexCoordinatesTr[mesh]],
     FaceNormals[mesh] := Transpose @ FaceNormalsTr[mesh],
     FaceNormalsTr[mesh, Xtr_] := CalculateFaceNormalsTr[FaceListTr[mesh], Xtr],
     FaceNormals[mesh, X_] := Transpose @ CalculateFaceNormalsTr[FaceListTr[mesh], Transpose @ X],
     
     (* #FaceAxes *)
     FaceAxesTr[mesh] :> CalculateFaceAxes3DTr[
       FaceListTr[mesh],
       VertexCoordinatesTr[mesh],
       FaceNormalsTr[mesh]],
     FaceAxes[mesh] := Transpose @ FaceAxesTr[mesh],
     FaceAxesTr[mesh, Xt_] := CalculateFaceAxesTr[FaceListTr[mesh], Xt],
     FaceAxes[mesh, Xx_] := Transpose @ CalculateFaceAxesTr[FaceListTr[mesh], Transpose @ Xx],

     (* #FaceBisectors *)
     FaceBisectorsTr[mesh] :> CalculateFaceBisectorsTr[FaceListTr[mesh], VertexCoordinatesTr[mesh]],
     FaceBisectors[mesh] := Transpose @ FaceBisectorsTr[mesh],
     FaceBisectorsTr[mesh, X_] := CalculateFaceBisectorsTr[FaceListTr[mesh], X],
     FaceBisectors[mesh, X_] := Transpose @ FaceBisectorsTr[mesh, X],

     (* #FaceCoordinates *)
     FaceCoordinatesTr[mesh] :> With[
       {Xt = VertexCoordinatesTr[mesh],
        Ft = VertexIndex[mesh, FaceListTr[mesh]]},
       Transpose @ {Xt[[All, Ft[[1]]]], Xt[[All, Ft[[2]]]], Xt[[All, Ft[[3]]]]}],
     FaceCoordinates[mesh] := Transpose[FaceCoordinatesTr[mesh], {3,2,1}],

     (* #FacePlaneCoordinates *)
     FacePlaneCoordinatesTr[mesh] :> CalculateFacePlaneCoordinatesTr[
       FaceListTr[mesh],
       VerteListTr[mesh],
       FaceAxesTr[mesh]],
     FacePlaneCoordinatesTr[mesh, Xt_] := CalculateFacePlaneCoordinatesTr[
       FaceListTr[mesh],
       Xt,
       FaceAxesTr[mesh, Xt]],
     FacePlaneCoordinates[mesh] := Transpose @ FacePlaneCoordinatesTr[mesh],
     FacePlaneCoordinates[mesh, Xx_] := Transpose @ FacePlaneCoordinatesTr[mesh, Transpose @ Xx],

     (* ------------------------------------------ Edges ---------------------------------------- *)

     (* #EdgeCoordinates *)
     EdgeCoordinatesTr[mesh] :> With[
       {Xt = VertexCoordinatesTr[mesh],
        EL = VertexIndex[mesh, EdgePairsTr[mesh]]},
       Transpose @ {Xt[[All, EL[[1]]]], Xt[[All, EL[[2]]]]}],
     EdgeCoordinates[mesh] := Transpose[EdgeCoordinatesTr[mesh], {3,2,1}],
     
     (* #EdgeLengths *)
     EdgeLengths[mesh] :> With[
       {EL = EdgePairsTr[mesh],
        Xt = VertexCoordinatesTr[mesh]},
       Sqrt @ Total[(Xt[[All, EL[[1]]]] - Xt[[All, EL[[2]]]])^2]],
     EdgeLengths[mesh, X_] := With[
       {EL = EdgePairsTr[mesh],
        Xt = Transpose @ X},
       Sqrt @ Total[(Xt[[All, EL[[1]]]] - Xt[[All, EL[[2]]]])^2]],
     EdgeLengthsTr[mesh, Xt_] := With[
       {EL = EdgePairsTr[mesh]},
       Sqrt @ Total[(Xt[[All, EL[[1]]]] - Xt[[All, EL[[2]]]])^2]],
     EdgeWeight[mesh] := EdgeLengths[mesh],
     EdgeWeight[mesh, e:(List|UndirectedEdge)[_Integer, _Integer]] := Part[
       EdgeLengths[mesh],
       EdgeIndex[mesh, e]],
     EdgeWeight[mesh, es:{{_Integer,_Integer}..}] := Part[
       EdgeLengths[mesh],
       EdgeIndex[mesh, es]],
     EdgeWeight[mesh, es:{(List|UndirectedEdge)[_Integer, _Integer]..}] := Part[
       EdgeLengths[mesh],
       EdgeIndex[mesh, ReplaceAll[es, UndirectedEdges -> List]]],

     (* ---------------------------------------- Vertices --------------------------------------- *)

     (* #VertexNormals *)
     VertexNormalsTr[mesh] :> NormalizeColumns @ SumOverFacesTr[mesh, FaceNormalsTr[mesh]],
     VertexNormalsTr[mesh, Xt_] := NormalizeColumns @ SumOverFacesTr[mesh, FaceNormalsTr[mesh, Xt]],
     VertexNormals[mesh] := Transpose @ VertexNormalsTr[mesh],
     VertexNormals[mesh, Xx_] := Transpose @ VertexNormalsTr[mesh, Transpose @ Xx],

     (* ------------------------------------- Neighborhoods ------------------------------------- *)

     (* #NeighborhoodList *)
     NeighborhoodList[mesh] :> With[
       {X = VertexCoordinates[mesh],
        E = EdgePairsTr[mesh]},
       With[
         {neivecs = NeighborhoodVectors[mesh],
          neis = Part[
            SplitBy[
              SortBy[
                Transpose[{Join[E[[1]], E[[2]]], Join[E[[2]], E[[1]]]}],
                First],
              First],
            All, All, 2]},
         MapThread[NeighborhoodSort3DCompiled, {X, neivecs, neis, X[[#]] & /@ neis}]]],
     (* #NeighborhoodVectors [private] *)
     NeighborhoodVectors[mesh] := With[
       {Xt = VertexNormalsTr[mesh]},
       With[
         {x = Xt[[1]], y = Xt[[2]], z = Xt[[3]],
          norms = Sqrt@Total[Xt^2]},
         Transpose[
           {(y + I*x*z/norms)/(I*x + y),
            (I*x + y*z/norms)/(x - I*y),
            -(x + I*y)/norms}]]],
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


     (* ======================================= Interfaces ====================================== *)
   
     (* #MeshRegion *)
     MeshRegion[mesh] :> MeshRegion[VertexCoordinates[mesh], Polygon[FaceList[mesh]]],
     
     (* #BoundaryMeshRegion *)
     BoundaryMeshRegion[mesh] :> BoundaryMeshRegion[
       VertexCoordinates[mesh],
       Polygon[FaceList[mesh]]],
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
     AdjacencyMatrix[mesh, opts___]             := AdjacencyMatrix[Graph[mesh], opts],
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
     VertexOutDegree[mesh, opts___]             := VertexDegree[mesh]}],

  SetSafe -> True,
  Symbol -> CorticalMesh3D];

(* We want to create our own version of the Clone constructor for modifying things like options.
 * Note that there are only 6 settables in the CorticalMesh immutable structure:
 *  - VertexProperties
 *  - EdgeProperties
 *  - FaceProperties
 *  - VertexCoordinatesTr
 *  - FaceListTr
 *  - Options
 *)
$CorticalMeshReconstructionOptions = Join[
  {VertexCoordinates, FaceList},
  Cases[Options[CorticalMesh][[All,1]], Except[Properties], {1}]];
Protect[$CorticalMeshReconstructionOptions];
CorticalMesh[mesh_?CorticalMeshQ, args___Rule] := Check[
  With[
    {optsarg = (
       If[Complement[#[[All,1]], $CorticalMeshReconstructionOptions] == {},
         #,
         Message[
           CorticalMesh::badarg,
           StringJoin[
             "CorticalMesh[] reconstructor may only be passed CorticalMesh options or",
             " VertexCoordinates or FaceList"]]]
       )& @ {args}},
    If[Length[optsarg] == 0,
      mesh,
      With[
        {vtx = Replace[VertexCoordinates, optsarg],
         fl = Replace[FaceList, optsarg],
         opts = Fold[
           Function @ With[
             {edit = Replace[#1, (Rule|RuleDelayed)[#2[[1]], _] :> #2, {1}]},
             If[SameQ[edit, #1], Append[#1, #2], edit]],
           Options[mesh],
           Select[optsarg, (#[[1]] =!= VertexCoordinates && #[[1]] =!= FaceList)&]]},
        Clone[
          mesh,
          Sequence @@ Join[
            If[vtx =!= VertexCoordinates, 
              {VertexCoordinatesTr -> If[Length[vtx] == 3, vtx, Transpose[vtx]]},
              {}],
            If[fl =!= FaceList,
              {FaceListTr -> If[Length[fl] == 3, fl, Transpose[fl]]},
              {}],
            If[Length[opts] > 0,
              {Options -> opts},
              {}]]]]]],
  $Failed];

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
   Method -> "Orthographic",
   Center -> Automatic,
   Exclusions -> Automatic,
   Radius -> Full}];
DefineImmutable[
  CorticalMap[mesh_?CorticalMeshQ, OptionsPattern[]] :> map,
  With[
    {optionalProperties = Last@Reap[
       Replace[
         OptionValue[Properties] /. None -> {},
         {r:Rule[_Integer, {_Rule..}] :> Sow[r, "Vertex"],
          r:Rule[{_Integer, _Integer}, {_Rule..}] :> Sow[r, "Edge"],
          Rule[UndirectedEdge[a_Integer, b_Integer], rs:{_Rule..}] :> Sow[{a,b} -> rs, "Edge"],
          r:Rule[{_Integer,_Integer,_Integer}, {_Rule..}] :> Sow[r, "Face"],
          x_ :> Message[CorticalMap::badarg, "badly formatted property option"]},
         {1}],
       {"Vertex", "Edge", "Face"},
       (Sequence@@#2)&]},
    {(* ==================================== Init Parameters ==================================== *)
     
     (* First we declare the simple constants for the projection *)
     SourceMesh[map] = mesh,
     
     (* Options is settable, but depends on nothing but the initial options
        (and should have nothing downstream but CorticalMapQ); 
        note that these are the graphics options only *)
     Options[map] = Cases[
       Options[CorticalMap][[All,1]],
       opt:Except[Properties] :> (opt -> OptionValue[opt]),
       {1}],
     
     (* Here, we setup all of the special map-specific parameters and translate them *)
     TranslatedCenter[map] -> CorticalMapTranslateCenter[
       SourceMesh[map],
       Center /. Join[Options[map], Options[CorticalMap]]],
     Inclusions[map] -> CorticalMapTranslateExclusions[
       SourceMesh[map],
       Method /. Join[Options[map], Options[CorticalMap]],
       TranslatedCenter[map],
       Exclusions /. Join[Options[map], Options[CorticalMap]],
       Radius /. Join[Options[map], Options[CorticalMap]]],
     TransformationFunction[map] -> CorticalMapTranslateMethod[
       SourceMesh[map],
       Method /. Join[Options[map], Options[CorticalMap]],
       TranslatedCenter[map],
       Inclusions[map],
       Radius /. Join[Options[map], Options[CorticalMap]]],
   

     (* ====================================== Properties ======================================= *)
     
     (* Note that these Properties literally just steal the mesh's properties *)
     (* #VertexProperties [private] *)
     VertexProperties[map] = With[
       {idx = Inclusions[map][[1]]},
       (#[[1]] -> #[[2]][[idx]]) & /@ VertexProperties[SourceMesh[map]]],
     (* #EdgeProperties [private] *)
     EdgeProperties[map] = With[
       {idx = Inclusions[map][[2]]},
       (#[[1]] -> #[[2]][[idx]]) & /@ EdgeProperties[SourceMesh[map]]],
     (* #FaceProperties [private] *)
     FaceProperties[map] = With[
       {idx = Inclusions[map][[3]]},
       (#[[1]] -> #[[2]][[idx]]) & /@ FaceProperties[SourceMesh[map]]],

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


     (* ======================================== Members ======================================== *)

     (* Now we have the actual vertex/face/edge data *)
     VertexList[map] -> Part[VertexList[SourceMesh[map]], Inclusions[map][[1]]],
     VertexCoordinatesTr[map] = With[
       {f = TransformationFunction[map]},
       f[VertexCoordinatesTr[SourceMesh[map]][[All, Inclusions[map][[1]]]]]],
     FaceListTr[map] -> Part[FaceListTr[SourceMesh[map]], All, Inclusions[map][[3]]],
     EdgePairsTr[map] -> Part[EdgePairsTr[SourceMesh[map]], All, Inclusions[map][[2]]],
     
     (* Untranslated / primary interface *)
     VertexCoordinates[map] := Transpose @ VertexCoordinatesTr[map],
     FaceList[map] := Transpose @ FaceListTr[map],
     EdgePairs[map] := Transpose @ EdgePairsTr[map],
     EdgeList[map] := ReplacePart[EdgePairs[map], {_,0} -> UndirectedEdge],
     
     (* pattern-matching *)
     VertexList[map, patt_] := Cases[VertexList[map], patt, {1}],
     FaceList[map, patt_] := Cases[FaceList[map], patt, {1}],
     EdgePairs[map, patt_] := Cases[EdgePairs[map], patt, {1}],
     EdgeList[map, patt_] := Cases[EdgeList[map], patt, {1}],

     (* #VertexCount *)
     VertexCount[map] -> Length @ VertexList[map],
     VertexCount[map, patt_] := Count[VertexList[map], patt, {1}],
     (* #EdgeCount *)
     EdgeCount[map] -> Length[EdgePairsTr[map][[1]]],
     EdgeCount[map, patt_] := Count[EdgePairs[map], patt, {1}],
     (* #FaceCount *)
     FaceCount[map] -> Length[FaceListTr[map][[1]]],
     FaceCount[map, patt_] := Count[FaceList[map], patt, {1}],

     (* #CorticalMapQ *)
     CorticalMapQ[map] -> Which[
       And[
         Options[map] =!= Automatic,
         Complement[
           Options[map][[All,1]],
           Cases[Options[CorticalMap][[All, 1]], Except[Properties], {1}]] != {}], Message[
             CorticalMap::badarg,
             "Unrecognized option given to CorticalMap"],
       True, True],


     (* ======================================== Indices ======================================== *)

     (* #FaceIndexArray [private] and #FaceIndex [public] *)
     FaceIndexArray[map] :> With[
       {Ft = FaceListTr[map],
        RF = Range[FaceCount[map]]},
       SparseArray @ Rule[
         Transpose @ MapThread[
           Join,
           {Ft,
            {Ft[[2]], Ft[[3]], Ft[[1]]},
            {Ft[[3]], Ft[[1]], Ft[[2]]},
            {Ft[[3]], Ft[[2]], Ft[[1]]},
            {Ft[[2]], Ft[[1]], Ft[[3]]},
            {Ft[[1]], Ft[[3]], Ft[[2]]}}],
         Join[RF, RF, RF, RF, RF, RF]]],
     FaceIndex[map, {a_Integer, b_Integer, c_Integer}] := With[
       {id = FaceIndexArray[map][[a,b,c]]},
       If[id == 0, $Failed, id]],
     FaceIndex[map, list_List] := With[
       {id = FaceIndexArray[map]},
       ReplaceAll[
         Map[Extract[id, #]&, list, {-2}],
         0 -> $Failed]],
     (* #EdgeIndexArray [private] and #EdgeIndex [public] *)
     EdgeIndexArray[map] :> With[
       {Et = EdgePairsTr[map],
        RE = Range[EdgeCount[map]]},
       SparseArray[Transpose[MapThread[Join, {Et, Reverse[Et]}]] -> Join[RE, RE]]],
     EdgeIndex[map, (List|UndirectedEdge)[a_Integer, b_Integer]] := With[
       {id = EdgeIndexArray[map][[a,b]]},
       If[id == 0, $Failed, id]],
     EdgeIndex[map, list_List] := With[
       {id = EdgeIndexArray[map]},
       ReplaceAll[
         Map[Extract[id, #]&, list, {-2}],
         0 -> $Failed]],
     (* #VertexIndexArray [private] and #VertexIndex [public] *)
     VertexIndexArray[map] :> If[VertexList[map] == Range[VertexCount[map]],
       VertexList[map],
       Normal @ SparseArray[
         VertexList[map] -> Range[VertexCount[map]],
         Max[VertexList[SourceMesh[map]]],
         0]],
     VertexIndex[map, i_Integer] := With[
       {id = If[i > Length[VertexIndexArray[map]], 0, VertexIndexArray[map][[i]]]},
       If[id == 0, $Failed, id]],
     VertexIndex[map, is_List] := With[
       {idx = VertexIndexArray[map]},
       Map[Part[idx, #]&, is, {-2}]],
     

     (* ======================================== Delays ========================================= *)

     (* These indices depend only on the face list and edge pairs and tell us to which faces a
        * vertex or an edge belongs.
        *)
     VertexEdgeList[map] :> With[
       {RE = Range[EdgeCount[map]],
        EP = VertexIndex[map, EdgePairsTr[map]]},
       Part[
         SplitBy[
           SortBy[
             Transpose[{Join[EP[[1]], EP[[2]]], Join[RE, RE]}],
             First],
           First],
         All, All, 2]],
     VertexFaceList[map] :> With[
       {T = VertexIndex[map, FaceListTr[map]],
        RT = Range[FaceCount[map]]},
       Part[
         SplitBy[
           SortBy[
             Transpose[{Join @@ T, Join[RT, RT, RT]}],
             First],
           First],
         All, All, 2]],
     EdgeFaceList[map] :> With[
       {Ft = VertexIndex[map, Transpose @ FaceList[map]],
        RT = Range[FaceCount[map]],
        eidx = EdgeIndexArray[map]},
       With[
         {edges = Extract[
            eidx,
            Transpose @ {
              Join[Ft[[1]], Ft[[2]], Ft[[3]]],
              Join[Ft[[2]], Ft[[3]], Ft[[1]]]}]},
         Part[
           SplitBy[
             SortBy[
               Transpose[{edges, Join[RT, RT, RT]}],
               First],
             First],
           All, All, 2]]],

     (* #SumOverFacesMatrix and #SumOverEdgesMatrix *)
     SumOverFacesMatrix[map] :> SparseArray[
       Transpose[{Range[3*FaceCount[map]], Join @@ FaceListTr[map]}] -> 1],
     SumOverEdgesMatrix[map] :> SparseArray[
       Transpose[{Range[2*EdgeCount[map]], Join @@ EdgePairsTr[map]}] -> 1],
     (* #SumOverFaces *)
     SumOverFacesTr[map, datat_] := Dot[
       MapThread[Join, {datat, datat, datat}],
       SumOverFacesMatrix[map]],
     SumOverFaces[map, data_] := Dot[Transpose @ Join[data, data, data], SumOverFacesMatrix[map]],
     SumOverFaceVerticesTr[map, datat_] := Dot[Join @@ datat, SumOverFacesMatrix[map]],
     SumOverFaceVertices[map, data_] := Dot[Join @@ Transpose[data], SumOverFacesMatrix[map]],
     (* #SumOverEdges *)
     SumOverEdgesTr[map, datat_] := Dot[
       MapThread[Join, {datat, datat}],
       SumOverEdgesMatrix[map]],
     SumOverEdges[map, data_] := Dot[
       Transpose @ Join[data, data],
       SumOverEdgesMatrix[map]],
     SumOverEdgeVerticesTr[map, datat_] := Dot[Join @@ datat, SumOverEdgesMatrix[map]],
     SumOverEdgeVertices[map, data_] := Dot[Join @@ Transpose[data], SumOverEdgesMatrix[map]],

     (* #VertexDegree *)
     VertexDegree[map] :> Length /@ NeighborhoodList[map],
     VertexInDegree[map] := VertexDegree[map],
     VertexOutDegree[map] := VertexDegree[map],
     VertexDegree[map, i_Integer] := VertexDegree[map][[VertexIndex[map, i]]],
     VertexInDegree[map, i_Integer] := VertexDegree[map][[VertexIndex[map, i]]],
     VertexOutDegree[map, i_Integer] := VertexDegree[map][[VertexIndex[map, i]]],


     (* ================================= Geometric Properties ================================== *)

     (* ------------------------------------------ Faces ---------------------------------------- *)

     (* #FaceAngles *)
     FaceAngleCosinesTr[map] :> CalculateFaceAngleCosinesTr[
       FaceListTr[map],
       VertexCoordinatesTr[map]],
     FaceAngleCosinesTr[map, Xtr_] := CalculateFaceAngleCosinesTr[FaceListTr[map], Xtr],
     FaceAnglesTr[map, Xtr_]       := CalculateFaceAnglesTr[FaceListTr[map], Xtr],
     FaceAngleCosines[map, Xx_]    := CalculateFaceAngleCosines[FaceList[map], Xx],
     FaceAngles[map, Xx_]          := CalculateFaceAngles[FaceList[map], Xx],
     FaceAnglesTr[map]             := ArcCos @ FaceAngleCosinesTr[map],
     FaceAngleCosines[map]         := Transpose @ FaceAngleCosinesTr[map],
     FaceAngles[map]               := Transpose @ ArcCos @ FaceAngleCosinesTr[map],

     (* #FaceAxes *)
     FaceAxesTr[map] :> CalculateFaceAxes2DTr[
       FaceListTr[map],
       VertexCoordinatesTr[map]],
     FaceAxes[map] := Transpose @ FaceAxesTr[map],
     FaceAxesTr[map, Xt_] := CalculateFaceAxesTr[FaceListTr[map], Xt],
     FaceAxes[map, Xx_] := Transpose @ CalculateFaceAxesTr[FaceListTr[map], Transpose @ Xx],

     (* #FaceBisectors *)
     FaceBisectorsTr[map] :> CalculateFaceBisectorsTr[FaceListTr[map], VertexCoordinatesTr[map]],
     FaceBisectors[map] := Transpose @ FaceBisectorsTr[map],
     FaceBisectorsTr[map, X_] := CalculateFaceBisectorsTr[FaceListTr[map], X],
     FaceBisectors[map, X_] := Transpose @ FaceBisectorsTr[map, X],

     (* #FaceCoordinates *)
     FaceCoordinatesTr[map] :> With[
       {Xt = VertexCoordinatesTr[map],
        Ft = VertedIndex[map, FaceListTr[map]]},
       Transpose @ {Xt[[All, Ft[[1]]]], Xt[[All, Ft[[2]]]], Xt[[All, Ft[[3]]]]}],
     FaceCoordinates[map] := Transpose[FaceCoordinatesTr[map], {3,2,1}],

     (* #FacePlaneCoordinates *)
     FacePlaneCoordinatesTr[map] :> CalculateFacePlaneCoordinatesTr[
       FaceListTr[map],
       VerteListTr[map],
       FaceAxesTr[map]],
     FacePlaneCoordinatesTr[map, Xt_] := CalculateFacePlaneCoordinatesTr[
       FaceListTr[map],
       Xt,
       FaceAxesTr[map, Xt]],
     FacePlaneCoordinates[map] := Transpose @ FacePlaneCoordinatesTr[map],
     FacePlaneCoordinates[map, Xx_] := Transpose @ FacePlaneCoordinatesTr[map, Transpose @ Xx],

     (* ------------------------------------------ Edges ---------------------------------------- *)

     (* #EdgeCoordinates *)
     EdgeCoordinatesTr[map] :> With[
       {Xt = VertexCoordinatesTr[map],
        EL = VertexIndex[map, EdgePairsTr[map]]},
       Transpose @ {Xt[[All, EL[[1]]]], Xt[[All, EL[[2]]]]}],
     EdgeCoordinates[map] := Transpose[EdgeCoordinatesTr[map], {3,2,1}],
     
     (* #EdgeLengths *)
     EdgeLengths[map] :> With[
       {EL = VertexIndex[map, EdgePairsTr[map]],
        Xt = VertexCoordinatesTr[map]},
       Sqrt @ Total[(Xt[[All, EL[[1]]]] - Xt[[All, EL[[2]]]])^2]],
     EdgeWeight[map] := EdgeLengths[map],
     EdgeWeight[map, e:(List|UndirectedEdge)[_Integer, _Integer]] := Part[
       EdgeLengths[map],
       EdgeIndex[map, e]],
     EdgeWeight[map, es:{{_Integer,_Integer}..}] := Part[
       EdgeLengths[map],
       EdgeIndex[map, es]],
     EdgeWeight[map, es:{(List|UndirectedEdge)[_Integer, _Integer]..}] := Part[
       EdgeLengths[map],
       EdgeIndex[map, ReplaceAll[es, UndirectedEdges -> List]]]

     (* ------------------------------------- Neighborhoods ------------------------------------- *)

     (* #NeighborhoodList *)
     NeighborhoodList[map] :> With[
       {X = VertexCoordinates[map],
        E = EdgePairsTr[map]},
       With[
         {neis = Part[
            SplitBy[
              SortBy[
                Transpose[{Join[E[[1]], E[[2]]], Join[E[[2]], E[[1]]]}],
                First],
              First],
            All, All, 2]},
         MapThread[NeighborhoodSort2DCompiled, {X, neis, X[[#]] & /@ neis}]]],
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
   

     (* ======================================= Interfaces ====================================== *)

     (* #BoundaryMeshRegion *)
     BoundaryMeshRegion[map] :> BoundaryMeshRegion[VertexCoordinates[map], Polygon[FaceList[map]]],

     (* #MeshRegion *)
     MeshRegion[map] :> MeshRegion[VertexCoordinates[map], Polygon[FaceList[map]]],
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
     AdjacencyMatrix[map, opts___]             := AdjacencyMatrix[Graph[map], opts],
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
     VertexOutDegree[map, opts___]             := VertexDegree[map]}],

  SetSafe -> True,
  Symbol -> CorticalMesh2D];

(* Like with CorticalMesh, we want to create our own version of the Clone constructor for modifying
 * things like options. Note that there are only 5 settables in the CorticalMap immutable structure:
 *  - VertexProperties
 *  - EdgeProperties
 *  - FaceProperties
 *  - SourceMesh
 *  - Options
 *)
$CorticalMapReconstructionOptions = Join[
  {SourceMesh},
  Cases[Options[CorticalMap][[All,1]], Except[Properties], {1}]];
Protect[$CorticalMapReconstructionOptions];
CorticalMap[map_?CorticalMapQ, args___Rule] := Check[
  With[
    {optsarg = (
       If[Complement[#[[All,1]], $CorticalMapReconstructionOptions] == {},
         #,
         Message[
           CorticalMap::badarg,
           StringJoin[
             "CorticalMap[] reconstructor may only be passed CorticalMesh options or SourceMesh"]]]
       )& @ {args}},
    If[Length[optsarg] == 0,
      map,
      With[
        {source = Replace[SourceMesh, optsarg],
         opts = Fold[
           Function @ With[
             {edit = Replace[#1, (Rule|RuleDelayed)[#2[[1]], _] :> #2, {1}]},
             If[SameQ[edit, #1], Append[#1, #2], edit]],
           Options[mesh],
           Select[optsarg, (#[[1]] =!= SourceMesh)&]]},
        Clone[
          map,
          Sequence @@ Flatten[
            {If[source =!= SourceMesh, SourceMesh -> source, {}],
             If[Length[opts] > 0, Options -> opts, {}]}]]]]],
  $Failed];

(* #CorticalMapQ *)
CorticalMapQ[_] := False;

(* #CorticalObjectQ *)
CorticalObjectQ[c_] := Or[CorticalMeshQ[c], CorticalMapQ[c]];

(* We want to make sure to have a nicely formatted output for a mesh or projection *)

MakeBoxes[mesh_CorticalMesh3D, form_] := RowBox[
  {"CorticalMesh3D","[",
   "<"<>ToString[VertexCount[mesh]]<>" vertices>", ",",
   "<"<>ToString[FaceCount[mesh]]<>" faces>","]"}];
MakeBoxes[mesh_CorticalMesh2D, form_] := RowBox[
  {"CorticalMesh2D","[",
   "<"<>ToString[VertexCount[mesh]]<>" vertices>", ",",
   "<"<>ToString[FaceCount[mesh]]<>" faces>","]"}];

(* Protect these functions... *)
Protect[CorticalMesh, CorticalMeshQ, Inclusions, VertexCoordinates, VertexCoordinatesTr, EdgePairs, 
        EdgePairsTr, FaceIndexArray, EdgeIndexArray, FaceIndex, EdgeIndex, EdgeCoordinates, 
        EdgeCoordinatesTr, OptionalProperties, VertexProperties, EdgeProperties, FaceProperties,
        FaceAxes, FaceAxesTr, FaceAngleCosines, FaceAngleCosinesTr, FaceAngles, FaceAnglesTr, 
        FaceCoordinates, FaceCoordinatesTr, FacePlaneCoordinates, FacePlaneCoordinatesTr, 
        FaceNormals, FaceNormalsTr, EdgeLengths, NeighborhoodList, NeighborhoodAngles, 
        NeighborhoodBisectors, NeighborhoodEdgeLengths, SourceImage, VertexEdgeList, VertexFaceList,
        EdgeFaceList];


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
        CorticalMesh[c, VertexCoordinates -> res]]]]];


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
          {Which[
             ffn === None, {},
             ffn === Automatic, {EdgeForm[], Gray, Polygon[F]},
             True, {
               EdgeForm[], Gray, 
               MapThread[
                 ffn,
                 {FaceCoordinates[mesh], F, GetProperties[FaceList], vprop[[#]]& /@ F}]}],
           If[efn === None || efn === Automatic, 
             {},
             MapThread[
               efn,
               {EdgeCoordinates[mesh], EdgePairs[mesh],
                GetProperties[EdgeList],
                Transpose[vprop[[#]]& /@ VertexIndex[mesh, EdgePairsTr[mesh]]]}]],
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
