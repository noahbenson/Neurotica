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
  * Properties, PropertyValue, SetProperty, RemoveProperty: All property functions supported by Graphs are also supported by cortical meshes. It is recommended that all data attached to nodes, faces, or edges be attached using properties. Note that the Properties option is accepted by CorticalMesh (see ?Properties) and Property wrappers in the vertex and face lists are parsed as well (see ?Property). Additionally, for programming convenience, the ReplaceAll (/.) and ReplaceRepeated (//.) operators are overloaded for meshes such that mesh /. name -> values is equivalent to SetProperty[{mesh, VertexList}, name -> values] (mesh //. data is the same, buf for edge properties), and name /. mesh is equivalent to PropertyValue[{mesh, VertexList}, name] (again, //. is for edge properties).
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
  * 3D Graphics: CorticalMesh[] accepts any option that can be passed to Graphics3D; these options will be used as default values when plotting the cortical mesh. The options may be accessed via Options[mesh] and may be changed (when cloning a mesh m via CorticalMesh[m, options...]) with Options -> {new-options}; Options -> Automatic will reset the options to the defaults accepted by Graphics3D with the following differences: Lighting -> \"Neutral\", ColorFunction -> None, ColorFunctionScaling -> False, Boxed -> False.
  * Rendering: CorticalMeshes (and CorticalMaps) are rendered by a combination of four functions: ColorFunction, EdgeRenderingFunction, FaceRenderingFunction, and VertexRenderingFunction. The ColorFunction argument is all that one usually needs to use; it determines the default color for each vertex, which is reflected in the rendering color of the faces. The edges and vertices are not rendered by default, but the vertex coloring from ColorFunction determines the VertexColors setting in the face rendering. For all of these options, the provided argument must be a function that accepts an association as its first argument, within which are keys for each property of the relevant vertex, face, or edge. For the ColorFunction and VertexRenderingFunction options, the association additionally contains the keys \"Vertex\" and \"Coordinate\". For EdgeRenderingFunction, the association additionally contains the keys \"Edge\" and \"Coordinates\". For FaceRenderingFunction, it additionally contains \"Face\" and \"Coordinates\".";
CorticalMesh::badarg = "Bad argument given to CorticalMesh constructor: `1`";
CorticalMesh3D::usage = "CorticalMesh3D is a form used to store data for 3D surface mesh objects (see CorticalMesh).";
CorticalMeshQ::usage = "CorticalMeshQ[mesh] yields True if and only if mesh is a CorticalMesh object and False otherwise.";
CorticalMesh::error = "Error in cortical mesh function `1`: `2`";

CorticalMap::usage = "CorticalMap[mesh] yields a 2D flattened cortical projection of the given cortical mesh. The following options may be given:
  * Method (default: \"Equirectangular\") specifies the projection method to be used in projecting the cortical surface to a flat map. Possible values for Method include:
    * \"Mollweide\", a Mollweide projection, parameters: Center
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
CorticalMap::error = "Error in cortical map function `1`: `2`";
SourceMesh::usage = "SourceMesh[map] yields the mesh object from which the given cortical projection, map, was constructed.";
SphericalMesh::usage = "SphericalMesh[map] yields the spherical mesh object from which coordinates were drawn for projection (if specified in the source mesh's MetaInformation; otherwise this is just the source mesh).";
NearestFaceCenters::usage = "NearestFaceCenters[mesh] yields a Nearest object f such that f[x] yields the face index whose center is nearest the point x.";
NearestFace::usage = "NearestFace[mesh, x] yields the face nearest to the point x. If x is a matrix, then the indices for all pointsi n x are yielded as a list.
The option NearestPoints may be set to true, in which case a list of {indices, points} is returned with the indices of the faces containing the points in x and the 
points that are themselves on the mesh in points.
NearestFace[mesh] yields a function f such that f[x] is equivalent to NearestFace[mesh, x].";
NearestPoints::usage = "NearestPoints is an option for the function NearestFace that indicates that {indices, points} should be returned instead of just the indices.";
CortexAddress::usage = "CortexAddress[mesh, x] yields the cortical face address of the given point or points x. A face address is a list {{a, b, c}, {t, r}} that allows the topologically equivalent point to be blooked up in any mesh with the same face topology as the given mesh. In the address, {a,b,c} are the vertex labels of the vertices of the triangle containing the point q in mesh nearest to the point x; t is the angle from the vector (a->b) to the vector (a->q) normalized by the angle from (a->b) to (a->c); and r is the distance from a to q normalized by the distance from a to the point on line-segment bc that is colinear with a and q.
CortexAddress[mesh] yields a function f such that f[x] is equivalent to CortexAddress[mesh, x].
Cortical addresses can be looked up using the function CortexLookup.";
CortexAddress::err = "Error in CortexAddress: `1`";
CortexLookup::usage = "CortexLookup[mesh, address] yields the point in the given mesh that is specified by the given address. See also CortexAddress.
CortexLookup[mesh] yields a function f such that f[x] is equivalent to CortexLookup[mesh, x].";
CortexLookup::err = "Error in CortexLookup: `1`";

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

IndexedEdgePairsTr::usage = "IndexedEdgePairsTr[mesh] is identical to Transpose@IndexedEdgePairs[mesh].";
IndexedEdgePairs::usage = "IndexedEdgePairs[mesh] is identical to VertexIndex[mesh, EdgePairs[mesh]].";
IndexedFaceListTr::usage = "IndexedFaceListTr[mesh] is identical to Transpose@IndexedFaceList[mesh].";
IndexedFaceList::usage = "IndexedFaceList[mesh] is identical to VertexIndex[mesh, FaceList[mesh]].";

VertexFaceMatrix::usage = "VertexFaceMatrix[M] yields a SparseArray matrix S such that each element S[[m*k + i,j]] of S is equal to 0 if vertex j is the k'th part of face i and 0 if not (where m is the number of faces in the mesh M).";
VertexEdgeMatrix::usage = "VertexEdgeMatrix[M] yields a SparseArray matrix S such that each element S[[m*k + i,j]] of S is equal to 0 if vertex j is the k'th part of edge i and 0 if not (where m is the number of edges in the mesh M).";
EdgeFaceMatrix::usage = "VertexFaceMatrix[M] yields a SparseArray matrix S such that each element S[[m*k + i,j]] of S is equal to 0 if edge j is the kth part of face i and 0 if not (where m is the number of faces in the mesh M).";

SumOverFacesTr::usage = "SumOverFacesTr[M, Q] yields the 3 x n result of summing over the faces given in the matrix Q whose first dimension must be equal to 3 and whose last dimension must be equal to the number of faces in cortical mesh or map M.";
SumOverFaceVerticesTr::usage = "SumOverFaceVerticesTr[M, Q] yields the vector result of summing over the vertices in the faces given in the matrix Q whose dimensions must be {3,m} where m must be equal to the number of faces in cortical mesh or map M; the result is a length n vector where n is the number of vertices.";
SumOverFaceVertices::usage = "SumOverFaceVertices[M, Q] is equivalent to Transpose @ SumOverFaceVerticesTr[M, Transpose @ Q].";
SumOverEdgesTr::usage = "SumOverEdgesTr[M, Q] yields the 3 x n result of summing over the edges given in the matrix Q whose first dimension must be equal to 2 and whose last dimension must be equal to the number of edges in cortical mesh or map M.";
SumOverFaces::usage = "SumOverFaces[M, Q] yields the result of summing over the faces given in the matrix Q whose last dimension must be equal to 3 and whose first dimension must be equal to the number of faces in cortical mesh or map M.";
SumOverEdges::usage = "SumOverEdges[M, Q] yields the result of summing over the edges given in the matrix Q whose last dimension must be equal to 2 and whose first dimension must be equal to the number of faces in cortical mesh or map M.";
SumOverEdgesDirected::usage = "SumOverEdgesDirected[M, Q] yields the equivalent of SumOverEdges[M, Q] except that the values in Q are subtracted from the second element of each edge instead of added to it.";
SumOverEdgesDirectedTr::usage = "SumOverEdgesDirectedTr[M, Q] yields the equivalent of SumOverEdgesTr[M, Q] except that the values in Q are subtracted from the second element of each edge instead of added to it.";

NeighborhoodList::usage = "NeighborhoodList[s] yields a list of length N (where N is the number of vertices in s) of the neighboring vertices of each vertex; each entry k of the list is a list of the integer id's of the neighbors of the kth vertex. The neighbor id's are sorted such that they are listed in a counter-clockwise order around vertex k starting from the x-axis. The argument s may be either a map or a surface.";
NeighborhoodAngles::usage = "NeighborhoodAngles[mesh] yields a list of the angles between the nodes in the NeighborhoodList; these angles are in the same order as the nodes in NeighborhoodList such that the first angle in a neighborhood is between the first vertex in the neighborhood, the central vertex, and second vertex in the neighborhood and the last angle in the neighborhood angles list is the angle between the first vertex in the neighborhood list, the center vertex, and the last vertex in the neighborhood list.
NeighborhoodAngles[s, X] yields the neighborhood angles for s if the vertices of s were replaced with the vertices in the list X.";
NeighborhoodBisectors::usage = "NeighborhoodBisectors[mesh] yeilds a list of the vectors that bisect each of the angles in NeighborhoodAngles[mesh].
NeighborhoodBisectors[s, X] yields the neighborhood bisecting vectors for the points given by the coordinate matrix X.";
NeighborhoodEdgeLengths::usage = "NeighborhoodEdgeLengths[s] yields a list of the edge lengths for the neighborhood of each vertex in the surface or map s. The results are in the same order as NeighborhoodList[s] such that for the neighborhood list L, and the neighborhood edge length list G, G[[i,j]] is the length of the edge from vertex i to the vertex L[[i,j]].";
SumOverNeighbors::usage = "SumOverNeighbors[mesh] yields a sparse matrix of size n x n where n is the number of vertices in the given mesh and whose entries are 0's everywhere except when, for row i and column j, vertex i and vertex j are neighbors.
SumOverNeighbors[mesh, vector] yields the sum over neighbors, for each vertex, of the given vector.
SumOverNeighbors[mesh, matrix] yields the sum over neighbors, for each column of the given matrix.";

FaceIndex::usage = "FaceIndex[mesh, f] yields the index in the FaceList[mesh] of the face f. The vertices in f may be in any order.";
VertexFaceList::usage = "VertexFaceList[mesh] yields a list of the faces (in the same order as VertexList[mesh]) that each vertex is a member of.
VertexFaceList[mesh, vertex] yields a list of just the faces that the given vertex is a member of.";
VertexEdgeList::usage = "VertexEdgeList[mesh] yields a list of the edges (in the same order as VertexList[mesh]) that each vertex is a member of.
VertexEdgeList[mesh, vertex] yields a list of just the edges that the given vertex is a member of.";
EdgeFaceList::usage = "EdgeFaceList[mesh] yields a list of, for each edge (in the same order as EdgeList and EdgePairs), the faces to which that edge belongs.";

MapBoundaryVertexList::usage = "MapBoundaryVertexList[map] yields a list of the vertices, in counter-clockwise order, that form the outer boundary of the given cortical map.";
MapBoundaryEdgeList::usage = "MapBoundaryEdgeList[map] yields a list of the edges, in counter-clockwise order, that form the outer boundary of the given cortical map.";
MapBoundaryEdgePairs::usage = "MapBoundaryEdgePairs[map] yields a list of the edge pairs, in counter-clockwise order, that form the outer boundary of the given cortical map.";
MapBoundaryFaceList::usage = "MapBoundaryFaceList[map] yields a list of the faces, in counter-clockwise order, that form the outer boundary of the given cortical map.";
MapBoundaryEdgePairsTr::usage = "MapBoundaryEdgePairsTr[map] is equivalent to Transpose @ MapBoundaryEdgePairs[map].";
MapBoundaryFaceListTr::usage = "MapBoundaryFaceListTr[map] is equivalent to Transpose @ MapBoundaryFaceList[map].";

SourceImage::usage = "SourceImage[mesh] yields the source volume of the given mesh, if the given mesh has specified a volume; otherwise None is yielded.";

FaceRenderingFunction::usage = "FaceRenderingFunction is an option that can be given to CorticalMesh or CortexPlot3D, which specifies how to render the faces of a cortical mesh. The function is organized as with VertexRenderingFunction and EdgeRenderingFunction; the arguments are in an association passed as the first argument; each property for the face is listed in the association as well as the keys \"Coordinates\" and \"Face\" which give the cartesian coordinates of the vertcices in the face and the vertices in the face, respectively.";

Inclusions::usage = "Inclusions[map] yields a list of the {vertexIndives, edgeIndices, faceIndices} of the vertices, edges, and faces from the source mesh of the given map that are included in map.";

CortexPlot3D::usage = "CortexPlot3D[mesh] yields a 3D Graphics form for the given CorticalMesh3D mesh. All options available to Graphics3D may be passed to CortexPlot3D. Note that some of the default options for Graphics3D have been altered in CortexPlot3D, and 3D graphics options that have been attached to the mesh will be used as well. See ?CorticalMesh for more details.";

CortexPlot::usage = "CortexPlot[mesh] yields a Graphics form for the given CorticalMesh2D or CorticalMesh3D mesh. If the given mesh is a 3D mesh, then the options accepted by CorticalMap are used to create the projection (Method, Center, etc.). All options available to Graphics may be passed to CortexPlot. Note that some of the default options for Graphics have been altered in CortexPlot, and 2D graphics options that have been attached to the mesh will be used as well. See ?CorticalMap for more details.

If a CorticalMesh3D is passed to CortexPlot, the function attempts to convert it into a map. If the mesh's MetaInformation includes a rule with the head \"CorticalMap\", then the tail must be a list and must contain valid options which are passed to CorticalMap.";

InterpretVertexColor::usage = "InterpretVertexColor[instruction, vertexData] yields the color that corresponds to the given color instruction for the given vertex data; this is used by CortexPlot and CortexPlot3D to translate instructions into actual color value when the ColorFunction option is something such as ColorFunction -> (\"Curvature\"&) or ColorFunction -> {{\"PolarAngle\", \"Curvature\"}, 0.65}. The vertexData should be an association of the vertex's properties from, i.e., VertexData[mesh].";

CorticalCurvatureColor::usage = "CorticalCurvatureColor[c] yields the appropriate color for cortical curvature c in a CortexPlot or CortexPlot3D; c may be a list or a single value.";
CorticalCurvatureVertexColors::usage = "CorticalCurvatureVertexColors[m] yields the colors for each vertex in the given mesh m according to CorticalCurvatureColor[c] for the curvature c of each vertex; if there is no curvature proeprty defined, then Gray is used for each vertex.";

ColorCortex::usage = "ColorCortex[forms...] yields a function that colors the cortical surface according to the given forms. Each forms is evaluated as in a CompoundExpression form in order to indicate some aspect of the coloring of the vertices in the cortical object that is being colored, based on the final value: 
  * CorticalColorSchema: If a registered cortical color schema is given (see CorticalColorData), then it is attempted as a coloring method.
  * RGBColor[...], Hue[...], or any form x for which ColorQ[x] yields True: assigns the given color to the vertex in question; note that opacity is interpreted and implemented with a blending operation.
In all cases and in all forms, any occurance of Slot[name] will result in that slot being replaced with the value of the property with the given name for the vertex being colored. For example, #Curvature may be used to refer to the Curvature property. Note that #Curvature is actually Slot[\"Curvature\"] rather than Slot[Curvature]; symbol names are interpreted as strings in the property list in such cases. The slots #1 and #2 are bound to the vertex id and vertex coordinate, respectively.
If no valid color can be found for a form, then the color RGBColor[1,1,1,0] (transparent) is used.";
CorticalColorData::usage = "CorticalColorData[name] yields the cortical coloring instructions registered to the given name. This will always be a CorticalColorSchema object, which, among other things, serves as a function that, given a vertex id, vertex coordinate, and set of properties, yields a color.
The form CorticalColorData[name] = schema; may be evaluated to construct a CorticalColorData object.  The schema may be a CorticalColorSchema object or any instruction that can be used to establish a cortical color schema (see CorticalColorSchema for more information).";
CorticalColorSchema::usage = "CorticalColorSchema[...] is a form that declares a cortical color schema object. CorticalColorSchema[] may be called with the following arguments to construct a color schema object:
CorticalColorSchema[property -> {range, colors}] indicates the vertices should be colored according to the given colors, blended over the given range of the given property.
CorticalColorSchema[property -> function] indicates that the given function should be passed the given property value and will return a color.
CorticalColorSchema[All -> function] indicates that the given function should accept three arguments: the vertex id, the vertex coordinate, and the vertex property list; the fun
ction must return a color or $Failed or None.";

EdgePropertyAssociation::usage = "EdgePropertyAssociation[mesh] yields an association whose keys are the property names of the edges of the given mesh and whose values are the lists of properties for each edge.";
FacePropertyAssociation::usage = "FacePropertyAssociation[mesh] yields an association whose keys are the property names of the faces of the given mesh and whose values are the lists of properties for each face.";

VertexDataset::usage = "VertexDataset[mesh] yields a dataset of the vertex properties for the given mesh.";
EdgeDataset::usage = "EdgeDataset[mesh] yields a dataset of the edge properties for the given mesh.";
FaceDataset::usage = "FaceDataset[mesh] yields a dataset of the face properties for the given mesh.";

VertexPropertyList::usage   = "VertexPropertyList[mesh, prop] is equivalent to PropertyList[{mesh, VertexList}, prop].";
EdgePropertyList::usage     = "EdgePropertyList[mesh, prop] is equivalent to PropertyList[{mesh, EdgeList}, prop].";
FacePropertyList::usage     = "FacePropertyList[mesh, prop] is equivalent to PropertyList[{mesh, FaceList}, prop].";
VertexPropertyValues::usage = "VertexPropertyValues[mesh, prop] is equivalent to PropertyValue[{mesh, VertexList}, prop].";
EdgePropertyValues::usage   = "EdgePropertyValues[mesh, prop] is equivalent to PropertyValue[{mesh, EdgeList}, prop].";
FacePropertyValues::usage   = "FacePropertyValues[mesh, prop] is equivalent to PropertyValue[{mesh, FaceList}, prop].";
SetVertexProperties::usage  = "SetVertexProperties[mesh, prop] is equivalent to SetProperty[{mesh, VertexList}, prop].";
SetEdgeProperties::usage    = "SetEdgeProperties[mesh, prop] is equivalent to SetProperty[{mesh, EdgeList}, prop].";
SetFaceProperties::usage    = "SetFaceProperties[mesh, prop] is equivalent to SetProperty[{mesh, faceList}, prop].";
RemoveVertexProperty::usage = "RemoveVertexProperty[mesh, prop] is equivalent to RemoveProperty[{mesh, VertexList}, prop].";
RemoveEdgeProperty::usage   = "RemoveEdgeProperty[mesh, prop] is equivalent to RemoveProperty[{mesh, EdgeList}, prop].";
RemoveFaceProperty::usage   = "RemoveFaceProperty[mesh, prop] is equivalent to RemoveProperty[{mesh, FaceList}, prop].";

EdgeVertexProperties::usage = "EdgeVertexProperties[mesh] yields a 2 x n matrix (n = number of edges in the given mesh) in which row i gives the pair of associations of vertex properties for the two vertices in the i'th edge of the mesh.";
FaceVertexProperties::usage = "FaceVertexProperties[mesh] yields a 3 x n matrix (n = number of faces in the given mesh) in which row i gives the three associations of vertex properties for the vertices in the i'th face of the mesh.";

MapVertices::usage = "MapVertices[f, mesh] is equivalent to Map[f, Normal@VertexDataset[mesh]].";
MapEdges::usage = "MapEdges[f, mesh] is equivalent to Map[f, Normal@EdgeDataset[mesh]].";
MapFaces::usage = "MapFaces[f, mesh] is equivalent to Map[f, Normal@FaceDataset[mesh]].";

SelectVertices::usage = "SelectVertices[mesh, f] is equivalent to Pick[VertexList[mesh], Normal@VertexDataset[map], x_ /; TrueQ[f[x]]].";
SelectEdges::usage = "SelectEdges[mesh, f] is equivalent to Pick[EdgeList[mesh], Normal@EdgeDataset[map], x_ /; TrueQ[f[x]]].";
SelectFaces::usage = "SelectFaces[mesh, f] is equivalent to Pick[FaceList[mesh], Normal@FaceDataset[map], x_ /; TrueQ[f[x]]].";

SelectVertexIndices::usage = "SelectVertexIndices[mesh, f] is equivalent to VertexIndex[mesh, SelectVertices[mesh, f]].";
SelectEdgeIndices::usage = "SelectEdgeIndices[mesh, f] is equivalent to EdgeIndex[mesh, SelectEdges[mesh, f]].";
SelectFaceIndices::usage = "SelectFaceIndices[mesh, f] is equivalent to FaceIndex[mesh, SelectFaces[mesh, f]].";
SelectIndexedEdges::usage = "SelectIndexedEdges[mesh, f] is equivalent to VertexIndex[mesh, SelectEdges[mesh, f]].";
SelectIndexedFaces::usage = "SelectIndexedFaces[mesh, f] is equivalent to VertexIndex[mesh, SelectFaces[mesh, f]].";

Reproject::usage = "Reproject[map, X] yields a map identical to the given map except that it reprojects its coordinates from the alternate coordinate list for the original mesh, given by X. If X is instead a mesh with the same number of elements as the original mesh, then its coordinates are used.";
ReporjectTr::usage = "ReprojectTr[map, Xt] is equivalent to Reproject[map, Transpose[Xt]].";
Reproject::badarg = "Bad argument given to Reproject: `1`";

InverseProject::usage = "InverseProject[map] attempts to reverse the projection of map back to the 3D space of the SourceMesh[map].
InverseProject[map, X] performs the same projection, but using the replacement map coordinates given in X.";
InverseProjectTr::usage = "InverseProjectTr[map] is equivalent to Transpose @ InverseProject[map].
InverseProjectTr[map, Xt] is equivalent to Transpose @ InverseProject[map, Transpose[Xt]].";

InverseProjectVectors::usage = "InverseProjectVectors[map, U] yields a list of 3D vectors, one for each vertex in the given cortical map, corresponding to the 2D vectors given in the matrix U. U is expected to have a single 2D vector for each vertex in map.";
InverseProjectVectorsTr::usage = "InverseProjectVectorsTr[map, Ut] is equivalent to Transpose @ InverseProjectVectors[map, Transpose[Ut]].";

CortexResample::usage = "CortexResample[surf1, surf2] yields a cortical equivalent to the CorticalMesh object surf2 but such that the field of the surface has been resampled from the cortical mesh object surf1.
CortexResample[surf2 -> surf1] is equivalent to CortexResample[surf1, surf2].
The following options may be provided:
  * Method: a Method option may specify \"Interpolation\" (default) for trilinear interpolation  or \"Nearest\" for nearest-neighbor interpolation. Note that trilinear interpolation of a mesh gives the interpolation of the nearest point in the mesh to the queried point.
  * Properties: a Properties argument specifies that the given property should be resampled; a list of properties may also be given, or All. If no property is given then All is the default value.
  * Indeterminate: a value to insert into the interpolated values in the case that non-numeric values were found nearby";
CortexResample::badarg = "Bad argument given to CortexResample: `1`";
CortexResampleOnto::usage = "CortexResampleOnto[args...] is identical to CortexResample[args...] but instead of returning the new properties, it places the resamples properties onto the new mesh and returns the updated mesh.";

CorticalLabelQ::usage = "CorticalLabelQ[mesh, label] yields True if and only if label is a valid label for the given mesh; otherwise yields false. A label is valid for the given mesh if any of the following conditions are met:
  * label is a list of valid vertices in mesh;
  * label is a list or sparse array with a length equal to the vertex count of mesh in which every value is in the set {True, False, 1, 0, 1.0, 0.0};
  * mesh is 2D and label is a list of edge pairs or edges that form a cycle;
  * mesh is 3D and label is a list whose first element is a vertex u and whose remaining elements are a list of edge pairs or edges that form a cycle such that u is not in the list of edges; in this case, u is by definition on the inside of the label;
  * mesh is 2D and label is a BoundaryMeshRegion whose RegionDimension is 2;
  * mesh is 3D and label is a BoundaryMeshRegion whose RegionDimension is 3.";
CorticalLabelMask::usage = "CorticalLabelMask[mesh, label] yields a mask of the vertices in the given mesh that are included in the label. The mask returned is in the form or a SparseArray or list with 0s indicating vertices that are not included and 1s indicating vertices that are. If label is not a valid cortical label, then $Failed is yielded and a message is generated.";
CorticalLabelMask::notlab = "Label argument given to CorticalLabelMask is not a valid cortical label.";
CorticalLabelVertexList::usage = "CorticalLabelVertexList[mesh, label] yields a list of the vertices in the given mesh that are included in the label. If label is not a valid cortical label, then $Failed is yielded and a message is generated.";
CorticalLabelVertexList::notlab = "Label argument given to CorticalLabelVertexList is not a valid cortical label.";

LabelVertexCoordinates::usage = "LabelVertexCoordinates[sub, mesh, hemi, name] yields the vertex coordinates for the given subject, mesh, and hemisphere of the vertices that lie in the label with the given name.
Note that if you are defining a new subject modality, you should define LabelVertexList[] instead of LabelVertexCoordinates.";
LabelVertexCoordinatesTr::usage = "LabelVertexCoordinatesTr[sub, mesh, hemi, name] is equivalent to Transpose[LabelVertexCoordinates[sub, mesh, hemi, name]].
Note that if you are defining a new subject modality, you should define LabelVertexList[] instead of LabelVertexCoordinatesTr.";

LabelEdgePairsTr::usage = "LabelEdgePairsTr[sub, hemi, name] yields the equivalent of Transpose @ LabelEdgePairs[sub, hemi, name].";
LabelEdgePairs::usage = "LabelEdgePairs[sub, hemi, name] yields a list of the edge pairs that compose the label with the given name.";
LabelEdgeList::usage = "LabelEdgeList[sub, hemi, name] is equivalent to LabelEdgePairs[sub, hemi, name] except that it yields UndirectedEdge's instead of lists of vertex pairs.";

LabelFaceListTr::usage = "LabelFaceListTr[sub, hemi, name] yields the equivalent of Transpose @ LabelFaceList[sub, hemi, name].";
LabelFaceList::usage = "LabelFaceList[sub, hemi, name] yields the list of faces that are part of the label with the given name in the given subject and hemisphere.";

LabelBoundaryVertexList::usage = "LabelBoundaryVertexList[sub, hemi, name] yields a list of the vertices, in counter-clockwise order, that are part of the boundary between the given label and the outside of the label.";
LabelBoundaryVertexCoordinates::usage = "LabelBoundaryVertexCoordinates[sub, mesh, hemi, name] yields a list of the coordinates of the vertices, in counter-clockwise order, that appear in on the boundary of the given label in the given subject, mesh, and hemisphere.";
LabelBoundaryVertexCoordinatesTr::usage = "LabelBoundaryVertexCoordinatesTr[sub, mesh, hemi, name] is equivalent to Transpose @ LabelBoundaryVertexCoordinates[sub, mesh, hemi, name].";

LabelBoundaryEdgeList::usage = "LabelBoundaryEdgeList[sub, hemi, name] yields a list of the edges that form an approximate boundary to the region with the given label name for the given subject and hemisphere.
Note that if you are defining a new subject modality, then you only need to define the LabelVertexList function.";
LabelBoundaryEdgePairsTr::usage = "LabelBoundaryEdgePairsTr[sub, hemi, name] yields a list equivalent to Transpose[LabelBoundaryEdgePairs[sub, hemi]].";
LabelBoundaryEdgePairs::usage = "LabelBoundaryEdgePairs[sub, hemi, name] yields a list of the edge pairs rather than the edges themselves that are returned by LabelEdgeList.";
LabelBoundaryEdgePairs::badarg = "Bad argument given to LabelBoundary function: `1`";

OccipitalPole::usage = "OccipitalPole[subject, mesh, hemisphere] is usually defined by subject modalities (e.g., FreeSurferSubject[]) such that the function yields the vertex coordinate for the occipital pole in the particular mesh and hemisphere requested.
Note that if you define a new subject modality, then defining OccipitalPoleIndex[] for the subject should be sufficient.";

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

(* #CorticalMapUnorientForMesh
 * Performs the inverse of the above function.
 *)
CorticalMapUnorientForMesh[Xt_, center_] := With[
  {RMain = If[center[[1]] === None || center[[1]] === Automatic, 
    IdentityMatrix[3],
    RotationMatrix[{center[[1]], {1,0,0}}]]},
  Dot[
    Inverse@Dot[
      If[center[[2]] === None || center[[2]] === Automatic,
        IdentityMatrix[3], 
        With[
          {orientPt = Dot[RMain, center[[2]]]},
          RotationMatrix[-ArcTan[orientPt[[1]], orientPt[[2]]], {1,0,0}]]],
      RMain],
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
    "mollweide"|"equirectangular"|"mercator", Pick[
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
        {incl = If[prad === Full || prad === All || prad === Automatic,
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
               Intersection[
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
 * Yields a list of {translation function, inverse translation function} given the Method argument
 * and the translated parameters: Center, vertex-exclusions, projection area, and projection radius.
 *)
CorticalMapTranslateMethod[mesh_, method_, center_, incl_, prad_] := Check[
  With[
    {projFn = Switch[
       ToLowerCase[method],
       (* Here we actually define the functions for transformation;
          these should accept a list of vertices (already centered such that the center lies at 
          (1,0,0) and the orient point lies in the <positive X, Y> half-plane) and should return
          the 2D coordinates (while ignoring the orient point). *)
       "mollweide", {
         Function[{X},
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
         None},
       "equirectangular", With[
         {radvar = {Mean[#], Variance[#]}& @ ColumnNorms[VertexCoordinatesTr[mesh]]},
         {Function[{X},
            With[
              {S = ConvertCoordinates[Transpose@X, Cartesian -> {Longitude, Latitude}]},
              Transpose[radvar[[1]] * S / Pi]]],
          Function[{S},
            With[
              {X = ConvertCoordinates[
                 Transpose[S * Pi / radvar[[1]]], 
                 {Longitude, Latitude} -> Cartesian]},
              Transpose[X * radvar[[1]]]]]}],
       "mercator", With[
         {radvar = {Mean[#], Variance[#]}& @ ColumnNorms[VertexCoordinatesTr[mesh]]},
         {Function[{X},
            With[
              {S = Transpose@ConvertCoordinates[Transpose@X, Cartesian -> {Longitude, Latitude}]},
              With[
                {sinPhi = Sin[S[[2]]]},
                radvar[[1]] * {
                  S[[1]],
                  (*0.5*Log[(1 + sinPhi) / (1 - sinPhi)]*)
                  Log @ Tan[0.25*Pi + 0.5*S[[2]]]}]]],
          Function[{S},
            radvar[[1]] * Transpose@ConvertCoordinates[
              Transpose[{S[[1]] * radvar[[1]], 2.0 * ArcTan[Exp[S[[2]] / radvar[[1]]]]}],
              {Longitude, Latitude} -> Certesian]]}],
       "orthographic", With[
         {radvar = {Mean[#], Variance[#]}& @ ColumnNorms[VertexCoordinatesTr[mesh]]},
         {Function[{Xt}, Xt[[2;;3]]],
          Function[{Xt}, radvar[[1]] * Prepend[#, Sqrt[1.0 - Total[#^2]]]&[Xt / radvar[[1]]]]}],
       "polarstretching", {
         Function[{X},
           With[
             {polar = Transpose @ ConvertCoordinates[
                Transpose @ X[[{2,3,1}]],
                Cartesian -> {Longitude, SphericalPolarAngle}]},
             {# * Cos[polar[[1]]], # * Sin[polar[[1]]]}& @ Tan[polar[[2]] / 2]]],
         None},
       "graph", {
         With[
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
         None},
       _, Message[CorticalMap::badarg, "Could not recognize projection type"]],
     RMtx = RotationMatrix[{center[[1]], {1,0,0}}]},
    With[
      {orientFn = If[center[[2]] === Automatic,
         Function[RotationMatrix[{First@Eigenvectors[Covariance[Transpose@#],1], {1,0}}] . #],
         With[
           {orientRMtx = RotationMatrix[
              {projFn[[1]][CorticalMeshOrientForMap[List /@ center[[2]], center]][[All, 1]],
               {1,0}}]},
           Function[orientRMtx . #]]],
       unorientFn = If[center[[2]] === Automatic,
         None,
         With[
           {unorientRMtx = Inverse @ RotationMatrix[
              {projFn[[1]][CorticalMeshOrientForMap[List /@ center[[2]], center]][[All, 1]],
               {1,0}}]},
           Function[unorientRMtx . #]]]},
      {Function[
         If[Length[#] == 3,
           orientFn @ projFn[[1]][CorticalMeshOrientForMap[#, center]],
           Transpose @ orientFn[projFn[[1]][CorticalMeshOrientForMap[Transpose @ #, center]]]]],
       If[unorientFn === None || projFn[[2]] === None,
         None,
         Function[
           If[Length[#] == 3,
             CorticalMapUnorientForMesh[projFn[[2]][unorientFn[#]], center],
             Transpose @ CorticalMapUnorientForMesh[
               projFn[[2]][unorientFn[#]],
               center]]]]}]],
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
   VertexColors -> Automatic,
   Lighting -> "Neutral",
   Boxed -> False,
   VertexRenderingFunction -> None,
   EdgeRenderingFunction -> None,
   FaceRenderingFunction -> Automatic}];
Protect[$CortexPlot3DOptions];

(* These properties of meshes are read-only: *)
$VertexReadOnlyProperties = {"VertexIndex", "VertexLabel", "VertexNormals"};
$EdgeReadOnlyProperties = {"VertexIndex", "VertexLabel", "VertexCoordinates",
                           "EdgeWeight", "EdgeLengths", "VertexProperties"};
$FaceReadOnlyProperties = {"VertexIndex", "VertexLabel", "VertexCoordinates",
                           "FaceNormals", "VertexProperties"};

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
       Transpose[{Range[3*FaceCount[mesh]], Join @@ VertexIndex[mesh, FaceListTr[mesh]]}] -> 1],
     SumOverEdgesMatrix[mesh] :> SparseArray[
       Transpose[{Range[2*EdgeCount[mesh]], Join @@ VertexIndex[mesh, EdgePairsTr[mesh]]}] -> 1],
     SumOverNeighborsMatrix[mesh] :> SparseArray[
       Join@@MapThread[
         Thread[{#1, #2}]&,
         {Range@VertexCount[mesh], VertexIndex[mesh, NeighborhoodList[mesh]]}] -> 1,
       {VertexCount[mesh], VertexCount[mesh]},
       0],


     (* ======================================= Settables ======================================= *)

     (* All data is actually stored transposed for efficiency *)
     VertexCoordinatesTr[mesh] = Transpose @ cleaned[[1]],
     FaceListTr[mesh]          = Transpose @ cleaned[[3]],
     Options[mesh]             = Cases[
       Options[CorticalMesh][[All,1]],
       opt:Except[Properties] :> (opt -> OptionValue[opt]),
       {1}],
     MetaInformation[mesh]    := Fold[
       Replace,
       MetaInformation,
       {Options[mesh], MetaInformation -> {}}],
     CortexToRASMatrix[mesh]  := Replace[
       "CortexToRASMatrix",
       Append[MetaInformation[mesh], "CortexToRASMatrix" -> None]],


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
       Length@Intersection[VertexProperties[mesh][[All,1]], $VertexReadOnlyProperties] > 0, Message[
         CorticalMesh::badarg,
         "Unable to set read-only vertex property"],
       Length@Intersection[EdgeProperties[mesh][[All,1]], $EdgeReadOnlyProperties] > 0, Message[
         CorticalMesh::badarg,
         "Unable to set read-only edge property"],
       Length@Intersection[FaceProperties[mesh][[All,1]], $FaceReadOnlyProperties] > 0, Message[
         CorticalMesh::badarg,
         "Unable to set read-only face property"],
       And[
         Options[mesh] =!= Automatic,
         Complement[
           Options[mesh][[All,1]], 
           Cases[Options[CorticalMesh][[All, 1]], Except[Properties], {1}]] != {}], Message[
           CorticalMesh::badarg,
           "Unrecognized option given to CorticalMesh"],
       True, True],



     (* ======================================== Delays ========================================= *)

     (* These indices are simple conversions of the vertex labels into vertex indices: *)
     IndexedEdgePairsTr[mesh] :> VertexIndex[mesh, EdgePairsTr[mesh]],
     IndexedEdgePairs[mesh] := Transpose@IndexedEdgePairsTr[mesh],
     IndexedEdgeList[mesh] := UndirectedEdge@@@IndexedEdgePairs[mesh],
     IndexedFaceListTr[mesh] :> VertexIndex[mesh, FaceListTr[mesh]],
     IndexedFaceList[mesh] :> Transpose@IndexedFaceListTr[mesh],

     (* We keep track of the nearest vertices *)
     Nearest[mesh] :> Nearest[VertexCoordinates[mesh] -> VertexList[mesh]],
     NearestFaceCenters[mesh] :> Nearest@Rule[
       Transpose@Total@Transpose@FaceCoordinatesTr[mesh] / 3.0,
       Automatic],
    
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
       {Ft = FaceListTr[mesh],
        RT = Range[FaceCount[mesh]]},
       With[
         {edges = {
            Join[
              Range[EdgeCount[mesh]],
              Extract[
                EdgeIndexArray[mesh],
                Transpose @ {
                  Join[Ft[[1]], Ft[[2]], Ft[[3]]],
                  Join[Ft[[2]], Ft[[3]], Ft[[1]]]}]],
            Join[ConstantArray[0, EdgeCount[mesh]], RT, RT, RT]}},
         With[
           {idx = SplitBy[
              Transpose @ edges[[All, Ordering[edges[[1]]]]],
              First]},
           idx[[All, 2;;All, 2]]]]],

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

     EdgeVertexProperties[mesh] :> With[
       {pp = Properties[mesh],
        p = Normal@VertexDataset[mesh]},
       Transpose[p[[#]]& /@ IndexedEdgePairsTr[mesh]]],
     FaceVertexProperties[mesh] :> With[
       {pp = Properties[mesh],
        p = Normal@VertexDataset[mesh],
        F = IndexedFaceListTr[mesh]},
       Transpose[p[[#]]& /@ F]],
     
     (* Vertex, Edge, and Face Property associations and datasets *)
     VertexPropertyAssociation[mesh] :> With[
       {pp = Properties[mesh]}, (* this is a dependency, basically *)
       Association@Table[
         prop -> VertexPropertyValues[mesh, prop],
         {prop, PropertyList[{mesh, VertexList}]}]],
     VertexDataset[mesh] :> With[
       {cols = VertexPropertyList[mesh],
        n = VertexCount[mesh],
        pp = Properties[mesh]},
       Dataset@Map[
         Function@Association@Thread[cols -> #],
         Transpose@Replace[
           Map[Function@Normal@VertexPropertyValues[mesh, #], cols],
           Except[_List] :> ConstantArray[$Failed, n],
           {1}]]],
     EdgePropertyAssociation[mesh] :> With[
       {dep1 = Properties[mesh],
        dep2 = EdgeVertexProperties[mesh]},
       Association@Table[
         prop -> EdgePropertyValues[mesh, prop],
         {prop, EdgePropertyList[mesh]}]],
     EdgeDataset[mesh] :> With[
       {cols = EdgePropertyList[mesh],
        n = EdgeCount[mesh],
        dep1 = Properties[mesh],
        dep2 = EdgeVertexProperties[mesh]},
       Dataset@Map[
         Function@Association@Thread[cols -> #],
         Transpose@Replace[
           Map[Function@EdgePropertyValues[mesh, #], cols],
           Except[_List] :> ConstantArray[$Failed, n],
           {1}]]],
     FacePropertyAssociation[mesh] :> With[
       {dep1 = Properties[mesh],
        dep2 = FaceVertexProperties[mesh]},
       Association@Table[
         prop -> FacePropertyValues[mesh, prop],
         {prop, FacePropertyList[mesh]}]],
     FaceDataset[mesh] :> With[
       {cols = FacePropertyList[mesh],
        n = FaceCount[mesh],
        dep1 = Properties[mesh],
        dep2 = FaceVertexProperties[mesh]},
       Dataset@Map[
         Function@Association@Thread[cols -> #],
         Transpose@Replace[
           Map[Function@FacePropertyValues[mesh, #], cols],
           Except[_List] :> ConstantArray[$Failed, n],
           {1}]]],


     (* ======================================== Methods ======================================== *)

     (* Extension of Options... *)
     Options[mesh, opt_] := Replace[opt, Options[mesh]],

     (* A few simple extensions of the immediates that count things *)
     VertexCount[mesh, patt_] := Count[VertexList[mesh], patt, {1}],
     FaceCount[mesh, patt_] := Count[FaceList[mesh], patt, {1}],
     EdgeCount[mesh, patt_] := Count[EdgePairs[mesh], patt, {1}],
     
     (* EdgeList is here just for the compatibility with graph-like operations (use EdgePairs) *)
     EdgeList[mesh] := UndirectedEdge@@@EdgePairs[mesh],
     EdgeList[mesh, patt_] := Cases[EdgeList[mesh], patt, {1}],
     
     (* Untransposed accessors to basic data *)
     VertexCoordinates[mesh] := Transpose @ VertexCoordinatesTr[mesh],
     FaceList[mesh] := Transpose @ FaceListTr[mesh],
     EdgePairs[mesh] := Transpose @ EdgePairsTr[mesh],
     EdgePairs[mesh, patt_] := Cases[Transpose @ EdgePairsTr[mesh], patt, {1}],
     FaceList[mesh, patt_] := Cases[Transpose @ FaceListTr[mesh], patt, {1}],
     VertexList[mesh, patt_] := Cases[VertexList[mesh], patt, {1}],
     
     (* #PropertyList method extension *)
     PropertyList[mesh] := Union[
       VertexPropertyList[mesh],
       EdgePropertyList[mesh],
       FacePropertyList[mesh]],
     
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
       ReplaceAll[Map[Part[idx, #]&, is, {-2}], 0 -> $Failed]],

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
         ListQ[e] || ArrayQ[e], Map[Part[idx, #]&, EdgeIndex[mesh, e], {-2}],
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
     (* SumOverNeighbors *)
     SumOverNeighbors[mesh] := SumOverNeighborsMatrix[mesh],
     SumOverNeighbors[mesh, u_] := Dot[SumOverNeighborsMatrix[mesh], u],


     (* ================================ Geometrical Properties ================================= *)

     (* ------------------------------------------ Faces ---------------------------------------- *)

     (* #FaceAngles *)
     FaceAngleCosinesTr[mesh] :> CalculateFaceAngleCosinesTr[
       VertexIndex[mesh, FaceListTr[mesh]],
       VertexCoordinatesTr[mesh]],
     FaceAngleCosinesTr[mesh, Xtr_] := CalculateFaceAngleCosinesTr[
       VertexIndex[mesh, FaceListTr[mesh]],
       Xtr],
     FaceAnglesTr[mesh, Xtr_] := CalculateFaceAnglesTr[
       VertexIndex[mesh, FaceListTr[mesh]],
       Xtr],
     FaceAngleCosines[mesh, Xx_] := CalculateFaceAngleCosines[
       VertexIndex[mesh, FaceList[mesh]],
       Xx],
     FaceAngles[mesh, Xx_] := CalculateFaceAngles[VertexIndex[mesh, FaceList[mesh]], Xx],
     FaceAnglesTr[mesh]       := ArcCos @ FaceAngleCosinesTr[mesh],
     FaceAngleCosines[mesh]   := Transpose @ FaceAngleCosinesTr[mesh],
     FaceAngles[mesh]         := Transpose @ ArcCos @ FaceAngleCosinesTr[mesh],

     (* #FaceNormals *)
     FaceNormalsTr[mesh] :> CalculateFaceNormalsTr[
       VertexIndex[mesh, FaceListTr[mesh]],
       VertexCoordinatesTr[mesh]],
     FaceNormals[mesh] := Transpose @ FaceNormalsTr[mesh],
     FaceNormalsTr[mesh, Xtr_] := CalculateFaceNormalsTr[
       VertexIndex[mesh, FaceListTr[mesh]],
       Xtr],
     FaceNormals[mesh, X_] := Transpose @ CalculateFaceNormalsTr[
       VertexIndex[mesh, FaceListTr[mesh]],
       Transpose @ X],
     
     (* #FaceAxes *)
     FaceAxesTr[mesh] :> CalculateFaceAxes3DTr[
       VertexIndex[mesh, FaceListTr[mesh]],
       VertexCoordinatesTr[mesh],
       FaceNormalsTr[mesh]],
     FaceAxes[mesh] := Transpose @ FaceAxesTr[mesh],
     FaceAxesTr[mesh, Xt_] := CalculateFaceAxesTr[VertexIndex[mesh, FaceListTr[mesh]], Xt],
     FaceAxes[mesh, Xx_] := Transpose @ CalculateFaceAxesTr[
       VertexIndex[mesh, FaceListTr[mesh]],
       Transpose @ Xx],

     (* #FaceBisectors *)
     FaceBisectorsTr[mesh] :> CalculateFaceBisectorsTr[
       VertexIndex[mesh, FaceListTr[mesh]],
       VertexCoordinatesTr[mesh]],
     FaceBisectors[mesh] := Transpose @ FaceBisectorsTr[mesh],
     FaceBisectorsTr[mesh, X_] := CalculateFaceBisectorsTr[
       VertexIndex[mesh, FaceListTr[mesh]],
       X],
     FaceBisectors[mesh, X_] := Transpose @ FaceBisectorsTr[mesh, X],

     (* #FaceCoordinates *)
     FaceCoordinatesTr[mesh] :> With[
       {Xt = VertexCoordinatesTr[mesh],
        Ft = VertexIndex[mesh, FaceListTr[mesh]]},
       Transpose @ {Xt[[All, Ft[[1]]]], Xt[[All, Ft[[2]]]], Xt[[All, Ft[[3]]]]}],
     FaceCoordinates[mesh] := Transpose[FaceCoordinatesTr[mesh], {3,2,1}],

     (* #FacePlaneCoordinates *)
     FacePlaneCoordinatesTr[mesh] :> CalculateFacePlaneCoordinatesTr[
       VertexIndex[mesh, FaceListTr[mesh]],
       VerteListTr[mesh],
       FaceAxesTr[mesh]],
     FacePlaneCoordinatesTr[mesh, Xt_] := CalculateFacePlaneCoordinatesTr[
       VertexIndex[mesh, FaceListTr[mesh]],
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
       {EL = VertexIndex[mesh, EdgePairsTr[mesh]],
        Xt = VertexCoordinatesTr[mesh]},
       Sqrt @ Total[(Xt[[All, EL[[1]]]] - Xt[[All, EL[[2]]]])^2]],
     EdgeLengths[mesh, X_] := With[
       {EL = VertexIndex[mesh, EdgePairsTr[mesh]],
        Xt = Transpose @ X},
       Sqrt @ Total[(Xt[[All, EL[[1]]]] - Xt[[All, EL[[2]]]])^2]],
     EdgeLengthsTr[mesh, Xt_] := With[
       {EL = VertexIndex[mesh, EdgePairsTr[mesh]]},
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
       EdgeIndex[mesh, ReplaceAll[es, UndirectedEdge -> List]]],

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
        E = VertexIndex[mesh, EdgePairsTr[mesh]]},
       With[
         {neivecs = NeighborhoodVectors[mesh],
          neis = Part[
            SplitBy[
              SortBy[
                Transpose[{Join[E[[1]], E[[2]]], Join[E[[2]], E[[1]]]}],
                First],
              First],
            All, All, 2]},
         MapThread[
           NeighborhoodSort3DCompiled,
           {X, neivecs, neis, X[[#]] & /@ VertexIndex[mesh, neis]}]]],
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

     (* Labels *)
     LabelVertexList[mesh, args___] := CorticalLabelVertexList[mesh, args],
   
     (* #MeshRegion *)
     MeshRegion[mesh] :> MeshRegion[
       Transpose@VertexCoordinatesTr[mesh],
       Polygon@Transpose@FaceListTr[mesh]],
     
     (* #BoundaryMeshRegion *)
     BoundaryMeshRegion[mesh] :> BoundaryMeshRegion[
       Transpose@VertexCoordinatesTr[mesh],
       Polygon@Transpose@FaceListTr[mesh]],
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

(* We want to be able to pass options to MeshRegion and BoundaryMeshRegion *)
Unprotect[CorticalMesh3D];
CorticalMesh3D /: MeshRegion[m_CorticalMesh3D, opts__] := MeshRegion[MeshRegion[m], opts];
CorticalMesh3D /: BoundaryMeshRegion[m_CorticalMesh3D, opts__] := BoundaryMeshRegion[
  BoundaryMeshRegion[m],
  opts];
CorticalMesh3D /: Graph[m_CorticalMesh3D, opts__] := Graph[Graph[m], opts];
Protect[CorticalMesh3D];

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
   VertexColors -> Automatic,
   VertexRenderingFunction -> None,
   EdgeRenderingFunction -> None,
   FaceRenderingFunction -> Automatic}];
Protect[$CortexPlotOptions];

(* $CortexProjectionOptions is a list of the options to CorticalMap that affect the map projection;
 * i.e., these options require the vertex coordinates to be changed if set.
 *)
$CortexProjectionOptions = {Method, Center, Exclusions, Radius, AffineTransform};
Protect[$CortexProjectionOptions];

Options[CorticalMap] = Join[
  $CortexPlotOptions,
  {MetaInformation -> {},
   Properties -> None,
   Method -> "Orthographic",
   Center -> Automatic,
   Exclusions -> Automatic,
   Radius -> Full,
   AffineTransform -> None}];
DefineImmutable[
  CorticalMap[mesh_?CorticalMeshQ, optsPatt:OptionsPattern[]] :> map,
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
     SphericalMesh[map] -> With[
       {smesh = SourceMesh[map]},
       With[
         {smeta = Association[Options[smesh, MetaInformation] /. None -> {}]},
         With[
           {Xsph = Which[
              KeyExistsQ[smeta, "SphericalCoordinates"], smeta["SphericalCoordinates"],
              KeyExistsQ[smeta, "SphericalMesh"], VertexCoordinates@smeta["SphericalMesh"],
              KeyExistsQ[smeta, SphericalMesh], VertexCoordinates@smeta[SphericalMesh],
              True, None]},
           If[Xsph === None, smesh, CorticalMesh[smesh, VertexCoordinates -> Xsph]]]]],

     (* These two options are private to the namespace and are used so that the the projection 
        data in the map doesn't depend on options like MetaInformation. We don't want 
        MetaInformation to trigger a recalculation of the vertex coordinates when changed.     *)
     InitialOptions[map] -> Normal@Join[
       Association[$CortexPlotOptions],
       Association@FilterRules[Options@SourceMesh[map], Options[CorticalMap]],
       Association@Select[#, #[[2]] =!= Automatic&]&@Replace[
         "CorticalMap",
         Append[Options[mesh, MetaInformation], "CorticalMap" -> {}]],
       Association@Select[{optsPatt}, #[[1]] =!= Properties && #[[2]] =!= Automatic&]],
     ProjectionOptions[map] = Select[
       InitialOptions[map],
       MemberQ[$CortexProjectionOptions, #[[1]]]&],
     NonprojectionOptions[map] = Select[
       InitialOptions[map],
       !MemberQ[$CortexProjectionOptions, #[[1]]]&],
     
     (* Options is settable, but depends on nothing but the initial options
        (and should have nothing downstream but CorticalMapQ); 
        note that these are the graphics options only *)
     Options[map] -> Join[ProjectionOptions[map], NonprojectionOptions[map]],
     MetaInformation[map] := Fold[
       Replace,
       MetaInformation,
       {NonprojectionOptions[map], MetaInformation -> {}}],

     
     (* Here, we setup all of the special map-specific parameters and translate them *)
     TranslatedCenter[map] -> CorticalMapTranslateCenter[
       SphericalMesh[map],
       Center /. Join[ProjectionOptions[map], Options[CorticalMap]]],
     Inclusions[map] -> CorticalMapTranslateExclusions[
       SphericalMesh[map],
       Method /. Join[ProjectionOptions[map], Options[CorticalMap]],
       TranslatedCenter[map],
       Exclusions /. Join[ProjectionOptions[map], Options[CorticalMap]],
       Radius /. Join[ProjectionOptions[map], Options[CorticalMap]]],
     TransformationFunctions[map] -> With[
       {mthd = CorticalMapTranslateMethod[
          SphericalMesh[map],
          Method /. Join[ProjectionOptions[map], Options[CorticalMap]],
          TranslatedCenter[map],
          Inclusions[map],
          Radius /. Join[ProjectionOptions[map], Options[CorticalMap]]],
        tx = (AffineTransform /. ProjectionOptions[map]) /. AffineTransform|None -> Identity},
       With[
         {itx = If[tx === Identity, tx, InverseFunction[tx]]},
         {Transpose@tx@Transpose[mthd[[1]][##]]&, mthd[[2]][Transpose@itx@Transpose[#]]&}]],
   

     (* ====================================== Properties ======================================= *)
     
     (* Note that these Properties literally just steal the mesh's properties *)
     (* #VertexProperties [private] *)
     VertexProperties[map] = With[
       {idx = Inclusions[map][[1]]},
       Map[
         Function@With[
           {sym = TemporarySymbol["prop"<>#[[1]]], name = #[[1]], rule = #},
           name :> SetSafe[
             sym,
             With[{dat = Replace[name, rule]}, If[dat === $Failed, dat, dat[[idx]]]]]],
       VertexProperties[SourceMesh[map]]]],
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

     EdgeVertexProperties[map] :> With[
       {pp = Properties[map],
        p = Normal@VertexDataset[map]},
       Transpose[p[[#]]& /@ IndexedEdgePairsTr[map]]],
     FaceVertexProperties[map] :> With[
       {pp = Properties[map],
        p = Normal@VertexDataset[map]},
       Transpose[p[[#]]& /@ IndexedFaceListTr[map]]],
     
     (* Vertex, Edge, and Face Property associations and datasets *)
     VertexPropertyAssociation[map] :> With[
       {pp = Properties[map]},
       Association@Table[
         prop -> VertexPropertyValues[map, prop],
         {prop, VertexPropertyList[map]}]],
     VertexDataset[map] :> With[
       {cols = VertexPropertyList[map],
        n = VertexCount[map],
        pp = Properties[map]},
       Dataset@Map[
         Function@Association@Thread[cols -> #],
         Transpose@Replace[
           Map[Function@Normal@VertexPropertyValues[map, #], cols],
           Except[_List] :> ConstantArray[$Failed, n],
           {1}]]],
     EdgePropertyAssociation[map] :> With[
       {dep1 = Properties[map],
        dep2 = EdgeVertexProperties[map]},
       Association@Table[
         prop -> EdgePropertyValues[map, prop],
         {prop, EdgePropertyList[map]}]],
     EdgeDataset[map] :> With[
       {cols = EdgePropertyList[map],
        n = EdgeCount[map],
        dep1 = Properties[map],
        dep2 = EdgeVertexProperties[map]},
       Dataset@Map[
         Function@Association@Thread[cols -> #],
         Transpose@Replace[
           Map[Function@EdgePropertyValues[map, #], cols],
           Except[_List] :> ConstantArray[$Failed, n],
           {1}]]],
     FacePropertyAssociation[map] :> With[
       {dep1 = Properties[map],
        dep2 = FaceVertexProperties[map]},
       Association@Table[
         prop -> FacePropertyValues[map, prop],
         {prop, FacePropertyList[map]}]],
     FaceDataset[map] :> With[
       {cols = FacePropertyList[map],
        n = FaceCount[map],
        dep1 = Properties[map],
        dep2 = FaceVertexProperties[map]},
       Dataset@Map[
         Function@Association@Thread[cols -> #],
         Transpose@Replace[
           Map[Function@FacePropertyValues[map, #], cols],
           Except[_List] :> ConstantArray[$Failed, n],
           {1}]]],

     (* #PropertyList *)
     PropertyList[map] := Union[
       VertexPropertyList[map],
       EdgePropertyList[map],
       FacePropertyList[map]],
     (* Note that property value and set property are below, outside the immutable definition *)


     (* ======================================== Members ======================================== *)

     (* Now we have the actual vertex/face/edge data *)
     VertexList[map] -> Part[VertexList[SourceMesh[map]], Inclusions[map][[1]]],
     VertexCoordinatesTr[map] = With[
       {f = TransformationFunctions[map][[1]]},
       f@VertexCoordinatesTr[SphericalMesh[map]][[All, Inclusions[map][[1]]]]],
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
       Length@Intersection[VertexProperties[map][[All,1]], $VertexReadOnlyProperties] > 0, Message[
         CorticalMap::badarg,
         "Unable to set read-only vertex property"],
       Length@Intersection[EdgeProperties[map][[All,1]], $EdgeReadOnlyProperties] > 0, Message[
         CorticalMap::badarg,
         "Unable to set read-only edge property"],
       Length@Intersection[FaceProperties[map][[All,1]], $FaceReadOnlyProperties] > 0, Message[
         CorticalMap::badarg,
         "Unable to set read-only face property"],
       {} != Complement[
         Options[map][[All,1]],
         Cases[Options[CorticalMap][[All, 1]], Except[Properties], {1}]],
       Message[CorticalMap::badarg, "Unrecognized option given to CorticalMap"],
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
       SparseArray[Transpose[{Join[Et[[1]], Et[[2]]], Join[Et[[2]], Et[[1]]]}] -> Join[RE, RE]]],
     EdgeIndex[map, (List|UndirectedEdge)[a_Integer, b_Integer]] := With[
       {id = EdgeIndexArray[map][[a,b]]},
       If[id == 0, $Failed, id]],
     EdgeIndex[map, list:{(List|UndirectedEdge)[a_Integer, b_Integer]...}] := Extract[
       EdgeIndexArray[map],
       list],
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
       ReplaceAll[Map[Part[idx, #]&, is, {-2}], 0 -> $Failed]],
     

     (* ======================================== Delays ========================================= *)

     (* These indices are simple conversions of the vertex labels into vertex indices: *)
     IndexedEdgePairsTr[map] :> VertexIndex[map, EdgePairsTr[map]],
     IndexedEdgePairs[map] := Transpose@IndexedEdgePairsTr[map],
     IndexedEdgeList[map] := UndirectedEdge@@@IndexedEdgePairs[map],
     IndexedFaceListTr[map] :> VertexIndex[map, FaceListTr[map]],
     IndexedFaceList[map] :> Transpose@IndexedFaceListTr[map],

     (* Keep track of nearest vertices... *)
     Nearest[map] :> Nearest[VertexCoordinates[map] -> VertexList[map]],
     NearestFaceCenters[map] :> Nearest@Rule[
       Transpose@Total@Transpose@FaceCoordinatesTr[map] / 3.0,
       Automatic],

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
       {Ft = FaceListTr[map],
        RT = Range[FaceCount[map]]},
       With[
         {edges = {
            Join[
              Range[EdgeCount[map]],
              Extract[
                EdgeIndexArray[map],
                Transpose @ {
                  Join[Ft[[1]], Ft[[2]], Ft[[3]]],
                  Join[Ft[[2]], Ft[[3]], Ft[[1]]]}]],
            Join[ConstantArray[0, EdgeCount[map]], RT, RT, RT]}},
         With[
           {idx = SplitBy[
              Transpose @ edges[[All, Ordering[edges[[1]]]]],
              First]},
           idx[[All, 2;;All, 2]]]]],
     (* extensions of the opposite indices; edge/face lists *)
     VertexEdgeList[map, i_Integer] := Part[VertexEdgeList[map], VertexIndex[map, i]],
     VertexEdgeList[map, l_List] := With[
       {idx = VertexEdgeList[map]},
       Map[Part[idx, #]&, VertexIndex[map, l], {-2}]],
     VertexFaceList[map, i_Integer] := Part[VertexFaceList[map], VertexIndex[map, i]],
     VertexFaceList[map, l_List] := With[
       {idx = VertexFaceList[map]},
       Map[Part[idx, #]&, VertexIndex[map, i], {-2}]],
     EdgeFaceList[map, e_] := With[
       {idx = EdgeFaceList[map]},
       Which[
         Head[e] === UndirectedEdge, Part[idx, EdgeIndex[map, e]],
         MatchQ[e, {_Integer, _Integer}], Part[idx, EdgeIndex[map, e]],
         ListQ[e] || ArrayQ[e], Map[Part[idx, #]&, EdgeIndex[map, e], {-2}],
         True, $Failed]],

     (* #SumOverFacesMatrix, #SumOverEdgesMatrix, and #SumOverNeighbors *)
     SumOverFacesMatrix[map] :> SparseArray[
       Transpose[{Range[3*FaceCount[map]], Join @@ VertexIndex[map, FaceListTr[map]]}] -> 1],
     SumOverEdgesMatrix[map] :> SparseArray[
       Transpose[{Range[2*EdgeCount[map]], Join @@ VertexIndex[map, EdgePairsTr[map]]}] -> 1],
     SumOverNeighborsMatrix[map] :> SparseArray[
       Join@@MapThread[
         Thread[{#1, #2}]&,
         {Range@VertexCount[map], VertexIndex[map, NeighborhoodList[map]]}] -> 1,
       {VertexCount[map], VertexCount[map]},
       0],


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
     SumOverEdgesDirectedTr[map, datat_] := Dot[
       MapThread[Join, {datat, -datat}],
       SumOverEdgesMatrix[map]],
     SumOverEdgesDirected[map, data_] := Dot[
       Transpose @ Join[data, -data],
       SumOverEdgesMatrix[map]],
     SumOverEdgeVerticesTr[map, datat_] := Dot[Join @@ datat, SumOverEdgesMatrix[map]],
     SumOverEdgeVertices[map, data_] := Dot[Join @@ Transpose[data], SumOverEdgesMatrix[map]],
     (* SumOverNeighbors *)
     SumOverNeighbors[map] := SumOverNeighborsMatrix[map],
     SumOverNeighbors[map, u_] := Dot[SumOverNeighborsMatrix[map], u],


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
       VertexIndex[map, FaceListTr[map]],
       VertexCoordinatesTr[map]],
     FaceAngleCosinesTr[map, Xtr_] := CalculateFaceAngleCosinesTr[
       VertexIndex[map, FaceListTr[map]],
       Xtr],
     FaceAnglesTr[map, Xtr_] := CalculateFaceAnglesTr[VertexIndex[map, FaceListTr[map]], Xtr],
     FaceAngleCosines[map, Xx_] := CalculateFaceAngleCosines[VertexIndex[map, FaceList[map]], Xx],
     FaceAngles[map, Xx_] := CalculateFaceAngles[VertexIndex[map, FaceList[map]], Xx],
     FaceAnglesTr[map] := ArcCos @ FaceAngleCosinesTr[map],
     FaceAngleCosines[map] := Transpose @ FaceAngleCosinesTr[map],
     FaceAngles[map] := Transpose @ ArcCos @ FaceAngleCosinesTr[map],

     (* #FaceAxes *)
     FaceAxesTr[map] :> CalculateFaceAxes2DTr[
       VertexIndex[map, FaceListTr[map]],
       VertexCoordinatesTr[map]],
     FaceAxes[map] := Transpose @ FaceAxesTr[map],
     FaceAxesTr[map, Xt_] := CalculateFaceAxesTr[VertexIndex[map, FaceListTr[map]], Xt],
     FaceAxes[map, Xx_] := Transpose @ CalculateFaceAxesTr[
       VertexIndex[map, FaceListTr[map]],
       Transpose @ Xx],

     (* #FaceBisectors *)
     FaceBisectorsTr[map] :> CalculateFaceBisectorsTr[
       VertexIndex[map, FaceListTr[map]],
       VertexCoordinatesTr[map]],
     FaceBisectors[map] := Transpose @ FaceBisectorsTr[map],
     FaceBisectorsTr[map, X_] := CalculateFaceBisectorsTr[VertexIndex[map, FaceListTr[map]], X],
     FaceBisectors[map, X_] := Transpose @ FaceBisectorsTr[map, X],

     (* #FaceCoordinates *)
     FaceCoordinatesTr[map] :> With[
       {Xt = VertexCoordinatesTr[map],
        Ft = VertexIndex[map, FaceListTr[map]]},
       Transpose @ {Xt[[All, Ft[[1]]]], Xt[[All, Ft[[2]]]], Xt[[All, Ft[[3]]]]}],
     FaceCoordinates[map] := Transpose[FaceCoordinatesTr[map], {3,2,1}],

     (* #FacePlaneCoordinates *)
     FacePlaneCoordinatesTr[map] :> CalculateFacePlaneCoordinatesTr[
       VertexIndex[map, FaceListTr[map]],
       VerteListTr[map],
       FaceAxesTr[map]],
     FacePlaneCoordinatesTr[map, Xt_] := CalculateFacePlaneCoordinatesTr[
       VertexIndex[map, FaceListTr[map]],
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
       EdgeIndex[map, ReplaceAll[es, UndirectedEdge -> List]]],
     EdgeLengths[map, X_] := With[
       {EL = VertexIndex[map, EdgePairsTr[map]],
        Xt = Transpose @ X},
       Sqrt @ Total[(Xt[[All, EL[[1]]]] - Xt[[All, EL[[2]]]])^2]],
     EdgeLengthsTr[map, Xt_] := With[
       {EL = VertexIndex[map, EdgePairsTr[map]]},
       Sqrt @ Total[(Xt[[All, EL[[1]]]] - Xt[[All, EL[[2]]]])^2]],


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
         MapThread[
           NeighborhoodSort2DCompiled,
           {X, neis, X[[#]] & /@ VertexIndex[map, neis]}]]],
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


     (* ======================================= Functions ======================================= *)

     (* Extension of Options... *)
     Options[map, opt_] := Replace[opt, Options[map], If[ListQ[opt], {1}, {0}]],

     (* #Reproject *)
     ReprojectTr[map, Xtr_List] := With[
       {f = TransformationFunctions[map][[1]],
        m = SphericalMesh[map]},
       If[Dimensions[Xtr] == Dimensions@VertexCoordinatesTr[m],
         Clone[map, VertexCoordinatesTr -> f[Xtr[[All, Inclusions[map][[1]]]]]],
         (Message[Reproject::badarg, "Dimensions of coordinates do not match spherical mesh"];
          $Failed)]],
     Reproject[map, X_List] := ReprojectTr[map, Transpose[X]],
     Reproject[map, mesh_?CoticalMeshQ] := If[
       Dimensions[VertexCoordinateTr[mesh]] == Dimensions[VertexCoordinatesTr[SourceMesh[map]]],
       Clone[map, SourceMesh -> mesh],
       (Message[Reproject::badarg, "Dimensions of new mesh coordinates do not match source mesh"];
        $Failed)],

     (* #InverseProject *)
     InverseProjectTr[map] := With[
       {f = TransformationFunctions[map][[2]]},
       If[f === None,
         (Message[InverseProject::unin]; $Failed),
         f[VertexCoordinatesTr[map]]]],
     InverseProjectTr[map, XTr_List] := With[
       {f = TransformationFunctions[map][[2]]},
       If[f === None,
         (Message[InverseProject::unin]; $Failed),
         f[Xtr]]],
     InverseProjectTr[map, m_?CorticalMapQ] := With[
       {f = TransformationFunctions[map][[2]]},
       If[f === None,
         (Message[InverseProject::unin]; $Failed),
         f[VertexCoordinatesTr[m]]]],
     InverseProject[map] := Transpose @ InverseProjectTr[map],
     InverseProject[map, X_List] := Transpose @ InverseProjectTr[map, Transpose[X]],
     InverseProject[map, m_?CorticalMapQ] := Transpose @ InverseProjectTr[
       map,
       VertexCoordinatesTr[m]],

     (* #InverseProjectVectors *)
     InverseProjectVectorsTr[map, Ut_List] := With[
       {sourceMesh = SphericalMesh[map]},
       With[
         {Xtp = InverseProjectTr[map],
          XUtp = InverseProjectTr[map, VertexCoordinatesTr[map] + Ut],
          norms = ColumnNorms[Ut],
          normals = Part[
            VectorNormalsTr[sourceMesh],
            All,
            VertexIndex[sourceMesh, VertexList[map]]]},
         With[
           {Utp0 = (Utp - XUtp)},
           {norms, norms, norms} * NormalizeColumns[
              Utp0 - (normals * Table[#, {3}])& @ Total[Utp0 * normals]]]]],
     InverseProjectVectors[map, U_List] := Transpose @ InverseProjectVectors[map, Transpose[U]],

     (* ======================================= Interfaces ====================================== *)

     (* Labels *)
     LabelVertexList[map, args___] := CorticalLabelVertexList[map, args],

     (* #BoundaryMeshRegion *)
     BoundaryMeshRegion[map] :> BoundaryMeshRegion[
       Transpose@VertexCoordinatesTr[map],
       Polygon@Transpose@IndexedFaceListTr[map]],

     (* #MeshRegion *)
     MeshRegion[map] :> MeshRegion[
       Transpose@VertexCoordinatesTr[map],
       Polygon@Transpose@IndexedFaceListTr[map]],
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
  {SourceMesh, SphericalMesh, VertexCoordinates},
  Cases[Options[CorticalMap][[All,1]], Except[Properties], {1}]];
Protect[$CorticalMapReconstructionOptions];
CorticalMap[map_?CorticalMapQ, args___Rule] := Check[
  With[
    {optsarg = {args} // Function[
       If[Complement[#[[All,1]], $CorticalMapReconstructionOptions] == {},
         #,
         Message[
           CorticalMap::badarg,
           "CorticalMap[] reconstructor may only be passed CorticalMap options or SourceMesh," <>
           " SphericalMesh, or VertexCoordinates"]]]},
    If[Length[optsarg] == 0,
      map,
      With[
        {source = Replace[SourceMesh, optsarg],
         sphere = Replace[SphericalMesh, optsarg],
         coords = Replace[VertexCoordinates, optsarg],
         sourceMesh = SourceMesh[map],
         opts = Select[
           optsarg,
           (#[[1]] =!= SourceMesh && #[[1]] =!= VertexCoordinates && #[[1]] =!= SphericalMesh)&]},
        With[
          {ss = Which[
             source === SourceMesh && sphere === SphericalMesh, None,
             source =!= SourceMesh && sphere === SphericalMesh, source,
             source === SourceMesh && sphere =!= SphericalMesh, With[
               {meta = Replace[Options[sourceMesh, MetaInformation], None -> {}]},
               CorticalMesh[
                 sourceMesh,
                 MetaInformation -> Normal@Append[
                   Association[meta],
                   "SphericalMesh" -> sphere]]],
             source =!= SourceMesh && sphere =!= SphericalMesh, With[
               {meta = Replace[Options[source, MetaInformation], None -> {}]},
               CorticalMesh[
                 source,
                 MetaInformation -> Normal@Append[
                   Association[meta],
                   "SphericalMesh" -> sphere]]]]},
          Clone[
            map,
            Sequence @@ With[
              {newOpts = Flatten[
                 {If[ss =!= None, SourceMesh -> source, {}],
                  If[coords =!= VertexCoordinates, VertexCoordinatesTr -> Transpose[coords], {}],
                  If[Length[opts] > 0,
                    With[
                      {opts = Fold[
                         Function @ With[
                           {edit = Replace[#1, (Rule|RuleDelayed)[#2[[1]], _] :> #2, {1}]},
                           If[SameQ[edit, #1], Append[#1, #2], edit]],
                         Options[map],
                         opts]},
                      {ProjectionOptions -> Select[
                         opts,
                         MemberQ[$CortexProjectionOptions, #[[1]]]&],
                       NonprojectionOptions -> Select[
                         opts,
                         !MemberQ[$CortexProjectionOptions, #[[1]]]&]}],
                    {}]}],
               nowOpts = Association@Options[map]},
              Map[
                If[0 == Length@Select[#[[2]], !SameQ[#[[2]], nowOpts[#[[1]]]]&], Nothing, #]&,
                newOpts]]]]]]],
  $Failed];

Unprotect[CorticalMesh2D];
CorticalMesh2D /: MeshRegion[m_CorticalMesh2D, opts__] := MeshRegion[MeshRegion[m], opts];
CorticalMesh2D /: BoundaryMeshRegion[m_CorticalMesh2D, opts__] := BoundaryMeshRegion[
  BoundaryMeshRegion[m],
  opts];
CorticalMesh2D /: Graph[m_CorticalMesh2D, opts__] := Graph[Graph[m], opts];
Protect[CorticalMesh2D];

(* #CorticalMapQ *)
CorticalMapQ[_] := False;

(* #CorticalObjectQ *)
CorticalObjectQ[c_] := Or[CorticalMeshQ[c], CorticalMapQ[c]];

(* We want to make sure to have a nicely formatted output for a mesh or projection *)

MakeBoxes[mesh_CorticalMesh3D, form_] := MakeBoxes[#]&[
  With[
    {style = {
       FontSize -> 11,
       FontColor -> Gray,
       FontFamily -> "Arial",
       FontWeight -> "Thin"}},
    Row[
      {"CorticalMesh"[
         Panel[
           Grid[
             MapThread[
               Function[
                 {Spacer[4], Style[#1, Sequence @@ style],
                  Spacer[2], #2, Spacer[4]}],
               {{"Vertex Count:", "Edge Count:", "Face Count:"},
                {VertexCount[mesh], EdgeCount[mesh], FaceCount[mesh]}}],
             Alignment -> Table[{Right, Right, Center, Left, Left}, {3}]]]]},
      BaseStyle -> Darker[Gray]]]];
MakeBoxes[mesh_CorticalMesh2D, form_] := MakeBoxes[#]&[
  With[
    {style = {
       FontSize -> 11,
       FontColor -> Gray,
       FontFamily -> "Arial",
       FontWeight -> "Thin"}},
    Row[
      {"CorticalMap"[
         Panel[
           Grid[
             MapThread[
               Function[
                 {Spacer[4], Style[#1, Sequence @@ style],
                  Spacer[2], #2, Spacer[4]}],
               {{"Vertex Count:", "Edge Count:", "Face Count:"},
                {VertexCount[mesh], EdgeCount[mesh], FaceCount[mesh]}}],
             Alignment -> 
             Table[{Right, Right, Center, Left, Left}, {3}]]]]},
      BaseStyle -> Darker[Gray]]]];

(* Protect these functions... *)
Protect[CorticalMesh, CorticalMeshQ, Inclusions, VertexCoordinates, VertexCoordinatesTr, EdgePairs, 
        EdgePairsTr, FaceIndexArray, EdgeIndexArray, FaceIndex, EdgeIndex, EdgeCoordinates, 
        EdgeCoordinatesTr, OptionalProperties, VertexProperties, EdgeProperties, FaceProperties,
        FaceAxes, FaceAxesTr, FaceAngleCosines, FaceAngleCosinesTr, FaceAngles, FaceAnglesTr, 
        FaceCoordinates, FaceCoordinatesTr, FacePlaneCoordinates, FacePlaneCoordinatesTr, 
        FaceNormals, FaceNormalsTr, EdgeLengths, NeighborhoodList, NeighborhoodAngles, 
        NeighborhoodBisectors, NeighborhoodEdgeLengths, SourceImage, VertexEdgeList, VertexFaceList,
        EdgeFaceList, IndexedEdgeList, IndexedFaceList, IndexedEdgeListTr, IndexedFaceListTr,
        IndexedEdgePairs, IndexedEdgePairsTr, VertexPropertyAssociation, EdgePropertyAssociation,
        FacePropertyAssociation, VertexDataset, EdgeDataset, FaceDataset];

(* #NearestFace ***********************************************************************************)
Options[NearestFace] = {NearestPoints -> False};
NearestFace[mesh_?CorticalObjectQ, x0_, OptionsPattern[]] := With[
  {doPoints = Replace[OptionValue[NearestPoints], Except[False] -> True],
   region = MeshRegion[mesh],
   dims = If[CorticalMapQ[mesh], 2, 3]},
  With[
    {x = RegionNearest[
       region,
       Which[
         VectorQ[x0, NumericQ] && Length[x0] == dims, {x0},
         !MatrixQ[x0, NumericQ], (
           Message[
             If[CorticalMeshQ[mesh], CorticalMesh::error, CorticalMap::error],
             NearestFace,
             "argument must be a numberic vector or matrix"];
           x0),
         Length@First[x0] == dims, x0,
         Length[x0] == dims, Transpose[x0], 
         True, (
           Message[
             CorticalMesh::error,
             NearestFace,
             "argument must have "<> ToString[dims] <>" dimensions"];
           x0)]]},
    With[
      {idcs = If[VectorQ[x0], First[#], #]&@Last@Transpose@Region`Mesh`MeshMemberCellIndex[
         region,
         x]},
      If[doPoints, {idcs, If[VectorQ[x0], x[[1]], x]}, idcs]]]];
Protect[NearestFace];

(* #CortexAddress *********************************************************************************)
CortexAddress[mesh_?CorticalObjectQ, X0_ /; MatrixQ[X0, NumericQ]] := With[
  {idcs = NearestFace[
     mesh,
     Which[
       Length@First[X0] == Length@VertexCoordinatesTr[mesh], X0,
       Length[X0] == Length@VertexCoordinatesTr[mesh], Transpose[X0],
       True, Message[CortexAddress::err, "Dimensionality of argument is incorrect"]],
     NearestPoints -> True]},
  With[
    {faces = FaceListTr[mesh][[All, idcs[[1]]]],
     coords = Transpose[FaceCoordinatesTr[mesh][[All, All, idcs[[1]]]]],
     triArea = Function[0.5 * (#1[[2]]*(#3[[1]] - #2[[1]]) + #1[[1]]*(#2[[2]] - #3[[2]])
                               + #2[[1]]*#3[[2]] - #2[[2]]*#3[[1]])],
     Xt = Transpose[idcs[[2]]]},
    With[
      {area = triArea@@coords,
       a1 = triArea@@ReplacePart[coords, 1 -> Xt],
       a2 = triArea@@ReplacePart[coords, 2 -> Xt],
       a3 = triArea@@ReplacePart[coords, 3 -> Xt]},
      With[
        {aZero = 1 - Unitize@Chop[area]},
        With[
          {bary = ConstantArray[aZero/3.0, 2] + (((1 - aZero)*# / (area + aZero))&/@{a1, a2})},
          If[Length[X0] == Length[Xt],
            {faces, bary},
            Transpose[{Transpose[faces], Transpose[bary]}]]]]]]];
CortexAddress[mesh_?CorticalMapQ, X0:{_?NumericQ, _?NumericQ}] := First@CortexAddress[mesh, {X0}];
CortexAddress[mesh_?CorticalMeshQ, X0:{_?NumericQ, _?NumericQ, _?NumericQ}] := First@CortexAddress[
  mesh, {X0}];
CortexAddress[mesh_?CorticalObjectQ] := Function@CortexAddress[mesh, #];
Protect[CortexAddress];

(* #CortexLookup **********************************************************************************)
CortexLookup[mesh_?CorticalObjectQ, {faces0_ /; MatrixQ[faces0, IntegerQ],
                                     bary0_  /; MatrixQ[bary0, NumericQ]  }] := With[
  {faces = VertexIndex[mesh, If[Length@First[faces0] == 3, Transpose[faces0], faces0]],
   bary = If[Length@First[bary0] == 2, Transpose[bary0], bary0],
   Xt = VertexCoordinatesTr[mesh]},
  With[
    {coords = (Xt[[All, #]])& /@ faces,
     b3 = 1 - Total[bary]},
    Total@MapThread[
      #1 * ConstantArray[#2, Length[#1]]&,
      {coords, Append[bary, b3]}]]];
CortexLookup[mesh_, a:{{{_Integer,_Integer,_Integer}, {_?NumericQ, _?NumericQ}}..}] := Transpose[
  CortexLookup[mesh, Transpose[a]]];
CortexLookup[mesh_, addr:{{_Integer, _Integer, _Integer}, {_?NumericQ, _?NumericQ}}] := First[
  CortexLookup[mesh, {addr}]];
CortexLookup[mesh_] := Function@CortexLookup[mesh, #];
Protect[CortexLookup];

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
       {VertexCoordinates|"VertexCoordinates" :> {
          VertexList -> VertexCoordinates[mesh],
          EdgeList -> EdgeCoordinates[mesh],
          FaceList -> FaceCoordinates[mesh]},
        VertexList|"VertexList" :> {
          VertexList -> $Failed,
          EdgeList -> EdgePairs[mesh],
          FaceList -> FaceList[mesh]},
        Label|"Label" :> {
          VertexList -> VertexList[mesh],
          EdgeList -> $Failed,
          FaceList -> $Failed},
        VertexIndex|"VertexIndex" :> {
          VertexList -> Range@VertexCount[mesh],
          EdgeList -> IndexedEdgePairs[mesh],
          FaceList -> IndexedFaceList[mesh]},
        VertexNormals|"VertexNormals" :> {
          VertexList -> VertexNormals[mesh],
          EdgeList -> $Failed,
          FaceList -> $Failed},
        EdgeWeight|"EdgeWeight"|EdgeLengths|"EdgeLengths" :> {
          VertexList -> $Failed,
          EdgeList -> EdgeLengths[mesh],
          FaceList -> $Failed},
        FaceNormals|"FaceNormals" :> {
          VertexList -> $Failed,
          EdgeList -> $Failed,
          FaceList -> FaceNormals[mesh]},
        VertexProperties|"VertexProperties" :> {
          VertexList -> $Failed,
          EdgeList -> EdgeVertexProperties[mesh],
          FaceList -> FaceVertexProperties[mesh]}}]]},
  If[list === prop || list === $Failed,
    $Failed,
    list]];
PropertyValue[{mesh_?CorticalObjectQ, VertexList}, prop:Except[_List]] := With[
  {list = Replace[prop, Properties[mesh]]},
  If[list === prop || list === $Failed,
    Switch[prop,
      VertexCoordinates|"VertexCoordinates", VertexCoordinates[mesh],
      VertexNormals|"VertexNormals", If[CorticalMeshQ[mesh], VertexNormals[mesh], $Failed],
      Label|"Label", VertexList[mesh],
      VertexIndex|"VertexIndex", Range@VertexCount[mesh],
      _, $Failed],
    list[[1,2]]]];
PropertyValue[{mesh_?CorticalObjectQ, EdgeList}, prop:Except[_List]] := With[
  {list = Replace[prop, Properties[mesh]]},
  If[list === prop || list === $Failed,
    Switch[prop,
      EdgeWeight|"EdgeWeight"|EdgeLengths|"EdgeLengths", EdgeLengths[mesh],
      VertexCoordinates|"VertexCoordinates", EdgeCoordinates[mesh],
      VertexList|"VertexList", EdgePairs[mesh],
      VertexIndex|"VertexIndex", IndexedEdgePairs[mesh],
      VertexProperties|"VertexProperties", EdgeVertexProperties[mesh],
      _, $Failed],
    list[[2,2]]]];
PropertyValue[{mesh_?CorticalObjectQ, FaceList}, prop:Except[_List]] := With[
  {list = Replace[prop, Properties[mesh]]},
  If[list === prop || list === $Failed,
    Switch[prop,
      FaceNormals|"FaceNormals", If[CorticalMeshQ[mesh], FaceNormals[mesh], $Failed],
      VertexCoordinates|"VertexCoordinates", FaceCoordinates[mesh],
      VertexList|"VertexList", FaceList[mesh],
      VertexIndex|"VertexIndex", IndexedFaceList[mesh],
      VertexProperties|"VertexProperties", FaceVertexProperties[mesh],
      _, $Failed],
    list[[3,2]]]];
PropertyValue[{mesh_?CorticalObjectQ, i_Integer}, prop:Except[_List]] := With[
  {list = PropertyValue[{mesh, VertexList}, prop]},
  If[list === $Failed || list[[1]] === $Failed,
    $Failed,
    list[[1, 2]][[i]]]];
PropertyValue[{mesh_?CorticalObjectQ, e:(List|UndirectedEdge)[_Integer, _Integer]}, 
              prop:Except[_List]] := With[
  {list = PropertyValue[{mesh, EdgeList}, prop]},
  If[list === prop || list === $Failed || list[[2]] === $Failed,
    $Failed,
    list[[2, 2]][[EdgeIndex[mesh, e]]]]];
PropertyValue[{mesh_?CorticalObjectQ, f:{_Integer, _Integer, _Integer}}, prop:Except[_List]] := With[
  {list = PropertyValue[{mesh, FaceList}, prop]},
  If[list === prop || list === $Failed || list[[3]] === $Failed,
    $Failed,
    list[[3, 2]][[FaceIndex[mesh, f]]]]];
PropertyValue[mesh_?CorticalObjectQ, prop_List] := Map[PropertyValue[mesh, #]&, prop];
PropertyValue[{mesh_?CorticalObjectQ, x_}, prop_List] := Map[PropertyValue[{mesh, x}, #]&, prop];

(* PropertyList stuff *)
PropertyList[{mesh_?CorticalMeshQ, type:(VertexList|EdgeList|FaceList)}] := With[
  {props = Properties[mesh],
   idx = Replace[type, {VertexList -> 1, EdgeList -> 2, FaceList ->3}]},
  Join[
    {If[type === VertexList, "Label", "VertexList"], "VertexIndex", "VertexCoordinates"},
    Select[props, (Hold @@ #[[2,idx]])[[{2}]] =!= Hold[$Failed]&][[All, 1]],
    Which[
      type === VertexList, {"VertexNormals"},
      type === EdgeList,   {"EdgeWeight", "VertexProperties"},
      True,                {"FaceNormals", "VertexProperties"}]]];
PropertyList[{mesh_?CorticalMapQ, type:(VertexList|EdgeList|FaceList)}] := With[
  {props = Properties[mesh],
   idx = Replace[type, {VertexList -> 1, EdgeList -> 2, FaceList ->3}]},
  Join[
    {If[type === VertexList, "Label", "VertexList"], "VertexIndex", "VertexCoordinates"},
    Select[props, (Hold @@ #[[2,idx]])[[{2}]] =!= Hold[$Failed]&][[All, 1]],
    Which[
      type === VertexList, {},
      type === EdgeList,   {"EdgeWeight", "VertexProperties"},
      True,                {"VertexProperties"}]]];
PropertyList[{mesh_?CorticalObjectQ, i_Integer}] := With[
  {props = VertexPropertyList[mesh]},
  Select[props, PropertyValue[{mesh, i}, #] =!= $Failed&]];
PropertyList[{mesh_?CorticalObjectQ, e:(List|UndirectedEdge)[_Integer, _Integer]}] := With[
  {props = EdgePropertyList[mesh]},
  Select[props, PropertyValue[{mesh, e}, #] =!= $Failed&]];
PropertyList[{mesh_?CorticalObjectQ, f:{_Integer, _Integer, _Integer}}] := With[
  {props = FacePropertyList[mesh]},
  Select[props, PropertyValue[{mesh, f}, #] =!= $Failed&]];

(* SetProperty stuff... *)

SetProperty[{mesh_?CorticalObjectQ, type:(VertexList|EdgeList|FaceList)},
            prop_ -> vals_List] := With[
  {allList = Switch[type, 
     VertexList, VertexProperties[mesh],
     EdgeList, EdgeProperties[mesh],
     FaceList, FaceProperties[mesh]],
   propType = Switch[type, 
     VertexList, VertexProperties, EdgeList, EdgeProperties, FaceList, FaceProperties]},
  With[
    {list = Replace[prop, allList]},
    Which[
      Length[vals] != Length[type[mesh]], (
        Message[CorticalMesh::error, "SetProperty", "Incorrect value list length"];
        $Failed),
      (* If the property doesn't yet exist, add it. *)
      list === prop || list === $Failed, Switch[prop,
        VertexCoordinates|"VertexCoordinates", If[type === VertexList,
          Clone[mesh, VertexCoordinatesTr -> Transpose[vals]],
          (Message[CorticalMesh::error, "SetProperty", "VertexCoordinates not set with VertexList"];
           $Failed)],
        _, Clone[mesh, propType -> Append[allList, prop -> Normal[vals]]]],
      (* If the property already exists, we overwrite it *)
      True, Clone[
        mesh,
        propType -> Replace[allList, (Rule|RuleDelayed)[prop, _] -> (prop -> vals), {1}]]]]];
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
      {res = Normal[vals]},
      If[!ArrayQ[res] || Length[res] != Length[type[mesh]],
        Table[$Failed, {Length[type[mesh]]}],
        (sym = res)]];
    If[list === prop || list === $Failed, 
      (* If the property doesn't yet exist, add it. *)
      Switch[prop,
        VertexCoordinates|"VertexCoordinates", If[type === VertexList,
          Clone[mesh, VertexCoordinatesTr -> Transpose[vals]],
          (Message[CorticalMesh::error, "SetProperty", "VertexCoordinates not set with VertexList"];
           $Failed)],
        _, Clone[mesh, propType -> Append[allList, (prop :> sym)]]],
      (* If the property already exists, we overwrite it *)
      Clone[
        mesh,
        propType -> Replace[allList, (Rule|RuleDelayed)[prop, _] -> (prop :> sym), {1}]]]]];

(* bulk property setting *)
SetProperty[{mesh_?CorticalObjectQ,
             type:(VertexList|EdgeList|FaceList)}, rs:{(_Rule|_RuleDelayed)..}] := Fold[
  If[#1 === $Failed, $Failed, SetProperty[{#1, type}, #2]]&,
  mesh,
  rs];
(* confused property setting (note that :> is not valid here) *)
SetProperty[mesh_?CorticalObjectQ, rs:{_Rule..}] := Fold[SetProperty, mesh, rs];
SetProperty[mesh_?CorticalObjectQ, r:(_ -> vals_List)] := Switch[
  Length[vals],
  VertexCount[mesh], SetProperty[{mesh, VertexList}, r],
  EdgeCount[mesh],   SetProperty[{mesh, EdgeList}, r],
  FaceCount[mesh],   SetProperty[{mesh, FaceList}, r],
  _, (
    Message[
      CorticalMesh::error,
      "SetProperty",
      "Could not deduce type for property of length " <> ToString@Length[vals]];
    $Failed)];

(* Set individual vertices/edges/faces *)
SetProperty[{mesh_?CorticalObjectQ, obj:(_Integer|_List|_UndirectedEdge)}, prop_ -> val_] := With[
  {type = Switch[obj,
     _Integer,                                  VertexProperties -> obj,
     (List|UndirectedEdge)[_Integer, _Integer], EdgeProperties -> EdgeIndex[mesh, obj],
     {_Integer, _Integer, _Integer},            FaceProperties -> FaceIndex[mesh, obj],
     _, $Failed]},
  If[type === $Failed,
    (Message[CorticalMesh::error, "SetProperty", "Invalid mesh part (" <> obj <> ")"];
     $Failed),
    Which[
      prop === VertexCoordinates || (StringQ[prop] && prop === "VertexCoordinates"), If[
        type[[1]] =!= VertexProperties,
        (Message[
           CorticalMesh::error,
           "SetProperty",
           "VertexCoordinates can only be set for a vertex"];
         $Failed),
        Clone[mesh, VertexCoordinates -> ReplacePart[VertexCoordinates[mesh], obj -> val]]],
      prop === EdgeWeight || (StringQ[prop] && prop == "EdgeWeight"), (
        Message[CorticalMesh::error, "SetProperty", "EdgeWeight is a read-only property"];
        $Failed),
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
    (Message[CorticalMesh::error, "SetProperty", "Invalid mesh part (" <> obj <> ")"];
     $Failed),
    Which[
      prop === VertexCoordinates || (StringQ[prop] && prop == "VertexCoordinates"), If[
        type[[1]] =!= VertexProperties, 
        (Message[
           CorticalMesh::error,
           "SetProperty",
           "VertexCoordinates can only be set for a vertex"];
         $Failed),
        Clone[mesh, VertexCoordinates -> ReplacePart[VertexCoordinates[mesh], obj :> val]]],
      prop === EdgeWeight || (StringQ[prop] && prop == "EdgeWeight"), (
        Message[CorticalMesh::error, "SetProperty", "EdgeWeight is a read-only property"];
        $Failed),
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
  {list = PropertyValue[mesh, prop]},
  If[list === $Failed, 
    mesh,
    Clone[
      mesh,
      VertexProperties -> DeleteCases[VertexProperties[mesh], (Rule|RuleDelayed)[prop, _]],
      EdgeProperties -> DeleteCases[EdgeProperties[mesh], (Rule|RuleDelayed)[prop, _]],
      FaceProperties -> DeleteCases[FaceProperties[mesh], (Rule|RuleDelayed)[prop, _]]]]];
RemoveProperty[mesh_?CorticalObjectQ, prop:{_Rule..}] := Fold[RemoveProperty, mesh, prop];
RemoveProperty[{mesh_?CorticalObjectQ, t:(VertexList|EdgeList|FaceList)},
               prop:Except[_List]] := With[
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
RemoveProperty[{mesh_?CorticalObjectQ,obj:Except[VertexList|EdgeList|FaceList]}, prop_List] := Fold[
  RemoveProperty[{#1, obj}, #2]&,
  mesh,
  prop];
RemoveProperty[{mesh_?CorticalObjectQ, t:(VertexList|EdgeList|FaceList)}] := Clone[
  mesh,
  Rule[
    Switch[t,
      VertexList, VertexProperties,
      EdgeList,   EdgeProperties,
      FaceList,   FaceProperties],
    {}]];
RemoveProperty[mesh_?CorticalObjectQ] := Clone[
  mesh,
  VertexProperties -> {},
  EdgeProperties -> {},
  FaceProperties -> {}];

(* We also setup a convenience syntax for properties: /. to add a property *)
Unprotect[CorticalMesh2D, CorticalMesh3D];
(* Add a Vertex property *)
CorticalMesh2D /: ReplaceAll[m_CorticalMesh2D, repl_] := SetProperty[{m, VertexList}, repl];
CorticalMesh3D /: ReplaceAll[m_CorticalMesh3D, repl_] := SetProperty[{m, VertexList}, repl];
(* Add a Face property *)
CorticalMesh2D /: ReplaceRepeated[m_CorticalMesh2D, repl_] := SetProperty[{m, FaceList}, repl];
CorticalMesh3D /: ReplaceRepeated[m_CorticalMesh3D, repl_] := SetProperty[{m, FaceList}, repl];
(* Get a Vertex property *)
CorticalMesh2D /: ReplaceAll[repl_, m_CorticalMesh2D] := PropertyValue[{m, VertexList}, repl];
CorticalMesh3D /: ReplaceAll[repl_, m_CorticalMesh3D] := PropertyValue[{m, VertexList}, repl];
(* Get an Face property *)
CorticalMesh2D /: ReplaceRepeated[repl_, m_CorticalMesh2D] := PropertyValue[{m, FaceList}, repl];
CorticalMesh3D /: ReplaceRepeated[repl_, m_CorticalMesh3D] := PropertyValue[{m, FaceList}, repl];
Protect[CorticalMesh2D, CorticalMesh3D];

(* Individualized property functions *)
VertexPropertyList[mesh_?CorticalObjectQ] := PropertyList[{mesh, VertexList}];
EdgePropertyList[mesh_?CorticalObjectQ] := PropertyList[{mesh, EdgeList}];
FacePropertyList[mesh_?CorticalObjectQ] := PropertyList[{mesh, FaceList}];
VertexPropertyValues[mesh_?CorticalObjectQ, prop_] := PropertyValue[{mesh, VertexList}, prop];
EdgePropertyValues[mesh_?CorticalObjectQ, prop_] := PropertyValue[{mesh, EdgeList}, prop];
FacePropertyValues[mesh_?CorticalObjectQ, prop_] := PropertyValue[{mesh, FaceList}, prop];
SetVertexProperties[mesh_?CorticalObjectQ, prop_] := SetProperty[{mesh, VertexList}, prop];
SetEdgeProperties[mesh_?CorticalObjectQ, prop_] := SetProperty[{mesh, EdgeList}, prop];
SetFaceProperties[mesh_?CorticalObjectQ, prop_] := SetProperty[{mesh, FaceList}, prop];
RemoveVertexProperty[mesh_?CorticalObjectQ, prop_] := RemoveProperty[{mesh, VertexList}, prop];
RemoveEdgeProperty[mesh_?CorticalObjectQ, prop_] := RemoveProperty[{mesh, EdgeList}, prop];
RemoveFaceProperty[mesh_?CorticalObjectQ, prop_] := RemoveProperty[{mesh, FaceList}, prop];

(* #MapVertices ***********************************************************************************)
MapVertices[f_, mesh_?CorticalObjectQ] := Map[f, Normal@VertexDataset[mesh]];

(* #MapEdges **************************************************************************************)
MapEdges[f_, mesh_?CorticalObjectQ] := Map[f, Normal@EdgeDataset[mesh]];

(* #MapFaces **************************************************************************************)
MapFaces[f_, mesh_?CorticalObjectQ] := Map[f, Normal@FaceDataset[mesh]];

(* #SelectVertices ********************************************************************************)
SelectVertices[mesh_?CorticalObjectQ, f_] := Pick[
  VertexList[mesh],
  f /@ Normal@VertexDataset[mesh],
  True];

(* #SelectEdges ***********************************************************************************)
SelectEdges[mesh_?CorticalObjectQ, f_] := Pick[
  EdgeList[mesh],
  f /@ Normal@EdgeDataset[mesh],
  True];

(* #SelectFaces ***********************************************************************************)
SelectFaces[mesh_?CorticalObjectQ, f_] := Pick[
  FaceList[mesh],
  f /@ Normal@FaceDataset[mesh],
  True];

(* #SelectVertexIndices ***************************************************************************)
SelectVertexIndices[mesh_?CorticalObjectQ, f_] := Pick[
  Range@VertexCount[mesh],
  f /@ Normal@VertexDataset[mesh],
  True];

(* #SelectEdgeIndices *****************************************************************************)
SelectEdgeIndices[mesh_?CorticalObjectQ, f_] := Pick[
  Range@EdgeCount[mesh],
  f /@ Normal@EdgeDataset[mesh],
  True];

(* #SelectFaceIndices *****************************************************************************)
SelectFaceIndices[mesh_?CorticalObjectQ, f_] := Pick[
  Range@FaceCount[mesh],
  f /@ Normal@FaceDataset[mesh],
  True];

(* #SelectIndexedEdges ****************************************************************************)
SelectIndexedEdges[mesh_?CorticalObjectQ, f_] := Pick[
  IndexedEdgeList[mesh],
  f /@ Normal@EdgeDataset[mesh],
  True];

(* #SelectIndexedFaces ****************************************************************************)
SelectIndexedFaces[mesh_?CorticalObjectQ, f_] := Pick[
  IndexedFaceList[mesh],
  f /@ Normal@FaceDataset[mesh],
  True];

Protect[PropertyValue, SetProperty, RemoveProperty, PropertyList,
        VertexPropertyList, EdgePropertyList, FacePropertyList, VertexPropertyValues,
        EdgePropertyValues, FacePropertyValues, SetVertexProperties, SetEdgePropertues,
        SetFaceProperties, RemoveVertexProperty, RemoveEdgeProperty, RemoveFaceProperty,
        MapVertices, MapEdges, MapFaces, SelectVertices, SelectEdges, SelectFaces,
        SelectVertexIndices, SelectEdgeIndices, SelectFaceIndices, SelectIndexedEdges,
        SelectIndexedFaces];


(* #CorticalCurvatureColor ************************************************************************)
CorticalCurvatureColor[c_?NumericQ] := If[c > 0, GrayLevel[0.2], Gray];
CorticalCurvatureColor[c_ /; ArrayQ[c, _, NumericQ]] := Map[CorticalCurvatureColor, c];
CorticalCurvatureColor[a_?AssociationQ] := CorticalCurvatureColor @ Which[
  KeyExistsQ[a, "Curvature"], a["Curvature"],
  KeyExistsQ[a, Curvature], a[Curvature],
  True, -0.5];
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

(* Simple vertex color instructions interpreter *)
InterpretVertexColor[r_, vtxData_?AssociationQ] := Which[
  (* A string like "PolarAngle" or "Curvature" *)
  StringQ[r] && Head[CorticalColorData[r]] === CorticalColorSchema, CorticalColorData[r][vtxData],
  (* A Blended color, like {{"Curvature", "PolarAngle"}, 0.6} *)
  ListQ[r] && Length[r] == 2 && VectorQ[r[[1]], StringQ] && NumericQ[r[[2]]], Blend[
    Table[CorticalColorData[t][vtxData], {t, r[[1]]}],
    r[[2]]],
  (* Anything else (like a color) *)
  True, r];

(* A Function for figuring out the vertex colors *)
GetVertexColors[mesh_, vcolorsOpt_, colorFnOpt_] := With[
  {allGray = ConstantArray[Gray, VertexCount[mesh]]},
  If[ListQ[vcolorsOpt] && VertexCount[mesh] == Length[vcolorsOpt],
    vcolorsOpt,
    Replace[
      colorFnOpt,
      {Automatic :> With[
         {sel = SelectFirst[
            $CortexPlotDefaultColorSchemas,
            VectorQ[VertexPropertyValues[mesh, #[[1]]], NumericQ]&,
            None]},
         If[sel === None, 
           allGray,
           (sel[[2]])[mesh]]],
       None -> allGray,
       name_ /; Head[CorticalColorData[name]] === CorticalColorSchema :> With[
         {f = CorticalColorData[name]},
         f /@ Normal@VertexDataset[mesh]],
       f_ :> With[
         {known = Replace[f, $CortexPlotDefaultColorSchemas]},
         If[f =!= known,
           known[mesh],
           Map[
             InterpretVertexColor[f[#], #]&,
             Normal@VertexDataset[mesh]]]]}]]];

Options[CortexPlot3D] = Map[(#[[1]] -> Automatic)&, $CortexPlot3DOptions];
CortexPlot3D[mesh_?CorticalMeshQ, opts:OptionsPattern[]] := With[
  {Opt = Function[{name},
     Fold[
       Replace,
       OptionValue[name],
       {Automatic :> Replace[name, Options[mesh]],
        (name|Automatic) :> Replace[name, $CortexPlot3DOptions]}]],
   U = VertexList[mesh],
   X = VertexCoordinates[mesh],
   F = Transpose@VertexIndex[mesh, FaceListTr[mesh]],
   vnorms = Replace[
     PropertyValue[{mesh, VertexList}, VertexNormals],
     $Failed :> Replace[
       PropertyValue[{mesh, VertexList}, "VertexNormals"],
       $Failed :> VertexNormals[mesh]]]},
  With[
    {vcolors = GetVertexColors[
       mesh,
       Opt[VertexColors],
       Opt[ColorFunction]]},
    With[
      {vfn = Opt[VertexRenderingFunction],
       efn = Opt[EdgeRenderingFunction],
       ffn = Opt[FaceRenderingFunction],
       vprop = If[Or[(vfn =!= Automatic && vfn =!= None),
                     (efn =!= Automatic && efn =!= None),
                     (ffn =!= Automatic && ffn =!= None)],
         GetProperties[VertexList],
         None]},
      WithOptions[
        Graphics3D[
          GraphicsComplex[
            VertexCoordinates[mesh],
            {Which[
               ffn === None, {},
               ffn === Automatic, {EdgeForm[], Gray, Polygon[F]},
               True, {
                 EdgeForm[], Gray,
                 MapThread[
                   ffn@Append[#1, "VertexColors" -> #2]&,
                   {Normal@FaceDataset[mesh], vcolors[[#]]& /@ F}]}],
             Which[
               efn === None || efn === Automatic, {},
               efn === Line, {
                 Thin,
                 Map[
                   Function@Line[#1, VertexColors -> vcolors[[#1]]],
                   Transpose@VertexIndex[mesh, EdgePairsTr[mesh]]]},
               True, efn /@ Normal@EdgeDataset[mesh]],
             If[vfn === None || vfn === Automatic,
               {},
               MapThread[
                 vfn@Append[#1,"VertexColor" -> #2]&,
                 {Normal@VertexDataset[mesh], vcolors}]]},
            VertexColors -> If[(ffn === Automatic || ffn === None) && vcolors =!= None, 
              vcolors,
              None],
            VertexNormals -> vnorms]],
        Join[{opts}, Options[mesh], $CortexPlot3DOptions]]]]];
Protect[CortexPlot3D, GetVertexColors];


(* #CortexPlot ************************************************************************************)
Options[CortexPlot] = Map[(#[[1]] -> Automatic)&, Options[CorticalMap]];
CortexPlot[mesh_?CorticalMapQ, opts:OptionsPattern[]] := With[
  {Opt = Function[{name},
     Fold[
       Replace,
       OptionValue[name],
       {Automatic :> Replace[name, Options[mesh]],
        (name|Automatic) :> Replace[name, $CortexPlotOptions]}]],
   U = VertexList[mesh],
   X = VertexCoordinates[mesh],
   F = IndexedFaceList[mesh],
   Q = IndexedEdgePairs[mesh]},
  With[
    {vcolors = GetVertexColors[
       mesh,
       Opt[VertexColors],
       Opt[ColorFunction]]},
    With[
      {vfn = Opt[VertexRenderingFunction],
       efn = Opt[EdgeRenderingFunction],
       ffn = Opt[FaceRenderingFunction]},
      WithOptions[
        Graphics[
          GraphicsComplex[
            VertexCoordinates[mesh],
            {Which[
               ffn === None, {},
               ffn === Automatic, {EdgeForm[], Gray, Polygon[F]},
               True, {
                 EdgeForm[], Gray,
                 MapThread[
                   ffn@Append[#1, "VertexColors" -> #2]&,
                   {Normal@FaceDataset[mesh], vcolors[[#]]& /@ F}]}],
             Which[
               efn === None || efn === Automatic, {},
               efn === Line, {
                 Thin,
                 Map[
                   Function@Line[#1, VertexColors -> vcolors[[#1]]],
                   Transpose@VertexIndex[mesh, EdgePairsTr[mesh]]]},
               True, {
                 Thin, Gray,
                 MapThread[
                   efn@Append[#1, "VertexColors" -> #2]&,
                   {Normal@EdgeDataset[mesh], vcolors[[#]]& /@ Q}]}],
             If[vfn === None || vfn === Automatic,
               {},
               {Gray, PointSize[Tiny],
                MapThread[
                  vfn@Append[#1, "VertexColor" -> #2]&,
                  {Normal@VertexDataset[mesh], vcolors}]}]},
            VertexColors -> If[(ffn === Automatic || ffn === None) && vcolors =!= None, 
              vcolors,
              None]]],
        Join[{opts}, Options[mesh], $CortexPlotOptions]]]]];
CortexPlot[mesh_?CorticalMeshQ, opts:OptionsPattern[]] := With[
  {mopts = "CorticalMap" /. Append[Options[mesh, MetaInformation], "CorticalMap" -> {}],
   usemesh = "SphericalMesh" /. Append[Options[mesh, MetaInformation], "SphericalMesh" -> mesh]},
  With[
    {map = CorticalMap @@ Prepend[
       Join[
         FilterRules[{opts}, Options[CorticalMap]],
         mopts],
       usemesh]},
    If[!CorticalMapQ[map], $Failed, CortexPlot[map, opts]]]];
Protect[CortexPlot];


(* #CorticalColorSchema ***************************************************************************)
CorticalColorSchema[property_ -> {range_, colors_}][assoc_?AssociationQ] := If[
  KeyExistsQ[assoc, property],
  With[
    {val = assoc[property]},
    If[!NumericQ[val],
      $Failed,
      If[MatrixQ[colors],
        Blend[colors, val], 
        Blend[colors, Clip[Rescale[val, range], {0,1}]]]]],
  $Failed];
CorticalColorSchema[property_ -> f:Except[{_, _}]][assoc_?AssociationQ] := If[
  KeyExistsQ[assoc, property],
  With[
    {val = assoc[property]},
    If[!NumericQ[val],
      $Failed,
      f[val]]],
  $Failed];
CorticalColorSchema[f:Except[_Rule]][arg_] := f[arg];
(schema:CorticalColorSchema[prop_ -> _])[arg:Except[_?AssociationQ]] := schema[<|prop -> arg|>];
Protect[CorticalColorSchema];

(* #CorticalColorData *****************************************************************************)
CorticalColorData[unknown_] := $Failed;
(* Some cortical color data schemas that we use... *)
polarAngleColors = <|
    -90 -> Magenta,
    0 -> Blue, 45 -> Darker[Cyan, 1/6], 90 -> Darker[Green],
    135 -> Darker[Yellow, 1/6], 180 -> Red|>;
CorticalColorData[Curvature] = CorticalColorSchema[Curvature -> CorticalCurvatureColor];
CorticalColorData["Curvature"] = CorticalColorSchema["Curvature" -> CorticalCurvatureColor];
CorticalColorData["PolarAngle"] = CorticalColorSchema[
  "PolarAngle" -> {
    {-360, 360},
    Table[
      {k, polarAngleColors@Min[{#, 360-#}]&@Mod[k, 360]},
      {k, -360, 360, 45}]}];
CorticalColorData["PolarAngleLH"] = CorticalColorSchema[
  "PolarAngle" -> {
    {-360, 360},
    Join[
      Table[{k, polarAngleColors[k + 360]}, {k, -360, -180, 45}],
      {{-90, polarAngleColors[-90]}},
      Table[{k, polarAngleColors[k]}, {k, 0, 180, 45}],
      {{270, polarAngleColors[-90]},
       {360, polarAngleColors[0]}}]}];
CorticalColorData["PolarAngleRH"] = CorticalColorSchema[
  "PolarAngle" -> {
    {-360, 360},
    Join[
      {{-360, polarAngleColors[0]},
       {-270, polarAngleColors[-90]}},
      Table[{k, polarAngleColors[-k]}, {k, -180, 0, 45}],
      {{90, polarAngleColors[-90]}},
      Table[{k, polarAngleColors[360 - k]}, {k, 180, 360, 45}]]}];
(*
CorticalColorData["Eccentricity"] = CorticalColorSchema[
  "Eccentricity" -> {
    {0, 90},
    Join[
      {Black, Purple, Red, Yellow, Green},
      Table[Blend[{Green, Cyan}, (u - 20.0)/20.0], {u, 25, 40, 5}],
      Table[Blend[{Cyan, White}, (u - 40.0)/50.0], {u, 45, 90, 5}]]}];
CorticalColorData["EccentricityReverse"] = CorticalColorSchema[
  "Eccentricity" -> {
    {0, 90},
    Join[
      {White, Cyan, Green, Yellow, Red},
      Table[Blend[{Red, Purple}, (u - 20.0)/20.0], {u, 25, 40, 5}],
      Table[Blend[{Purple, Black}, (u - 40.0)/50.0], {u, 45, 90, 5}]]}];
*)
CorticalColorData["Eccentricity"] = CorticalColorSchema[
  "Eccentricity" -> Function[
    Blend[
      {{0/90.0,   Black},
       {2.5/90.0, Purple},
       {5/90.0,   Red},
       {10/90.0,  Yellow},
       {20/90.0,  Green},
       {40/90.0,  Cyan},
       {90/90.0,  White}},
      #/90]]];
CorticalColorData["EccentricityReverse"] = CorticalColorSchema[
  "Eccentricity" -> Function[
    Blend[
      {{0/90.0,   White},
       {2.5/90.0, Cyan},
       {5/90.0,   Green},
       {10/90.0,  Yellow},
       {20/90.0,  Red},
       {40/90.0,  Purple},
       {90/90.0,  Black}},
      #/90]]];
CorticalColorData["PRFSize"] = CorticalColorSchema[
  "PRFSize" -> Function[
    Blend[
      {{0/90.0,   Black},
       {1.25/90,  Purple},
       {2.5/90.0, Blue},
       {5/90.0,   Cyan},
       {10/90.0,  Green},
       {25/90.0,  Yellow},
       {50/90.0,  White}},
      #/90]]];

(* #ColorCortex ***********************************************************************************)
SetAttributes[ColorCortex, HoldAll];
ColorCortex[instructions___] := With[
  {slots = Cases[Hold[instructions], Slot[s_ /; !IntegerQ[s]] :> s, Infinity]},
  With[
    {syms = Table[TemporarySymbol[], {Length@slots}]},
    ReplacePart[
      {Hold[
         Evaluate @ MapThread[
           Function[{sym, slot},
             If[StringQ[slot],
               Hold[sym, Replace[slot, #3]],
               Hold[sym, Replace[Replace[slot, #3], slot :> Replace[ToString[slot], #3]]]]],
           {syms, slots}],
         Evaluate @ With[
           {instr = ReplaceAll[
              Hold[instructions],
              MapThread[
                Function[{slot, sym}, Slot[slot] -> sym],
                {slots, syms}]]},
           Block[
             {result},
             Hold[
               {result = CompoundExpression @@ instr},
               Which[
                 ColorQ[result] || Head[result] === Opacity, result,
                 CorticalColorData[result] =!= $Failed, CorticalColorData[result][#1, #2, #3],
                 result === $Failed || result === None, RGBColor[1,1,1,0],
                 True, RGBColor[1,1,1,0]]]]]]},
        {{1,1,_,0} -> Set,
         {1,2,0} -> With,
         {1,0} -> With,
         0 -> Function}]]];
Protect[ColorCortex];

(* #CortexResample ********************************************************************************)
Options[CortexResample] = {Method -> Interpolation, Properties -> All, Indeterminate -> None};
(* These private functions return a pair of matrices: one of indices and one of weights, for the
 * given matrix and mesh.
 *)
CortexResampleNearest[X_?MatrixQ, from_?CorticalObjectQ] := With[
  {near = Nearest[from]},
  With[
    {idcs = near[If[Length[X] < 4 && Length@First[X] > 3, Transpose[X], X], 1]},
    {VertexIndex[from, #]&/@Transpose[idcs], ConstantArray[1, Length[idcs]]}]];
CortexResampleTrilinear[X_?MatrixQ, from_?CorticalObjectQ] := With[
  {dat = NearestFace[from, X, NearestPoints -> True],
   Fx = Transpose[FaceCoordinatesTr[from], {3, 1, 2}]},
  With[
    {px = dat[[2]],
     Ix = Fx[[All, dat[[1]]]]},
    With[
      {Px = {Ix[[1]] - px, Ix[[2]] - px, Ix[[3]] - px} // Function@If[CorticalMapQ[from],
         Map[Transpose@Append[Transpose[#], ConstantArray[0.0, Length[#]]]&, #],
         #]},
      With[
        {weights = (# / ConstantArray[Total[#], Length[#]])&@Map[
           RowNorms,
           MapThread[Cross, Px[[#]]]& /@ {{2,3}, {1,3}, {1,2}}]},
        {IndexedFaceListTr[from][[All, dat[[1]]]], weights}]]]];
CortexResample[a_ /; MatrixQ[a, NumericQ], b_?CorticalObjectQ, opts:OptionsPattern[]] := Check[
  With[
    {propNames = Select[
       Replace[
         OptionValue[Properties],
         {All :> VertexPropertyList[b],
          p:Except[_List] :> {p}}],
       # =!= VertexCoordinates &],
     fill = OptionValue[Indeterminate],
     method = Replace[
       OptionValue[Method],
       {Nearest|"Nearest" -> CortexResampleNearest,
        Interpolation|"Interpolation"|"Trilinear" -> CortexResampleTrilinear,
        _ :> Message[
          CortexResample::badarg,
          "Method arg should be \"Nearest\" or \"Trilinear\""]}],
     X0 = Which[
       Length[a] == Length@VertexCoordinatesTr[b], a,
       Length@First[a] == Length@VertexCoordinatesTr[b], Transpose[a],
       True, Message[CortexResample::badarg, "destination-matrix dimensions are incorrect"]]},
    With[
      {dat = method[X0, b],
       props = (# /. b)& /@ propNames},
      With[
        {res = If[Length@First[dat] == 1,
           MapThread[
             #1 -> #2[[dat[[1,1]]]] &,
             {propNames, props}],
           MapThread[
             Function@With[
               {name = #1, vals0 = Table[#2[[idcs]], {idcs, dat[[1]]}]},
               With[
                 {where = Boole@Map[NumericQ, vals0, {2}],
                  vals = Replace[vals0, Except[_?NumericQ] -> 0.0, {2}]},
                 With[
                   {weights = dat[[2]] * where},
                   With[
                     {tot = Total[weights]},
                     With[
                       {zeros = 1 - Unitize[tot]},
                       name -> ReplacePart[
                         Total[vals * weights] / (zeros + tot),
                         Position[zeros, 1, {1}] -> fill]]]]]],
             {propNames, props}]]},
        If[StringQ@OptionValue[Properties], res[[1,2]], res]]]],
  $Failed];
CortexResample[a_?CorticalObjectQ, b_?CorticalObjectQ, opts:OptionsPattern[]] := CortexResample[
  VertexCoordinatesTr[a], b, opts];
CortexResample[Rule[a_?CorticalObjectQ, b_], opts:OptionsPattern[]] := CortexResample[
  b, a, opts];
Protect[CortexResample];

(* #CortexResampleOnto ****************************************************************************)
Options[CortexResampleOnto] = Options[CortexResample];
CortexResampleOnto[a_?CorticalObjectQ, b_?CorticalObjectQ, opts:OptionsPattern[]] := With[
  {resamp = Normal@CortexResample[a, b, opts],
   p = OptionValue[Properties]},
  If[StringQ[p],
    SetVertexProperties[a, p -> resamp],
    Fold[SetVertexProperties, a, resamp]]];
CortexResampleOnto[Rule[a_, b_], opts:OptionsPattern[]] := CortexResampleOnto[
  b, a, opts];
Protect[CortexResampleOnto];

(* #CorticalLabelQ ********************************************************************************)
CorticalLabelQ[mesh_?CorticalObjectQ, label_] := False;
CorticalLabelQ[map_?CorticalMapQ, label_] := Which[
  VectorQ[label, IntegerQ], Which[
    Length[Complement[label, VertexList[map]]] == 0, True,
    Length[label] == VertexCount[map] && Length[Complement[label, {1, 0}]] == 0, True,
    True, False],
  ArrayQ[label, 1|2], Which[
    Length[label] == VertexCount[map] && Complement[label, {0,1,True,False,0.0,1.0}] == 0, True,
    AllTrue[label, Or[Head[#]===UndirectedEdge, VectorQ[label] && Length[label] == 2]&], With[
      {fix = ReplaceAll[label, UndirectedEdge -> List]},
      And[Complement[Join @@ fix, VertexList[map]] == {},
          FindCycle[Graph[fix], {Length[fix]}] != {}]],
    True, False],
  BoundaryMeshQ[label] && RegionDimension[label] == 2, True,
  True, False];
CorticalLabelQ[mesh_?CorticalMeshQ, label_] := Which[
  VectorQ[label, IntegerQ], Which[
    Length[Complement[label, VertexList[mesh]]] == 0, True,
    Length[label] == VertexCount[mesh] && Length[Complement[label, {1, 0}]] == 0, True,
    True, False],
  ArrayQ[label, 1|2], Which[
    Length[label] == VertexCount[mesh] && Complement[label, {0,1,True,False,0.0,1.0}] == 0, True,
    And[
      AllTrue[Rest[label], Or[Head[#]===UndirectedEdge, VectorQ[label] && Length[label] == 2]&],
      IntegerQ[label[[1]]],
      IntegerQ[VertexIndex[mesh, label[[1]]]]], With[
        {fix = ReplaceAll[label, UndirectedEdge -> List]},
        And[Complement[Join @@ fix, VertexList[mesh]] == {},
            FindCycle[Graph[fix], {Length[fix]}] != {}]],
    True, False],
  BoundaryMeshQ[label] && RegionDimension[label] == 3, True,
  True, False];
Protect[CorticalLabelQ];

(* #CorticalLabelMask *****************************************************************************)
CorticalLabelMask[map_?CorticalMapQ, label_] := Which[
  VectorQ[label, IntegerQ], Which[
    Length[Complement[label, VertexList[map]]] == 0, SparseArray[
      Thread[VertexIndex[map, label] -> 1],
      {VertexCount[map]},
      0],
    Length[label] == VertexCount[map] && Length[Complement[label, {1, 0}]] == 0, label,
    True, (Message[CorticalLabelMask::notlab]; $Failed)],
  ArrayQ[label, 1|2], Which[
    Length[label] == VertexCount[map] && Complement[label, {0,1,True,False,0.0,1.0}] == 0, label,
    AllTrue[label, Or[Head[#]===UndirectedEdge, VectorQ[label] && Length[label] == 2]&], With[
      {fix = ReplaceAll[label, UndirectedEdge -> List]},
      If[And[Complement[Join @@ fix, VertexList[map]] == {},
             FindCycle[Graph[fix], {Length[fix]}] != {}],
        With[
          {pgonfn = SignedRegionDistance[
             Polygon[Part[VertexCoordinates[map], #]& /@ VertexIndex[map, fix]]],
           X = VertexCoordinates[map]},
          SparseArray[
            Thread[Select[Range[Length@X], pgonfn[#] <= 0&] -> 1],
            {Length@X},
            0]],
        (Message[CorticalLabelMask::notlab]; $Failed)]],
    True, (Message[CorticalLabelMask::notlab]; $Failed)],
  BoundaryMeshQ[label] && RegionDimension[label] == 2, With[
    {pgonfn = SignedRegionDistance[label],
     X = VertexCoordinates[map]},
    SparseArray[
      Thread[Select[Range[Length@X], pgonfn[#] <= 0&] -> 1],
      {Length@X},
      0]],
  True, (Message[CorticalLabelMask::notlab]; $Failed)];
CorticalLabelMask[mesh_?CorticalMeshQ, label_] := Which[
  VectorQ[label, IntegerQ], Which[
    Length[Complement[label, VertexList[mesh]]] == 0, SparseArray[
      Thread[VertexIndex[mesh, label] -> 1],
      {VertexCount[mesh]},
      0],
    Length[label] == VertexCount[mesh] && Length[Complement[label, {1, 0}]] == 0, label,
    True, (Message[CorticalLabelMask::notlab]; $Failed)],
  ArrayQ[label, 1|2], Which[
    Length[label] == VertexCount[mesh] && Complement[label, {0,1,True,False,0.0,1.0}] == 0, label,
    AllTrue[label, Or[Head[#]===UndirectedEdge, VectorQ[label] && Length[label] == 2]&], With[
      {fix = ReplaceAll[label, UndirectedEdge -> List]},
      If[And[Complement[Join @@ fix, VertexList[mesh]] == {},
             FindCycle[Graph[fix], {Length[fix]}] != {}],
        With[
          {pgonfn = SignedRegionDistance[
             Polygon[Part[VertexCoordinates[mesh], #]& /@ VertexIndex[mesh, fix]]],
           X = VertexCoordinates[mesh]},
          SparseArray[
            Thread[Select[Range[Length@X], pgonfn[#] <= 0&] -> 1],
            {Length@X},
            0]],
        (Message[CorticalLabelMask::notlab]; $Failed)]],
    True, (Message[CorticalLabelMask::notlab]; $Failed)],
  BoundaryMeshQ[label] && RegionDimension[label] == 2, With[
    {pgonfn = SignedRegionDistance[label],
     X = VertexCoordinates[mesh]},
    SparseArray[
      Thread[Select[Range[Length@X], pgonfn[#] <= 0&] -> 1],
      {Length@X},
      0]],
  True, (Message[CorticalLabelMask::notlab]; $Failed)];

(* #CorticalLabelVertexList ***********************************************************************)
CorticalLabelVertexList[map_?CorticalMapQ, label_] := Which[
  VectorQ[label, IntegerQ], Which[
    Length[Complement[label, VertexList[map]]] == 0, label,
    Length[label] == VertexCount[map] && Length[Complement[label, {1, 0}]] == 0, Pick[
      VertexList[map], label, 1],
    True, (Message[CorticalLabelVertexList::notlab]; $Failed)],
  ArrayQ[label, 1|2], Which[
    Length[label] == VertexCount[map] && Complement[label, {0,1,True,False,0.0,1.0}] == 0, Pick[
      VertexList[map], label, 1|1.0|True],
    AllTrue[label, Or[Head[#]===UndirectedEdge, VectorQ[label] && Length[label] == 2]&], With[
      {fix = ReplaceAll[label, UndirectedEdge -> List]},
      If[And[Complement[Join @@ fix, VertexList[map]] == {},
             FindCycle[Graph[fix], {Length[fix]}] != {}],
        With[
          {pgonfn = SignedRegionDistance[
             Polygon[Part[VertexCoordinates[map], #]& /@ VertexIndex[map, fix]]],
           X = VertexCoordinates[map]},
          Pick[VertexList[map], X, x_ /; pgonfn[x] <= 0]],
        (Message[CorticalLabelVertexList::notlab]; $Failed)]],
    True, (Message[CorticalLabelVertexList::notlab]; $Failed)],
  BoundaryMeshQ[label] && RegionDimension[label] == 2, With[
    {pgonfn = SignedRegionDistance[label],
     X = VertexCoordinates[map]},
    Pick[VertexList[map], X, x_ /; pgonfn[x] <= 0]],
  True, (Message[CorticalLabelVertexList::notlab]; $Failed)];
CorticalLabelVertexList[mesh_?CorticalMeshQ, label_] := Which[
  VectorQ[label, IntegerQ], Which[
    Length[Complement[label, VertexList[mesh]]] == 0, label,
    Length[label] == VertexCount[mesh] && Length[Complement[label, {1, 0}]] == 0, Pick[
      VertexList[mesh], label, 1],
    True, (Message[CorticalLabelVertexList::notlab]; $Failed)],
  ArrayQ[label, 1|2], Which[
    Length[label] == VertexCount[mesh] && Complement[label, {0,1,True,False,0.0,1.0}] == {}, Pick[
      VertexList[mesh], label, 1|1.0|True],
    AllTrue[label, Or[Head[#]===UndirectedEdge, VectorQ[label] && Length[label] == 2]&], With[
      {fix = ReplaceAll[label, UndirectedEdge -> List]},
      If[And[Complement[Join @@ fix, VertexList[mesh]] == {},
             FindCycle[Graph[fix], {Length[fix]}] != {}],
        With[
          {pgonfn = SignedRegionDistance[
             Polygon[Part[VertexCoordinates[mesh], #]& /@ VertexIndex[mesh, fix]]],
           X = VertexCoordinates[mesh]},
          Pick[VertexList[mesh], X, x_ /; pgonfn[x] <= 0]],
        (Message[CorticalLabelVertexList::notlab]; $Failed)]],
    True, (Message[CorticalLabelVertexList::notlab]; $Failed)],
  BoundaryMeshQ[label] && RegionDimension[label] == 2, With[
    {pgonfn = SignedRegionDistance[label],
     X = VertexCoordinates[mesh]},
    Pick[VertexList[mesh], X, x_ /; pgonfn[x] <= 0]],
  True, (Message[CorticalLabelVertexList::notlab]; $Failed)];

(* #LabelVertexCoordinatesTr **********************************************************************)
LabelVertexCoordinatesTr[cortex_?CorticalObjectQ, name_] := Check[
  Part[
    VertexCoordinatesTr[cortex],
    All,
    VertexIndex[cortex, CorticalLabelVertexList[cortex, name]]],
  $Failed];
LabelVertexCoordinatesTr[sub_, mesh_, hemi_, name_] := Check[
  With[
    {cortex = Cortex[sub, hemi, mesh],
     U = Normal@LabelVertexList[sub, hemi, name]},
    If[MissingQ[U],
      U,
      Part[
        VertexCoordinatesTr[cortex],
        All,
        VertexIndex[cortex, U]]]],
  $Failed];
Protect[LabelVertexCoordinatesTr];

(* #LabelVertexCoordinates ************************************************************************)
LabelVertexCoordinates[cortex_?CorticalObjectQ, name_] := Check[
  Transpose[LabelVertexCoordinatesTr[cortex, name]],
  $Failed];
LabelVertexCoordinates[sub_, mesh_, hemi_, name_] := Check[
  Transpose[LabelVertexCoordinatesTr[sub, mesh, hemi, name]],
  $Failed];
Protect[LabelVertexCoordinates];

(* #LabelEdgePairsTr ******************************************************************************)
LabelEdgePairsTr[cortex_?CorticalObjectQ, name_] := Check[
  With[
    {U = CorticalLabelVertexList[cortex, name]},
    With[
      {UE = VertexEdgeList[cortex][[VertexIndex[cortex, U]]],
       Et = EdgePairsTr[cortex]},
      Et[[All, Select[Tally[Join@@UE], Last[#] == 2&][[All, 1]]]]]],
  $Failed];
LabelEdgePairsTr[sub_, hemi_, name_] := Check[
  LabelEdgePairsTr[Cortex[sub, hemi, Automatic], name],
  $Failed];
Protect[LabelEdgePairsTr];

(* #LabelEdgePairs ********************************************************************************)
LabelEdgePairs[cortex_?CorticalObjectQ, name_] := Check[
  Transpose @ LabelEdgePairsTr[cortex, name],
  $Failed];
LabelEdgePairs[sub_, hemi_, name_] := Check[
  Transpose @ LabelEdgePairsTr[sub, hemi, name],
  $Failed];
Protect[LabelEdgePairs];

(* #LabelEdgeList *********************************************************************************)
LabelEdgeList[sub_, hemi_, name_] := Check[
  Apply[UndirectedEdge, #]& /@ Transpose[LabelEdgePairsTr[sub, hemi, name]],
  $Failed];
Protect[LabelEdgeList];

(* #LabelFaceListTr ******************************************************************************)
LabelFaceListTr[cortex_?CorticalObjectQ, name_] := Check[
  With[
    {U = CorticalLabelVertexList[cortex, name]},
    With[
      {UF = VertexFaceList[cortex][[VertexIndex[cortex, U]]],
       Ft = FaceListTr[cortex]},
      Ft[[All, Select[Tally[Join@@UF], Last[#] == 3&][[All, 1]]]]]],
  $Failed];
LabelFaceListTr[sub_, hemi_, name_] := Check[
  With[
    {U = Normal@LabelVertexList[sub, hemi, name],
     cortex = Cortex[sub, hemi, Automatic]},
    If[MissingQ[U],
      U,
      With[
        {UF = VertexFaceList[cortex][[VertexIndex[cortex, U]]],
         Ft = FaceListTr[cortex]},
        Ft[[All, Select[Tally[Join@@UF], Last[#] == 3&][[All, 1]]]]]]],
  $Failed];
Protect[LabelFaceListTr];

(* #LabelFaceList *********************************************************************************)
LabelFaceList[cortex_?CorticalObjectQ, name_] := Check[
  Replace[LabelFaceListTr[cortex, name], m:Except[_?MissingQ] :> Transpose[m]],
  $Failed];
LabelFaceList[sub_, hemi_, name_] := Check[
  Replace[LabelFaceListTr[sub, hemi, name], m:Except[_?MissingQ] :> Transpose[m]],
  $Failed];
Protect[LabelFaceList];

(* #LabelBoundaryEdgePairsTr **********************************************************************)
Options[LabelBoundaryEdgePairsTr] = {Method -> "Longest"};
LabelBoundaryEdgePairsTrFromFaceList[Ft_, name_, methodOpt_] := Check[
  With[
    {accum = Switch[methodOpt,
       "Longest", If[#1 == {} || Length@First[#2] > Length@First[#1], #2, #1]&,
       All, Append[#1, #2]&,
       _, Message[
         LabelBoundaryEdgePairs::badarg,
         "Method must be All or \"Longest\""]]},
    If[Length[Ft] == 0,
      {},
      With[
        {allE = Sort /@ Transpose[{Join@@Ft, Join[Ft[[2]], Ft[[3]], Ft[[1]]]}]},
        With[
          {pairs = Select[Tally[allE], Last[#] == 1&][[All, 1]],
           sym = Unique["tag"]},
          First@NestWhile[
            Function@With[
              {oldpath = #[[1]], edges = #[[2]]},
              With[
                {newpath = {#, RotateLeft[#]}&@FindShortestPath[
                   Graph@Rest[edges],
                   edges[[1,1]],
                   edges[[1,2]]]},
                If[Length@First[newpath] < 2,
                  {oldpath, Rest[edges]},
                  {accum[oldpath, newpath], Complement[edges, Sort /@ Transpose[newpath]]}]]],
            {{}, pairs},
            Length[#[[2]]] > 2&]]]]],
  $Failed];
LabelBoundaryEdgePairsTr[cortex_?CorticalObjectQ, name_, OptionsPattern[]] := Check[
  LabelBoundaryEdgePairsTrFromFaceList[LabelFaceListTr[cortex, name], name, OptionValue[Method]],
  $Failed];
LabelBoundaryEdgePairsTr[sub_, hemi_, name_, OptionsPattern[]] := Check[
  LabelBoundaryEdgePairsTrFromFaceList[LabelFaceListTr[sub, hemi, name], name, OptionValue[Method]],
  $Failed];
Protect[LabelBoundaryEdgePairsTr];

(* #LabelBoundaryEdgePairs ************************************************************************)
Options[LabelBoundaryEdgePairs] = Options[LabelBoundaryEdgePairsTr];
LabelBoundaryEdgePairs[cortex_?CorticalObjectQ, name_, opts:OptionsPattern[]] := Check[
  With[
    {p = LabelBoundaryEdgePairsTr[cortex, name, opts]},
    If[p == {}, p, Transpose[p]]],
  $Failed];
LabelBoundaryEdgePairs[sub_, hemi_, name_, opts:OptionsPattern[]] := Check[
  With[
    {p = LabelBoundaryEdgePairsTr[sub, hemi, name, opts]},
    If[p == {}, p, Transpose[p]]],
  $Failed];
Protect[LabelBoundaryEdgePairs];

(* #LabelBoundaryEdgeList *************************************************************************)
Options[LabelBoundaryEdgeList] = Options[LabelBoundaryEdgePairsTr];
LabelBoundaryEdgeList[cortex_?CorticalObjectQsub_, name_, opts:OptionsPattern[]] := Check[
  Apply[UndirectedEdge, #]& /@ LabelBoundaryEdgePairs[cortex, name, opts],
  $Failed];
LabelBoundaryEdgeList[sub_, hemi_, name_, opts:OptionsPattern[]] := Check[
  Apply[UndirectedEdge, #]& /@ LabelBoundaryEdgePairs[sub, hemi, name, opts],
  $Failed];
Protect[LabelBoundaryEdgeList];

(* #LabelBoundaryVertexList ***********************************************************************)
Options[LabelBoundaryVertexList] = Options[LabelBoundaryEdgePairsTr];
LabelBoundaryVertexList[cortex_?CorticalObjectQ, name_, opts:OptionsPattern[]] := Check[
  With[
    {e = LabelBoundaryEdgePairsTr[cortex, name, opts]},
    Which[
      !ListQ[e], $Failed,
      Length[e] == 0, {},
      ArrayQ[e, 2], e[[1]],
      True, e[[All, 1]]]],
  $Failed];
LabelBoundaryVertexList[sub_, hemi_, name_, opts:OptionsPattern[]] := Check[
  With[
    {e = LabelBoundaryEdgePairsTr[sub, hemi, name, opts]},
    Which[
      !ListQ[e], $Failed,
      Length[e] == 0, {},
      ArrayQ[e, 2], e[[1]],
      True, e[[All, 1]]]],
  $Failed];
Protect[LabelBoundaryVertexList];

(* #LabelBoundaryVertexCoordinatesTr **************************************************************)
Options[LabelBoundaryVertexCoordinatesTr] = Options[LabelBoundaryEdgePairsTr];
LabelBoundaryVertexCoordinatesTr[cortex_?CorticalObjectQ, name_, opts:OptionsPattern[]] := Check[
  With[
    {vs = VertexIndex[cortex, LabelBoundaryVertexList[cortex, name, opts]],
     X = VertexCoordinatesTr[cortex]},
    If[ArrayQ[vs, 1],
      X[[All, vs]],
      X[[All, #]]& /@ vs]],
  $Failed];
LabelBoundaryVertexCoordinatesTr[sub_, mesh_, hemi_, name_, opts:OptionsPattern[]] := Check[
  LabelBoundaryVertexCoordinatesTr[Cortex[sub, hemi, mesh], name, opts],
  $Failed];
Protect[LabelBoundaryVertexCoordinatesTr];

(* #LabelBoundaryVertexCoordinates ****************************************************************)
Options[LabelBoundaryVertexCoordinates] = Options[LabelBoundaryEdgePairsTr];
LabelBoundaryVertexCoordinates[cortex_?CorticalObjectQ, name_, opts:OptionsPattern[]] := Check[
  Transpose @ LabelBoundaryVertexCoordinatesTr[cortex, name, opts],
  $Failed];
LabelBoundaryVertexCoordinates[sub_, mesh_, hemi_, name_, opts:OptionsPattern[]] := Check[
  Transpose @ LabelBoundaryVertexCoordinatesTr[sub, mesh, hemi, name, opts],
  $Failed];
Protect[LabelBoundaryVertexCoordinates];

(* #OccipitalPole *********************************************************************************)
OccipitalPole[sub_, hemi_?HemiQ, mesh_] := Check[
  VertexCoordinatesTr[Cortex[sub, hemi, mesh]][[All, OccipitalPoleIndex[sub, hemi]]],
  $Failed];
OccipitalPole[sub_, m:Except[_?HemiQ], h_?HemiQ] := OccipitalPole[sub, h, m];
OccipitalPole[sub_, hemi_] := OccipitalPole[sub, hemi, Automatic];
Protect[OccipitalPole];

(* #MapBoundaryVertexList *************************************************************************)
MapBoundaryVertexList[map_?CorticalMapQ] := With[
  {bounds = Indices[EdgeFaceList[map], {x_Integer}],
   edges = EdgePairs[map]},
  FindShortestPath[
   Graph[edges[[Rest@bounds]]],
   edges[[bounds[[1]], 2]],
   edges[[bounds[[1]], 1]]]];
Protect[MapBoundaryVertexList];

(* #MapBoundaryEdgePairsTr ************************************************************************)
MapBoundaryEdgePairsTr[map_?CorticalMapQ] := With[
  {vlist = MapBoundaryVertexList[map]},
  {vlist, RotateLeft[vlist]}];
Protect[MapBoundaryEdgePairsTr];

(* #MapBoundaryEdgePairs **************************************************************************)
MapBoundaryEdgePairs[map_?CorticalMapQ] := Transpose @ MapBoundaryEdgePairsTr[map];
Protect[MapBoundaryEdgePairs];

(* #MapBoundaryEdgeList ***************************************************************************)
MapBoundaryEdgeList[map_?CorticalMapQ] := Map[
  Apply[UndirectedEdge, #]&,
  MapBoundaryEdgePairs[map]];
Protect[MapBoundaryEdgeList];

(* #MapBoundaryFaceList ***************************************************************************)
MapBoundaryFaceList[map_?CorticalMapQ] := With[
  {bounds = EdgeIndex[map, MapBoundaryEdgePairs[map]],
   EFL = EdgeFaceList[map]},
  Join@@EFL[[bounds]]];
Protect[MapBoundaryFaceList];

(* #MapBoundaryFaceListTr *************************************************************************)
MapBoundaryFaceListTr[map_?CorticalMapQ] := Transpose @ MapBoundaryFaceList[map];
Protect[MapBoundaryFaceListTr];

End[];
EndPackage[];
