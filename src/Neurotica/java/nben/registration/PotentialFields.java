////////////////////////////////////////////////////////////////////////////////////////////////////
// PotentialFields.java
//
// The nben.registration namespace contains functions related to the registration of surface
// meshes to models defined on the cortical surface; it is designed to work with the Mathematica
// Neurotica`Registration namespace.
//
// Copyright (C) 2016 by Noah C. Benson.
// This file is part of the Neurotica library.
//
// This program is free software: you can redistribute it and/or modify it under the terms of the
// GNU General Public License as published by the Free Software Foundation, either version 3 of the
// License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
// without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See
// the GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along with this program.
// If not, see <http://www.gnu.org/licenses/>.

package nben.registration;

import java.util.concurrent.ExecutorService;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import java.util.HashSet;
import java.util.Iterator;

/** The PotentialFields class is basically a static container for certain data useful to the
 *  various registration potential field classes.
 *
 *  @author Noah C. Benson
 */
public final class PotentialFields {

   /** PotentialFields.newHarmonicEdgePotential(s, q, E, X) yields an EdgePotential object with a 
    *  HarmonicFunction form using scale parameter s/m and shape parameter q where where m is the
    *  number of edges. If the parameter E is actually a list of faces, it is interpreted 
    *  automatically.
    */
   public static EdgePotential newHarmonicEdgePotential(double[] scale, double[] shape,
                                                        int[][] E, double[][] X) {
      if (E.length == 3)
         return newHarmonicEdgePotential(scale, shape, Util.facesToEdges(E), X);
      int n = scale.length;
      HarmonicFunction[] fs = new HarmonicFunction[n];
      for (int i = 0; i < n; ++i)
         fs[i] = new HarmonicFunction(scale[i] / n, shape[i]);
      return new EdgePotential(fs, E, X);
   }
   public static EdgePotential newHarmonicEdgePotential(double[] scale, double shape,
                                                        int[][] E, double[][] X) {
      if (E.length == 3)
         return newHarmonicEdgePotential(scale, shape, Util.facesToEdges(E), X);
      int n = scale.length;
      HarmonicFunction[] fs = new HarmonicFunction[n];
      for (int i = 0; i < n; ++i)
         fs[i] = new HarmonicFunction(scale[i] / n, shape);
      return new EdgePotential(fs, E, X);
   }
   public static EdgePotential newHarmonicEdgePotential(double scale, double[] shape,
                                                        int[][] E, double[][] X) {
      if (E.length == 3)
         return newHarmonicEdgePotential(scale, shape, Util.facesToEdges(E), X);
      int n = shape.length;
      HarmonicFunction[] fs = new HarmonicFunction[n];
      for (int i = 0; i < n; ++i)
         fs[i] = new HarmonicFunction(scale / n, shape[i]);
      return new EdgePotential(fs, E, X);
   }
   public static EdgePotential newHarmonicEdgePotential(double[] scale, 
                                                        int[][] E, double[][] X) {
      if (E.length == 3)
         return newHarmonicEdgePotential(scale, 2.0, Util.facesToEdges(E), X);
      int n = scale.length;
      HarmonicFunction[] fs = new HarmonicFunction[n];
      for (int i = 0; i < n; ++i)
         fs[i] = new HarmonicFunction(scale[i] / n, 2.0);
      return new EdgePotential(fs, E, X);
   }   public static EdgePotential newHarmonicEdgePotential(double scale, double shape,
                                                        int[][] E, double[][] X) {
      if (E.length == 3) {
         int[][] tmp = Util.facesToEdges(E);
         return newHarmonicEdgePotential(scale, shape, tmp, X);
      } else
         return new EdgePotential(new HarmonicFunction(scale / E[0].length, shape), E, X);
   }
   public static EdgePotential newHarmonicEdgePotential(double scale, int[][] E, double[][] X) {
      if (E.length == 3) {
         int[][] tmp = Util.facesToEdges(E);
         return newHarmonicEdgePotential(scale, tmp, X);
      } else
         return new EdgePotential(new HarmonicFunction(scale / E[0].length), E, X);
   }
   public static EdgePotential newHarmonicEdgePotential(int[][] E, double[][] X) {
      if (E.length == 3)
         return newHarmonicEdgePotential(Util.facesToEdges(E), X);
      else
         return new EdgePotential(new HarmonicFunction(1.0 / E[0].length), E, X);
   }
   /** PotentialFields.newLJEdgePotential(s, q, E, X) yields an EdgePotential object with a 
    *  LennardJonesFunction form using scale parameter s/m and shape parameter q where where m is
    *  the number of edges.
    */
   public static EdgePotential newLJEdgePotential(double scale, double shape,
                                                  int[][] E, double[][] X) {
      if (E.length == 3) 
         return newLJEdgePotential(scale, shape, Util.facesToEdges(E), X);
      else
         return new EdgePotential(new LennardJonesFunction(scale / E[0].length, shape), E, X);
   }
   /** PotentialFields.newHarmonicAnglePotential(s, q, E, X) yields an EdgePotential object with a 
    *  LennardJonesFunction form using scale parameter s/m and shape parameter q where where m is
    *  the number of edges.
    */
   public static AnglePotential newHarmonicAnglePotential(double scale, double shape,
                                                          int[][] T, double[][] X) {
      return new AnglePotential(new HarmonicFunction(scale / (3.0 * T[0].length), shape), 
                                Util.facesToAngles(T), X);
   }
   /** PotentialFields.newLJAnglePotential(s, q, T, X) yields an AnglePotential object with a 
    *  LennardJonesFunction form using scale parameter s/m and shape parameter q where where m is
    *  the number of angles (triangles times 3).
    */
   public static AnglePotential newLJAnglePotential(double scale, double shape,
                                                    int[][] T, double[][] X) {
      return new AnglePotential(new LennardJonesFunction(scale / (3.0 * T[0].length), shape), 
                                Util.facesToAngles(T), X);
   }
   public static AnglePotential newLJAnglePotential(double scale,
                                                    int[][] T, double[][] X) {
      return new AnglePotential(new LennardJonesFunction(scale / (3.0 * T[0].length), 2.0), 
                                Util.facesToAngles(T), X);
   }
   /** PotentialFields.newAngleWellPotential(s, q, T, X) yields an AnglePotential object with an
    *  InfinteWellFunction form using scale parameter s/m and shape parameter q where where m is
    *  the number of angles (triangles times 3). Angle wells always have a center of Pi/4 and a 
    *  width of Pi/4.
    */
   public static AnglePotential newAngleWellPotential(double scale, 
                                                      double min, double max,
                                                      double shape,
                                                      int[][] T, double[][] X) {
      return new AnglePotential(new InfiniteWellFunction(scale / (3.0 * T[0].length), 
                                                         min, max, shape),
                                Util.facesToAngles(T), X);
   }
   public static AnglePotential newAngleWellPotential(double scale, 
                                                      double min, double max,
                                                      int[][] T, double[][] X) {
      return new AnglePotential(new InfiniteWellFunction(scale / (3.0 * T[0].length), 
                                                         min, max),
                                Util.facesToAngles(T), X);
   }
   public static AnglePotential newAngleWellPotential(double scale, 
                                                      int[][] T, double[][] X) {
      return new AnglePotential(new InfiniteWellFunction(scale / (3.0 * T[0].length)), 
                                Util.facesToAngles(T), X);
   }
   public static AnglePotential newAngleWellPotential(int[][] T, double[][] X) {
      return new AnglePotential(new InfiniteWellFunction(1.0 / (3.0 * T[0].length)), 
                                Util.facesToAngles(T), X);
   }

   /** newHarmonicAnchorPotential(s, q, vertices, anchorPoints, X) yields an AnchorPotential object
    *  that operates over the given list of vertices, each of which is attracted to the
    *  corresponding anchor point in the (dims x m)-sized matrix, anchorPoints, where m is the
    *  length of the vertices list. The scale used for the harmonic is s/vertices.length and the
    *  shape is q.
    */
   public static AnchorPotential newHarmonicAnchorPotential(double scale, double shape,
                                                            int[] vertices, double[][] points, 
                                                            double[][] X) {
      return new AnchorPotential(new HarmonicFunction(scale/vertices.length, shape),
                                 vertices, points, X);
   }
   public static AnchorPotential newHarmonicAnchorPotential(double scale, 
                                                            int[] vertices, double[][] points, 
                                                            double[][] X) {
      return new AnchorPotential(new HarmonicFunction(scale/vertices.length, 2.0),
                                 vertices, points, X);
   }
   public static AnchorPotential newHarmonicAnchorPotential(int[] vertices, double[][] points, 
                                                            double[][] X) {
      return new AnchorPotential(new HarmonicFunction(1.0/vertices.length, 2.0),
                                 vertices, points, X);
   }

   /** newGaussianAnchorPotentials(k, s, q, vertices, anchorPoints, X) yields an AnchorPotential
    *  object that operates over the given list of vertices, each of which is attracted to the 
    *  corresponding anchor point in the (dims x m)-sized matrix, anchorPoints, where m is the
    *  length of the vertices list. The scale used for the inverted Gaussian potential shape is 
    *  k/vertices.length, the Gaussian's standard deviation is s, and the shape is q.
    */
   public static AnchorPotential newGaussianAnchorPotential(double[] scale, double[] sig, 
                                                            double[] shape,
                                                            int[] vertices, double[][] points,
                                                            double[][] X) {
      GaussianFunction[] fs = new GaussianFunction[vertices.length];
      for (int i = 0; i < fs.length; ++i)
         fs[i] = new GaussianFunction(scale[i]/vertices.length, sig[i], shape[i]);
      return new AnchorPotential(fs, vertices, points, X);
   }
   public static AnchorPotential newGaussianAnchorPotential(double[] scale, double[] sig, 
                                                            double shape,
                                                            int[] vertices, double[][] points,
                                                            double[][] X) {
      GaussianFunction[] fs = new GaussianFunction[vertices.length];
      for (int i = 0; i < fs.length; ++i)
         fs[i] = new GaussianFunction(scale[i]/vertices.length, sig[i], shape);
      return new AnchorPotential(fs, vertices, points, X);
   }
   public static AnchorPotential newGaussianAnchorPotential(double[] scale, double sig, 
                                                            double[] shape,
                                                            int[] vertices, double[][] points,
                                                            double[][] X) {
      GaussianFunction[] fs = new GaussianFunction[vertices.length];
      for (int i = 0; i < fs.length; ++i)
         fs[i] = new GaussianFunction(scale[i]/vertices.length, sig, shape[i]);
      return new AnchorPotential(fs, vertices, points, X);
   }
   public static AnchorPotential newGaussianAnchorPotential(double scale, double[] sig, 
                                                            double[] shape,
                                                            int[] vertices, double[][] points,
                                                            double[][] X) {
      GaussianFunction[] fs = new GaussianFunction[vertices.length];
      for (int i = 0; i < fs.length; ++i)
         fs[i] = new GaussianFunction(scale/vertices.length, sig[i], shape[i]);
      return new AnchorPotential(fs, vertices, points, X);
   }
   public static AnchorPotential newGaussianAnchorPotential(double scale, double sig, 
                                                            double[] shape,
                                                            int[] vertices, double[][] points,
                                                            double[][] X) {
      GaussianFunction[] fs = new GaussianFunction[vertices.length];
      for (int i = 0; i < fs.length; ++i)
         fs[i] = new GaussianFunction(scale/vertices.length, sig, shape[i]);
      return new AnchorPotential(fs, vertices, points, X);
   }
   public static AnchorPotential newGaussianAnchorPotential(double scale, double[] sig, 
                                                            double shape,
                                                            int[] vertices, double[][] points,
                                                            double[][] X) {
      GaussianFunction[] fs = new GaussianFunction[vertices.length];
      for (int i = 0; i < fs.length; ++i)
         fs[i] = new GaussianFunction(scale/vertices.length, sig[i], shape);
      return new AnchorPotential(fs, vertices, points, X);
   }
   public static AnchorPotential newGaussianAnchorPotential(double[] scale, double sig, 
                                                            double shape,
                                                            int[] vertices, double[][] points,
                                                            double[][] X) {
      GaussianFunction[] fs = new GaussianFunction[vertices.length];
      for (int i = 0; i < fs.length; ++i)
         fs[i] = new GaussianFunction(scale[i]/vertices.length, sig, shape);
      return new AnchorPotential(fs, vertices, points, X);
   }
   public static AnchorPotential newGaussianAnchorPotential(double scale, double sig, double shape,
                                                            int[] vertices, double[][] points, 
                                                            double[][] X) {
      return new AnchorPotential(new GaussianFunction(scale/vertices.length, sig, shape),
                                 vertices, points, X);
   }

   public static AnchorPotential newGaussianAnchorPotential(double[] scale, double[] sig,
                                                            int[] vertices, double[][] points, 
                                                            double[][] X) {
      GaussianFunction[] fs = new GaussianFunction[vertices.length];
      for (int i = 0; i < fs.length; ++i)
         fs[i] = new GaussianFunction(scale[i]/vertices.length, sig[i], 2.0);
      return new AnchorPotential(fs, vertices, points, X);
   }
   public static AnchorPotential newGaussianAnchorPotential(double[] scale, double sig,
                                                            int[] vertices, double[][] points, 
                                                            double[][] X) {
      GaussianFunction[] fs = new GaussianFunction[vertices.length];
      for (int i = 0; i < fs.length; ++i)
         fs[i] = new GaussianFunction(scale[i]/vertices.length, sig, 2.0);
      return new AnchorPotential(fs, vertices, points, X);
   }
   public static AnchorPotential newGaussianAnchorPotential(double scale, double[] sig,
                                                            int[] vertices, double[][] points, 
                                                            double[][] X) {
      GaussianFunction[] fs = new GaussianFunction[vertices.length];
      for (int i = 0; i < fs.length; ++i)
         fs[i] = new GaussianFunction(scale/vertices.length, sig[i], 2.0);
      return new AnchorPotential(fs, vertices, points, X);
   }
   public static AnchorPotential newGaussianAnchorPotential(double scale, double sig,
                                                            int[] vertices, double[][] points, 
                                                            double[][] X) {
      return new AnchorPotential(new GaussianFunction(scale/vertices.length, sig, 2.0),
                                 vertices, points, X);
   }

   public static AnchorPotential newGaussianAnchorPotential(double[] sig,
                                                            int[] vertices, double[][] points, 
                                                            double[][] X) {
      GaussianFunction[] fs = new GaussianFunction[vertices.length];
      for (int i = 0; i < fs.length; ++i)
         fs[i] = new GaussianFunction(1.0/vertices.length, sig[i], 2.0);
      return new AnchorPotential(fs, vertices, points, X);
   }

   public static AnchorPotential newGaussianAnchorPotential(double sig,
                                                            int[] vertices, double[][] points, 
                                                            double[][] X) {
      return new AnchorPotential(new GaussianFunction(1.0/vertices.length, sig, 2.0),
                                 vertices, points, X);
   }

   public static AnchorPotential newGaussianAnchorPotential(int[] vertices, double[][] points, 
                                                            double[][] X) {
      return new AnchorPotential(new GaussianFunction(1.0/vertices.length, 1.0, 2.0),
                                 vertices, points, X);
   }

   /** newHarmonicPerimeterPotential(scale, shape, faces, X0) yields a potential function that
    *  includes a harmonic perimeter potential function with the given scale and shape (see
    *  HarmonicFunction). This perimeter prevents vertices around the perimeter of the mesh from
    *  moving very much. The scale passed to the HarmonicFunction is scale / q where q is the number
    *  of vertices in the perimeter.
    */
   public static AnchorPotential newHarmonicPerimeterPotential(double scale, double shape,
                                                               int[][] faces, double[][] X) {
      int[] perim = Util.perimeter(faces);
      double[][] perimX0 = new double[X.length][perim.length];
      for (int j = 0; j < X.length; ++j)
         for (int i = 0; i < perim.length; ++i)
            perimX0[j][i] = X[j][perim[i]];
      return new AnchorPotential(new HarmonicFunction(scale/perim.length, shape), 
                                 perim, perimX0, X);
   }
   public static AnchorPotential newHarmonicPerimeterPotential(double scale,
                                                               int[][] faces, double[][] X) {
      return newHarmonicPerimeterPotential(scale, 2.0, faces, X);
   }
   public static AnchorPotential newHarmonicPerimeterPotential(int[][] faces, double[][] X) {
      return newHarmonicPerimeterPotential(1.0, 2.0, faces, X);
   }

   /** newStandardMeshPotential(faces, X0) yields a SumPotential that includes:
    *    (a) an angle potential for the angles in the mesh [infinite-well: scale = 1/p, min = 0,
    *        max = pi/2, shape = 0.5]
    *    (b) an edge potential for all edges in the given set of faces [harmonic: scale = 1/m, 
    *        shape = 2], 
    */
   public static PotentialSum newStandardMeshPotential(int[][] faces, double[][] X) {
      return new PotentialSum(newAngleWellPotential(faces, X),
                              newHarmonicEdgePotential(10.0, 2.0, Util.facesToEdges(faces), X),
                              newHarmonicPerimeterPotential(1.0, 2.0, faces, X));
   }

   /** newSum(fields) yields a new potential field object that is the sum of the given lsit of
    *  fields.
    */
   public static PotentialSum newSum() {
      return new PotentialSum();
   }

}
