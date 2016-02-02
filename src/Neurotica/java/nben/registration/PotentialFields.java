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

import nben.registration.APotentialField;
import nben.registration.IDifferentiatedFunction;
import nben.registration.HarmonicFunction;
import nben.registration.EdgePotential;

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
    *  number of edges.
    */
   public static EdgePotential newHarmonicEdgePotential(double scale, double shape,
                                                        int[][] E, double[][] X) {
      return new EdgePotential(new HarmonicFunction(scale / E[0].length, shape), E, X);
   }
   public static EdgePotential newHarmonicEdgePotential(double scale, int[][] E, double[][] X) {
      return new EdgePotential(new HarmonicFunction(scale / E[0].length), E, X);
   }
   public static EdgePotential newHarmonicEdgePotential(int[][] E, double[][] X) {
      return new EdgePotential(new HarmonicFunction(1.0 / E[0].length), E, X);
   }
   /** PotentialFields.newLJEdgePotential(s, q, E, X) yields an EdgePotential object with a 
    *  LennardJonesFunction form using scale parameter s/m and shape parameter q where where m is
    *  the number of edges.
    */
   public static EdgePotential newLJEdgePotential(double scale, double shape,
                                                  int[][] E, double[][] X) {
      return new EdgePotential(new LennardJonesFunction(scale / E[0].length, shape), E, X);
   }
   /** PotentialFields.newHarmonicAnglePotential(s, q, E, X) yields an EdgePotential object with a 
    *  LennardJonesFunction form using scale parameter s/m and shape parameter q where where m is
    *  the number of edges.
    */
   public static AnglePotential newHarmonicAnglePotential(double scale, double shape,
                                                          int[][] T, double[][] X) {
      return new AnglePotential(new HarmonicFunction(scale / (3.0 * T[0].length), shape), T, X);
   }
   /** PotentialFields.newLJAnglePotential(s, q, T, X) yields an AnglePotential object with a 
    *  LennardJonesFunction form using scale parameter s/m and shape parameter q where where m is
    *  the number of angles (triangles times 3).
    */
   public static AnglePotential newLJAnglePotential(double scale, double shape,
                                                    int[][] T, double[][] X) {
      return new AnglePotential(new LennardJonesFunction(scale / (3.0 * T[0].length), shape), T, X);
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
                                T, X);
   }
   public static AnglePotential newAngleWellPotential(double scale, 
                                                      double min, double max,
                                                      int[][] T, double[][] X) {
      return new AnglePotential(new InfiniteWellFunction(scale / (3.0 * T[0].length), 
                                                         min, max),
                                T, X);
   }
   public static AnglePotential newAngleWellPotential(double scale, 
                                                      int[][] T, double[][] X) {
      return new AnglePotential(new InfiniteWellFunction(scale / (3.0 * T[0].length)), T, X);
   }
   public static AnglePotential newAngleWellPotential(int[][] T, double[][] X) {
      return new AnglePotential(new InfiniteWellFunction(1.0 / (3.0 * T[0].length)), T, X);
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

   /** newStandardMeshPotential(faces, X0) yields a SumPotential that includes:
    *    (a) an edge potential for all edges in the given set of faces [harmonic: scale = 1/m, 
    *        shape = 2], 
    *    (b) an angle potential for the angles in the mesh [infinite-well: scale = 1/p, min = 0,
    *        max = pi/2, shape = 0.5]
    */
   public static PotentialSum newStandardMeshPotential(int[][] faces, double[][] X) {
      return new PotentialSum(newAngleWellPotential(faces, X),
                              newHarmonicEdgePotential(Util.facesToEdges(X[0].length, faces), X));
   }
}
