////////////////////////////////////////////////////////////////////////////////////////////////////
// AnglePotential.java
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

import nben.registration.Util;
import nben.registration.APotentialField;
import nben.registration.IDifferentiatedFunction;

import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;

/** The AnglePotential class defines the code computations of potential fields based on the
 *  interactions between triples of neighboring vertices in the mesh. To construct an 
 *  AnglePotential, one must explicitly or implicitly choose a differentiated function f; the
 *  potential calculated by the edge potential class will then be the value f(a, a0) where a and a0
 *  are the angles and the reference angles of the mesh, respectively.
 *
 *  @author Noah C. Benson
 */
public class AnglePotential extends APotentialField {
   // the form of the potential function
   public final IDifferentiatedFunction form;
   // The angles we manage; angle i consists of the vectors from angles[0][i] to angles[1][i] and 
   // angles[2][i]
   public final int[][] angles;
   // angleIndex[i] is a list of the edge indices for vertex i; i.e., if q = angleIndex[i][j],
   // either angles[0][q] or angles[1][q] or angles[2][q] will be equal to i
   public final int[][] angleIndex;
   // T0 is an array of the original angles in the mesh
   public final double[] T0;
   // these are for convenience/use with the edge worker classes below
   private final int[] allVertices;
   private final int[] allAngles;

   public static final double[] cross(double[] a, double[] b) {
      double[] c = new double[3];
      c[0] = a[1]*b[2] - a[2]*b[1];
      c[1] = a[2]*b[0] - a[0]*b[2];
      c[2] = a[0]*b[1] - a[1]*b[0];
      return c;
   }

   /** anglePotential.calculateAngle(id, X, G) calculates the angle labeled by the given id,
    *  looking up the X coordinaets in the given coordinate matrix (size dims x vertices) X.
    *  If the final argument G is non-null, it places the gradient value in the appropriate
    *  entries of the (3 x dims x n)-sized matrix G where n is the number of angles.
    *
    *  @returns the value of angle id at the position X
    */
   public final double calculateAngle(int id, double[][] X, double[][][] G) {
      double abx, aby, acx, acy, theta;
      int i, dims = X.length;
      double[] nAB = new double[dims];
      double[] nAC = new double[dims];
      double dAC, dAB;
      if (dims == 2) {
         abx = X[0][angles[1][id]] - X[0][angles[0][id]];
         aby = X[1][angles[1][id]] - X[1][angles[0][id]];
         acx = X[0][angles[2][id]] - X[0][angles[0][id]];
         acy = X[1][angles[2][id]] - X[1][angles[0][id]];
         theta = Math.atan2(acy, acx) - Math.atan2(aby, abx);
         dAB = Math.sqrt(abx*abx + aby*aby);
         nAB[0] = abx/dAB;
         nAB[1] = aby/dAB;
         dAC = Math.sqrt(acx*acx + acy*acy);
         nAC[0] = acx/dAC;
         nAC[1] = acy/dAC;
      } else {
         // we have to convert to a flattened 2d form...
         int A = angles[0][id];
         int B = angles[1][id];
         int C = angles[2][id];
         double[] axisZ = new double[3];
         double tmp = Math.sqrt(X[0][A]*X[0][A] + X[1][A]*X[1][A] + X[2][A]*X[2][A]);
         axisZ[0] = X[0][A]/tmp;
         axisZ[1] = X[1][A]/tmp;
         axisZ[2] = X[2][A]/tmp;
         nAB[0] = X[0][B] - X[0][A];
         nAB[1] = X[1][B] - X[1][A];
         nAB[2] = X[2][B] - X[2][A];
         dAB = Math.sqrt(nAB[0]*nAB[0] + nAB[1]*nAB[1] + nAB[2]*nAB[2]);
         nAB[0] /= dAB;
         nAB[1] /= dAB;
         nAB[2] /= dAB;
         double[] axisY = cross(axisZ, nAB);
         tmp = Math.sqrt(axisY[0]*axisY[0] + axisY[1]*axisY[1] + axisY[2]*axisY[2]);
         axisY[0] /= tmp;
         axisY[1] /= tmp;
         axisY[2] /= tmp;
         double[] axisX = cross(axisY, axisZ);
         // now we can get the 2D coordinates out...
         nAC[0] = X[0][C] - X[0][A];
         nAC[1] = X[1][C] - X[1][A];
         nAC[2] = X[2][C] - X[2][A];
         acx = nAC[0]*axisX[0] + nAC[1]*axisX[1] + nAC[2]*axisX[2];
         acy = nAC[0]*axisY[0] + nAC[1]*axisY[1] + nAC[2]*axisY[2];
         dAC = Math.sqrt(nAC[0]*nAC[0] + nAC[1]*nAC[1] + nAC[2]*nAC[2]);
         nAC[0] /= dAC;
         nAC[1] /= dAC;
         nAC[2] /= dAC;
         // ab is lined up with the x-axis
         aby = 0.0;
         abx = dAB * (nAB[0]*axisX[0] + nAB[1]*axisX[1] + nAB[2]*axisX[2]);
         // and theta...
         theta = Math.atan2(acy, acx) - Math.atan2(aby, abx);
      }
      // make sure the vectors didn't cross the -x axis
      if (theta < -Math.PI) theta += 2.0*Math.PI;
      else if (theta > Math.PI) theta -= 2.0*Math.PI;
      // set G to the gradient of theta in terms of X
      if (G != null) {
         double cos = Math.cos(theta);
         double sin = Math.sqrt(1.0 - cos*cos);
         double[][] g0 = G[0];
         double[][] g1 = G[1];
         double[][] g2 = G[2];
         for (i = 0; i < dims; ++i) {
            // start with corner 1:
            g1[i][id] = (cos * nAB[i] - nAC[i]) / (sin*dAB);
            // then corner 2:
            g2[i][id] = (cos * nAC[i] - nAB[i]) / (sin*dAC);
            // finally corner 0:
            g0[i][id] = -(g1[i][id] + g2[i][id]);
         }
      }
      return theta;
   }
   

   /** Constructs an AnglePotential object.
    *  
    *  @param f an IDifferentiatedFunction object indicating the form of the potential landscape
    *  @param faces the (3 x n) list of triangles; must be 0-indexed and must have triangles listed
    *               all in the same ordering (clockwise or counter-clockwise)
    *  @param X0 a (dims x n) array of the starting coordinates of the vertices
    */
   public AnglePotential(IDifferentiatedFunction f, int[][] angles, double[][] X0) {
      this.form = f;
      int n = angles[0].length; // number of angles

      // make the angles matrix... we dup this just in case
      this.angles = new int[3][n];
      for (int i = 0; i < 3; ++i) System.arraycopy(angles[i], 0, this.angles[i], 0, n);
      // build the index...
      this.angleIndex = Util.buildSimplexIndex(X0[0].length, this.angles);
      // save the original angles
      int q;
      T0 = new double[n];
      for (int i = 0; i < n; ++i) {
         T0[i] = calculateAngle(i, X0, null);
         if (T0[i] < 0) {
            // we want all angles to be positive...
            q = angles[1][i];
            angles[1][i] = angles[2][i];
            angles[2][i] = q;
            T0[i] = -T0[i];
         }
      }
      // fill these in for convenience
      allAngles = new int[n];
      for (int i = 0; i < n; ++i) allAngles[i] = i;
      allVertices = new int[X0[0].length];
      for (int i = 0; i < allVertices.length; ++i) allVertices[i] = i;
   }
   
   // Here we have the code/subclasses that handle the workers
   private final class AngleCalculation extends AInPlaceCalculator {
      // the scratch space we need...
      private int[] asubset; // the subset of angles we operate over
      private double[] adat; // angle values; n where n = num angles
      private double[] gdat; // angle derivative values; n where n = num angles
      // the normalized gradient directions; 3 x d x n; first dim is the vertex number of the angle,
      // (0 is the corner, 1 is the next in order, and 2 is the next), second dim d is the 
      // dimension, and last dim is the angle number.
      private double[][][] dirs;

      public double[] getADat() {return adat;}
      public double[] getGDat() {return gdat;}
      public double[][][] getDirs() {return dirs;}

      // constructor
      public AngleCalculation(int[] ss, double[][] X0, double[][] G) {
         super(ss, X0, G);
         this.dirs = new double[3][X0.length][T0.length];
         this.adat = new double[T0.length];
         this.gdat = new double[T0.length];
         // must build the edge subset from the vertex subset
         if (ss == null) {
            this.asubset = allAngles;
         } else {
            this.asubset = Util.subsampleIndex(subset, angleIndex);
         }
      }
      // this class calculates the first stage (fill in AB and edat) -- operate over edges
      public final class AngleWorker1 implements Runnable {
         int id;
         int workers;
         public AngleWorker1(int myid, int nworkers) {id = myid; workers = nworkers;}
         public void run() {
            int n = asubset.length, i, j, a;
            double tmp, cos, sin;
            double[] x;
            for (i = id; i < n; i += workers) {
               a = asubset[i];
               // calculate the angle and gradient...
               adat[a] = calculateAngle(a, X, dirs);
               gdat[a] = form.dy(adat[a], T0[a]);
               adat[a] = form.y(adat[a], T0[a]);
            }
         }
      }
      // this class calculates the second stage (fill in G and PE) -- operate over vertices
      public final class AngleWorker2 implements Runnable {
         int id;
         int workers;
         double workerPE;
         public AngleWorker2(int myid, int nworkers) {id = myid; workers = nworkers; workerPE = 0;}
         public void run() {
            int n = subset.length, i, j, k, u, a, d0;
            int[] idx;
            for (i = id; i < n; i += workers) {
               u = subset[i];
               idx = angleIndex[u];
               for (j = 0; j < idx.length; ++j) {
                  a = idx[j];
                  workerPE += adat[a] / 3.0;
                  d0 = (angles[0][a] == u? 0 : (angles[1][a] == u? 1 : 2));
                  for (k = 0; k < gradient.length; ++k) 
                     gradient[k][u] += gdat[a] * dirs[d0][k][a];
               }
            }
         }
      }

      // Actually run the calculations...
      public void calculate(int nworkers, ExecutorService exc) throws Exception {
         int i;
         AngleWorker1[] aw1 = new AngleWorker1[nworkers];
         AngleWorker2[] aw2 = new AngleWorker2[nworkers];
         for (i = 0; i < nworkers; ++i) {
            aw1[i] = new AngleWorker1(i, nworkers);
            aw2[i] = new AngleWorker2(i, nworkers);
         }
         // submit the per-edge calculations now...
         potential = 0.0;
         if (runThreads(aw1, exc) && runThreads(aw2, exc)) {
            for (i = 0; i < aw2.length; ++i)
               potential += aw2[i].workerPE;
         } else {
            potential = Double.NaN;
         }
      }
   }

   public final AngleCalculation potentialCalculator(int[] subset, double[][] X, double[][] G) {
      return new AngleCalculation(subset, X, G);
   }

}