////////////////////////////////////////////////////////////////////////////////////////////////////
// AnchorPotential.java
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

/** The AnchorPotential class defines the code computations of potential fields based on the
 *  interactions between a vertex's position in the mesh and a fixed point in space. To construct an
 *  AnchorPotential, one must explicitly or implicitly choose a differentiated function f; the
 *  potential calculated by the anchor potential class will then be the value f(r, r0) where r and
 *  r0 are the distance from the vertex to the fixed point and the reference distance, respectively.
 *  Note that this class does not currently support a reference value (all references are 0) for the
 *  the IDifferentiatedFunction form of the potential; this may change in the future.
 *
 *  @author Noah C. Benson
 */
class AnchorPotential extends APotentialField {
   // the form of the potential function
   public final IDifferentiatedFunction[] forms;
   // The vertices we manage (length m)
   public final int[] vertices;
   // If vertexIndex[i] = {k1, k2,...}, then vertices[k1] == vertices[k2] == ... == i;
   // if vertex i does not appear in vertices, then vertexIndex[i] == null.
   public final int[][] vertexIndex;
   // The points to which the vertices are fixed (size dims x m) where m is vertices.length
   public final double[][] anchors;
   // R0 is an array of the original distances between the vertices[] and their fixed points
   // in the reference mesh (length = vertices.length)
   //public final double[] R0;
   
   // these are for convenience/use with the edge worker classes below
   private final int[] allVertices;
   private final int[] allAnchors;

   public int[] getAnchorVertices() {return vertices;}
   public int[] getAnchorVertexIndex(int i) {return vertexIndex[i];}
   public double[][] getAnchors() {return anchors;}
   //public double[] getR0() {return R0;}

   /** Constructs an AnchorPotential object.
    *  
    *  @param f an IDfferentiatedFunction that specifies the shape of the potential landscape
    *  @param ids an array of the indices of the vertices that are drawn toward an anchor point
    *  @param anchors a (dims x m) array of the anchor points for the given vertices.
    *  @param X0 a (dims x n) array of the starting coordinates of the vertices
    */
   public AnchorPotential(IDifferentiatedFunction[] f, 
                          int[] vertexIDs, double[][] anchorsX,
                          double[][] X0) {
      int n = X0[0].length;
      int m = vertexIDs.length;
      int dims = X0.length;

      this.forms = f;

      vertices = new int[m];
      System.arraycopy(vertexIDs, 0, vertices, 0, m);

      anchors = new double[X0.length][m];
      for (int i = 0; i < dims; ++i)
         System.arraycopy(anchorsX[i], 0, anchors[i], 0, m);

      vertexIndex = new int[n][];
      int u, count = 0; // we track count so we know how many vertices are involved...
      int[] ar, artmp;
      for (int i = 0; i < m; ++i) {
         u = vertices[i];
         ar = vertexIndex[u];
         if (ar == null) {
            artmp = new int[1];
            artmp[0] = i;
            ++count;
         } else {
            artmp = new int[ar.length + 1];
            System.arraycopy(ar, 0, artmp, 0, ar.length);
            artmp[ar.length] = i;
         }
         vertexIndex[u] = artmp;
      }

      // fill these in for convenience
      allAnchors = new int[m];
      for (int i = 0; i < m; ++i) allAnchors[i] = i;
      allVertices = new int[count];
      u = 0;
      for (int i = 0; i < n; ++i) {
         if (vertexIndex[i] != null)
            allVertices[u++] = i;
      }
   }
   private final static IDifferentiatedFunction[] fillDiffFns(IDifferentiatedFunction f, int n) {
      IDifferentiatedFunction[] fs = new IDifferentiatedFunction[n];
      for (int i = 0; i < n; ++i) fs[i] = f;
      return fs;
   }
   public AnchorPotential(IDifferentiatedFunction f,
                          int[] vertexIDs, double[][] anchorsX,
                          double[][] X0) {
      this(fillDiffFns(f, vertexIDs.length), vertexIDs, anchorsX, X0);
   }
   
   // Here we have the code/subclasses that handle the workers
   private final class AnchorCalculation extends AInPlaceCalculator {
      // the scratch space we need...
      private double[] rdat;
      private double[] gdat;
      private double[][] AB; // the normalized vectors from anchor to point (dims x vertices)
      private int[] asubset; // the subset of anchors we operate over
      private final int m_dims;

      public double[] getRDat() {return rdat;}
      public double[] getGDat() {return gdat;}
      public double[][] getAB() {return AB;}
      public int[] getss() {return subset;}

      // constructor
      public AnchorCalculation(int[] ss, double[][] X0, double[][] G, double[] Gn) {
         super((ss == null? allVertices : ss), X0, G, Gn);
         this.AB = new double[X.length][vertices.length];
         this.rdat = new double[vertices.length];
         this.gdat = new double[vertices.length];
         m_dims = X0.length;
         // must build the edge subset from the vertex subset
         if (ss == null)
            this.asubset = allAnchors;
         else
            this.asubset = Util.subsampleIndex(subset, vertexIndex);
      }
      // this class calculates the first stage (fill in AB and rdat) -- operate over anchors
      public final class AnchorWorker1 implements Runnable {
         int id;
         int workers;
         public AnchorWorker1(int myid, int nworkers) {id = myid; workers = nworkers;}
         public void run() {
            int n = asubset.length, i, j, a, u;
            double tmp;
            double[] x;
            // fill in AB and rdat with direction vectors and distances...
            for (i = id; i < n; i += workers)
               rdat[asubset[i]] = 0;
            for (j = 0; j < anchors.length; ++j) {
               for (i = id; i < n; i += workers) {
                  a = asubset[i];
                  u = vertices[a];
                  tmp = X[j][u] - anchors[j][a];
                  rdat[a] += tmp*tmp;
                  AB[j][a] = tmp;
               }
            }
            for (i = id; i < n; i += workers) {
               a = asubset[i];
               tmp = Math.sqrt(rdat[a]);
               // normalize vectors...
               if (!Util.zeroish(tmp))
                  for (j = 0; j < m_dims; ++j) AB[j][a] /= tmp;
               // get potential and gradient data...
               gdat[a] = forms[a].dy(tmp, 0.0);
               rdat[a] = forms[a].y(tmp, 0.0);
            }
         }
      }
      // this class calculates the second stage (fill in G and PE) -- operate over vertices
      public final class AnchorWorker2 implements Runnable {
         int id;
         int workers;
         double workerPE;
         public AnchorWorker2(int myid, int nworkers) {id = myid; workers = nworkers; workerPE = 0;}
         public void run() {
            int n = subset.length, i, j, k, u, a;
            int[] idx;
            double tmp;
            for (i = id; i < n; i += workers) {
               u = subset[i];
               idx = vertexIndex[u];
               if (idx != null) {
                  for (j = 0; j < idx.length; ++j) {
                     a = idx[j];
                     workerPE += rdat[a];
                     for (k = 0; k < gradient.length; ++k) gradient[k][u] += gdat[a] * AB[k][a];
                  }
               }
            }
         }
      }

      // Actually run the calculations...
      public void calculate(int nworkers, ExecutorService exc) throws Exception {
         int i;
         AnchorWorker1[] aw1 = new AnchorWorker1[nworkers];
         AnchorWorker2[] aw2 = new AnchorWorker2[nworkers];
         for (i = 0; i < nworkers; ++i) {
            aw1[i] = new AnchorWorker1(i, nworkers);
            aw2[i] = new AnchorWorker2(i, nworkers);
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

   public final AnchorCalculation potentialCalculator(int[] subset, double[][] X, 
                                                      double[][] G, double[] Gn) {
      return new AnchorCalculation(subset, X, G, Gn);
   }
}