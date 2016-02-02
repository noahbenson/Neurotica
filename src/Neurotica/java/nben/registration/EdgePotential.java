////////////////////////////////////////////////////////////////////////////////////////////////////
// EdgePotential.java
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

import nben.registration.PotentialFields;
import nben.registration.APotentialField;
import nben.registration.IDifferentiatedFunction;

import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;

/** The EdgePotential class defines the code computations of potential fields based on the
 *  interactions between neighboring vertices in the mesh. To construct an EdgePotential, one must
 *  explicitly or implicitly choose a differentiated function f; the potential calculated by the
 *  edge potential class will then be the value f(l, l0) where l and l0 are the edge-length and the
 *  reference edge-length respectively.
 *
 *  @author Noah C. Benson
 */
class EdgePotential extends APotentialField {
   // the form of the potential function
   public final IDifferentiatedFunction form;
   // The edges we manage
   public final int[][] edges;
   // edgeIndex[i] is a list of the edge indices for vertex i; i.e., if q = m_edgeIndex[i][j],
   // either edges[0][q] or edges[1][q] will be equal to i
   public final int[][] edgeIndex;
   // m_D0 is an array of the original distances between edges in the mesh
   public final double[] D0;
   // these are for convenience/use with the edge worker classes below
   private final int[] allVertices;
   private final int[] allEdges;

   public int[][] getEdges() {return edges;}
   public int[] getEdgeIndex(int i) {return edgeIndex[i];}
   public double[] getD0() {return D0;}

   /** Constructs an EdgePotential object.
    *  
    *  @param f an IDfferentiatedFunction that specifies the shape of the potential landscape
    *  @param edges the (2 x n) list of edges; must be 0-indexed
    *  @param X0 a (dims x n) array of the starting coordinates of the vertices
    */
   public EdgePotential(IDifferentiatedFunction f, int[][] edges, double[][] X0) {
      this.form = f;
      int n = edges[0].length;
      this.edges = new int[2][n];
      this.edgeIndex = new int[X0[0].length][];
      // copy the array, build the index...
      int[] u, v;
      for (int j = 0; j < 2; ++j) {
         for (int i = 0; i < n; ++i) {
            this.edges[j][i] = edges[j][i];
            u = edgeIndex[edges[j][i]];
            if (u == null) {
               u = new int[1];
               u[0] = i;
            } else {
               v = new int[u.length + 1];
               System.arraycopy(u, 0, v, 0, u.length);
               v[u.length] = i;
               u = v;
            }
            edgeIndex[edges[j][i]] = u;
         }
      }
      // save the original distances
      D0 = new double[n];
      for (int j = 0; j < X0.length; ++j) {
         double[] x0 = X0[j];
         double tmp;
         for (int i = 0; i < n; ++i) {
            tmp = x0[edges[0][i]] - x0[edges[1][i]];
            D0[i] += tmp*tmp;
         }
      }
      for (int i = 0; i < n; ++i) D0[i] = Math.sqrt(D0[i]);
      // fill these in for convenience
      allEdges = new int[n];
      for (int i = 0; i < n; ++i) allEdges[i] = i;
      allVertices = new int[X0[0].length];
      for (int i = 0; i < allVertices.length; ++i) allVertices[i] = i;
   }

   // Here we have the code/subclasses that handle the workers
   private final class EdgeCalculation extends AInPlaceCalculator {
      // the scratch space we need...
      private double[] edat;
      private double[] gdat;
      private double[][] AB; // the normalized vectors from a to b (dims x edges)
      private int[] esubset; // the subset of edges we operate over

      // constructor
      public EdgeCalculation(int[] ss, double[][] X0, double[][] G) {
         super(ss, X0, G);
         this.AB = new double[X.length][D0.length];
         this.edat = new double[D0.length];
         this.gdat = new double[D0.length];
         // must build the edge subset from the vertex subset
         if (ss == null) {
            this.esubset = allEdges;
         } else {
            this.esubset = PotentialFields.subsampleIndex(subset, edgeIndex);
         }
      }
      // this class calculates the first stage (fill in AB and edat) -- operate over edges
      public final class EdgeWorker1 implements Runnable {
         int id;
         int workers;
         public EdgeWorker1(int myid, int nworkers) {id = myid; workers = nworkers;}
         public void run() {
            int n = esubset.length, i, j, e;
            double tmp;
            double[] x;
            for (j = 0; j < X.length; ++j) {
               x = X[j];
               for (i = id; i < n; i += workers) {
                  e = esubset[i];
                  tmp = x[edges[1][e]] - x[edges[0][e]];
                  AB[j][e] = tmp;
                  edat[i] += tmp*tmp;
               }
            }
            for (i = id; i < n; i += workers) {
               e = esubset[i];
               edat[e] = Math.sqrt(edat[e]);
               for (j = 0; j < AB.length; ++j)
                  AB[j][e] /= edat[e];
               // now get potentials and gradient lengths
               gdat[e] = form.dy(edat[e], D0[e]);
               edat[e] = form.y(edat[e], D0[e]);
            }
         }
      }
      // this class calculates the second stage (fill in G and PE) -- operate over vertices
      public final class EdgeWorker2 implements Runnable {
         int id;
         int workers;
         double workerPE;
         public EdgeWorker2(int myid, int nworkers) {id = myid; workers = nworkers; workerPE = 0;}
         public void run() {
            int n = subset.length, i, j, k, u, e;
            int[] idx;
            double tmp;
            for (i = id; i < n; i += workers) {
               u = subset[i];
               idx = edgeIndex[u];
               for (j = 0; j < idx.length; ++j) {
                  e = idx[j];
                  workerPE += 0.5 * edat[e];
                  if (edges[0][e] == u)
                     for (k = 0; k < gradient.length; ++k) gradient[k][u] -= gdat[e] * AB[k][e];
                  else
                     for (k = 0; k < gradient.length; ++k) gradient[k][u] += gdat[e] * AB[k][e];
               }
            }
         }
      }

      // Actually run the calculations...
      public void calculate(int nworkers, ExecutorService exc) throws Exception {
         int i;
         EdgeWorker1[] ew1 = new EdgeWorker1[nworkers];
         EdgeWorker2[] ew2 = new EdgeWorker2[nworkers];
         for (i = 0; i < nworkers; ++i) {
            ew1[i] = new EdgeWorker1(i, nworkers);
            ew2[i] = new EdgeWorker2(i, nworkers);
         }
         // submit the per-edge calculations now...
         potential = 0.0;
         if (runThreads(ew1, exc) && runThreads(ew2, exc)) {
            for (i = 0; i < ew2.length; ++i)
               potential += ew2[i].workerPE;
         } else {
            potential = Double.NaN;
         }
      }
   }

   public final EdgeCalculation potentialCalculator(int[] subset, double[][] X, double[][] G) {
      return new EdgeCalculation(subset, X, G);
   }

}