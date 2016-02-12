////////////////////////////////////////////////////////////////////////////////////////////////////
// AInPlaceCalculator.java
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

import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;

/** AInPlaceCalculator is an abstract class that is used by low-level highly-optimized routines
 *  in the nben.registration minimization engine to ensure fast parallel access to gradient 
 *  calculations with minimal reallocation or copying of memory required.
 *  Every potential field must be able to produce an in-place-calculator, which is basically a
 *  workspace for the gradient calculation. All calculators must place the gradient into the public
 *  double array gradient, which is a member of AInPlaceCalculation. They should return Infinity or
 *  negative Infinity to indicate a singularity and they should return NaN to indicate an error of
 *  some kind.
 *
 *  @author Noah C. Benson
 */
public abstract class AInPlaceCalculator {

   /** subset is an array of vertex id's that should be considered by this calculator.
    */
   public int[] subset;
   /** X is a (dims x vertices)-sized matrix from which the calculator obtains the coordinates for
    *  the potential and gradient calculations each tiem calculate() is called. This matrix should
    *  be set or overwritten each time calculate() is called.
    */
   public double[][] X;
   /** potential is the value of the potential that is set by the calculate() function, according
    *  to the vertex configuration of the matrix X at the time of calling.
    */
   public double potential;
   /** gradient is a (dims x vertices)-sized matrix into which to write the gradient vectors for
    *  each vertex as column vectors whenever the potential() function is called. 
    */
   public double[][] gradient;
   /** gradientNorms is a (vertices)-sized vector into gradient norms are written by the gradNorms
    *  function.
    */
   public double[] gradientNorms;
   /** maxGradientNorm is the maximum gradient norm of all the vertices; it is filled in by the
    *  gradNorms function.
    */
   public double maxGradientNorm;
   /** gradientLength is the total length of the gradient vector
    */
   public double gradientLength;

   public AInPlaceCalculator(int[] ss, double[][] x, double[][] g, double[] gnorm) {
      if (ss == null) {
         ss = new int[x[0].length];
         for (int i = 0; i < ss.length; ++i) ss[i] = i;
      }
      if (g == null) g = new double[x.length][x[0].length];
      subset = ss;
      X = x;
      gradient = g;
      gradientNorms = (gnorm == null? new double[x[0].length] : gnorm);
      potential = Double.NaN;
   }
   public AInPlaceCalculator(int[] ss, double[][] x, double[][] g) {
      this (ss, x, g, null);
   }
   public AInPlaceCalculator(int[] ss, double[][] x) {
      this (ss, x, null, null);
   }

   /** The calculate function updates the gradient array and the potential value based on the 
    *  current contents of the matrix X. Note that this adds the gradient to the gradient vector;
    *  it does not overwrite the contents. Throws an exception on error.
    *
    *  @param exc the ThreadPoolExecutor service to use for multi-threading
    *  @param nworkers the number of workers to use in the calculation
    */
   public abstract void calculate(int nworkers, ExecutorService exc) throws Exception;
   /** Calls calculate(n, e) for a temporary executor e and a number of workers n detected according
    *  to the number of processors on the system.
    */
   public void calculate() throws Exception {
      calculate(Util.workers(), Util.pool());
   }

   private double[][] copyAr(double[][] a) {
      double[][] b = new double[a.length][];
      for (int i = 0; i < a.length; ++i) {
         b[i] = new double[a[i].length];
         System.arraycopy(a[i], 0, b[i], 0, a[0].length);
      }
      return b;         
   }
   private double[][] copyToAr(double[][] a, double[][] b) {
      for (int i = 0; i < a.length; ++i) {
         System.arraycopy(a[i], 0, b[i], 0, a[0].length);
      }
      return b;         
   }

   /** The calculateAt function calculates the gradient and potential at a particular position and
    *  yields the result without modifying the rest of the calculator; note that this funciton is
    *  note remotely thread-safe (as with most of the AInPlaceCalculator functions).
    *  Note that this function clears the calculator before running the calculations.
    *  If the argument grad is given and is not null, then the gradient is placed in that array;
    *  otherwise it is placed in the calculator's gradient array.
    */
   public double calculateAt(double[][] atX, double[][] grad) throws Exception {
      if (grad == null)
         grad = gradient;
      // swap out the potential, X, and the gradient..
      double[][] Xtmp = copyAr(X);
      double[][] Gtmp = copyAr(gradient);
      double Ptmp = potential;
      copyToAr(atX, X);
      copyToAr(grad, gradient);
      // run the calculation...
      clear();
      calculate();
      // now replace things
      if (atX != X) copyToAr(Xtmp, X);
      if (gradient != grad) copyToAr(Gtmp, gradient);
      double result = potential;
      potential = Ptmp;
      // return the result...
      return result;
   }
   public double calculateAt(double[][] atX) throws Exception {return calculateAt(atX, gradient);}

   /** runThreads(rs, exc) runs the given set of workers over the threads of the given executor 
    *  service and returns true on success or throws an exception on error. If the executor is
    *  null, then one is chosen automatically.
    */
   protected final boolean runThreads(Runnable[] rs, ExecutorService exc) throws Exception {
      if (exc == null) exc = Util.pool();
      if (rs == null || rs.length == 0) {
         return false;
      } else if (rs.length == 1) {
         rs[0].run();
      } else {
         int i;
         Future[] fut = new Future[rs.length];
         for (i = 0; i < rs.length; ++i) {
            fut[i] = exc.submit(rs[i]);
         }
         for (i = 0; i < rs.length; ++i) {
            fut[i].get();
         }
      }
      return true;
   }
   protected final boolean runThreads(Runnable[] rs) throws Exception {
      return runThreads(rs, null);
   }

   private class ClearWorker implements Runnable {
      int id;
      int workers;
      public ClearWorker(int i, int ws) {id = i; workers = ws;}
      public void run() {
         for (int j = 0; j < gradient.length; ++j) {
            for (int i = id; i < subset.length; i += workers)
               gradient[j][subset[i]] = 0;
         }
      }
   }

   /** clear() resets the gradient at the given subset to be 0 and sets the potential value to 0.
    */
   public void clear(int nworkers, ExecutorService exc) {
      if (nworkers == 0) nworkers = Util.workers();
      if (exc == null) exc = Util.pool();
      potential = 0;
      ClearWorker[] cw = new ClearWorker[nworkers];
      for (int i = 0; i < nworkers; ++i)
         cw[i] = new ClearWorker(i, nworkers);
      try {
         runThreads(cw);
      } catch (Exception e) {
         for (int j = 0; j < gradient.length; ++j) {
            for (int i = 0; i < subset.length; ++i)
               gradient[j][subset[i]] = 0;
         }
      }
   }
   public void clear(int nworkers) {clear(nworkers, null);}
   public void clear() {clear(0, null);}

   private class GradWorker implements Runnable {
      int id;
      int workers;
      double maxgrad;
      double gradlen;
      public GradWorker(int i, int ws) {id = i; workers = ws;}
      public void run() {
         double tmp;
         int i, j, u;
         if (subset == null) {
            for (i = id; i < gradientNorms.length; i += workers)
               gradientNorms[i] = 0;
            for (j = 0; j < gradient.length; ++j) {
               for (i = id; i < gradientNorms.length; i += workers) {
                  tmp = gradient[j][i];
                  gradientNorms[i] += tmp*tmp;
                  gradlen += tmp*tmp;
               }
            }
            for (i = id; i < gradientNorms.length; i += workers) {
               gradientNorms[i] = Math.sqrt(gradientNorms[i]);
               if (i == id || gradientNorms[i] > maxGradientNorm) 
                  maxgrad = gradientNorms[i];
            }
         } else {
            for (i = id; i < subset.length; i += workers)
               gradientNorms[subset[i]] = 0;
            for (j = 0; j < gradient.length; ++j) {
               for (i = id; i < subset.length; i += workers) {
                  u = subset[i];
                  tmp = gradient[j][u];
                  gradientNorms[u] += tmp*tmp;
                  gradlen += tmp*tmp;
               }
            }
            for (i = id; i < subset.length; i += workers) {
               u = subset[i];
               gradientNorms[u] = Math.sqrt(gradientNorms[u]);
               if (i == id || gradientNorms[u] > maxGradientNorm) 
                  maxgrad = gradientNorms[u];
            } 
         }
      }
   }
   /** gradNorms() fills in the gradientNorms member of the given calculator and returns the maximum
    *  gradient across all vertices in the calculator's subset.
    */
   public double gradNorms(int nworkers, ExecutorService exc) throws Exception {
      maxGradientNorm = 0;
      GradWorker[] gw = new GradWorker[nworkers];
      for (int i = 0; i < nworkers; ++i)
         gw[i] = new GradWorker(i, nworkers);
      runThreads(gw);
      maxGradientNorm = gw[0].maxgrad;
      gradientLength = 0;
      for (int i = 1; i < nworkers; ++i) {
         gradientLength += gw[i].gradlen;
         if (gw[i].maxgrad > maxGradientNorm)
            maxGradientNorm = gw[i].maxgrad;
      }
      gradientLength = Math.sqrt(gradientLength);
      return maxGradientNorm;
   }
   public double gradNorms(int nw) throws Exception {return gradNorms(nw, Util.pool());}
   public double gradNorms() throws Exception {return gradNorms(Util.workers(), Util.pool());}
}

