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

   public AInPlaceCalculator(int[] ss, double[][] x, double[][] g) {
      if (ss == null) {
         ss = new int[x[0].length];
         for (int i = 0; i < ss.length; ++i) ss[i] = i;
      }
      if (g == null) g = new double[x.length][x[0].length];
      subset = ss;
      X = x;
      gradient = g;
      potential = Double.NaN;
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
}

