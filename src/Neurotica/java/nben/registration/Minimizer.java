////////////////////////////////////////////////////////////////////////////////////////////////////
// Minimizer.java
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
import java.util.concurrent.Future;

import java.util.Arrays;

/** Minimizer is a class that handles the registration minimization. Minimizer, by default, uses the
 *  Util class's executor and worker suggestions.
 */
public class Minimizer {

   ////////////////////////////////////////////////////////////////////////////////
   // Private Classes
   private class StepWorker implements Runnable {
      public final int id;
      public final int workers;
      public int[] subset;
      public boolean direction;
      public StepWorker(int i, int nwork) {
         id = i;
         workers = nwork;
         subset = null;
      }
      public void run() {
         if (direction) { // direction = true: step forward
            if (subset == null) {
               int j;
               for (int i = id; i < m_gradientNorms.length; i += workers) {
                  for (j = 0; j < m_gradient.length; ++j) {
                     m_X0[j][i] = m_X[j][i];
                     m_gradient0[j][i] = m_gradient[j][i];
                     m_X[j][i] -= m_gradient[j][i]*m_stepSize;
                  }
               }
            } else {
               int j, u;
               for (int i = id; i < subset.length; i += workers) {
                  u = subset[i];
                  for (j = 0; j < m_gradient.length; ++j) {
                     m_X0[j][u] = m_X[j][u];
                     m_gradient0[j][u] = m_gradient[j][u];
                     m_X[j][u] -= m_gradient[j][u]*m_stepSize;
                  }
               }
            }
         } else { // direction = false: step back (unroll a step)
            if (subset == null) {
               int j;
               for (int i = id; i < m_gradientNorms.length; i += workers) {
                  for (j = 0; j < m_gradient.length; ++j) {
                     m_X[j][i] = m_X0[j][i];
                     m_gradient[j][i] = m_gradient0[j][i];
                  }
               }
            } else {
               int j, u;
               for (int i = id; i < subset.length; i += workers) {
                  u = subset[i];
                  for (j = 0; j < m_gradient.length; ++j) {
                     m_X[j][u] = m_X0[j][u];
                     m_gradient[j][u] = m_gradient0[j][u];
                  }
               }
            }
         }
      }
   }
   /** The Minimizer.Report class stores the data from a minimization trajectory and is returned
    *  by the step function.
    */
   public class Report {
      public double[] stepSizes;
      public double[] stepLengths;
      public double[] steepestVertexGradientNorms;
      public double[] potentialChanges;
      public double initialPotential;
      public double finalPotential;
      public int steps;

      private Report(double pe0) {
         steps = 0;
         stepSizes = new double[16];
         stepLengths = new double[16];
         steepestVertexGradientNorms = new double[16];
         potentialChanges = new double[16];
         initialPotential = pe0;
         finalPotential = 0;
         steps = 0;
      }
      private double[] extend(double[] a) {
         if (a == null || a.length == 0)
            return new double[16];
         double[] tmp = new double[2*a.length];
         System.arraycopy(a, 0, tmp, 0, a.length);
         return tmp;
      }
      private double[] trim(double[] a, int mx) {
         double[] b = new double[mx];
         System.arraycopy(a, 0, b, 0, mx);
         return b;
      }
      private void freeze(double pe) {
         finalPotential = pe;
         stepSizes = trim(stepSizes, steps);
         stepLengths = trim(stepLengths, steps);
         steepestVertexGradientNorms = trim(steepestVertexGradientNorms, steps);
         potentialChanges = trim(potentialChanges, steps);
      }
      private void push(double dt, double dx, double maxnorm, double dpe) {
         if (steps == stepSizes.length) {
            stepSizes = extend(stepSizes);
            stepLengths = extend(stepLengths);
            steepestVertexGradientNorms = extend(steepestVertexGradientNorms);
            potentialChanges = extend(potentialChanges);
         }
         stepSizes[steps] = dt;
         stepLengths[steps] = dx;
         steepestVertexGradientNorms[steps] = maxnorm;
         potentialChanges[steps] = dpe;
         steps++;
      }
   }

   ////////////////////////////////////////////////////////////////////////////////
   // Private Data
   private APotentialField    m_field;
   private AInPlaceCalculator m_calculator;
   private double[][]         m_X0;
   private double[][]         m_X;
   private double             m_potential;
   private double[][]         m_gradient;
   private double[][]         m_gradient0;
   // workspace used by this object
   private double[]           m_gradientNorms;
   private double             m_gradientTotalNorm;
   private double             m_stepSize;
   private StepWorker[]       m_stepWorkers;
   // public data meant to be examined by outer callers
   public Report report;
   

   ////////////////////////////////////////////////////////////////////////////////
   // Accessors
   public APotentialField getPotentialField() {return m_field;}
   public AInPlaceCalculator getPotentialCalculator() {return m_calculator;}
   public double[][] getX0() {return m_X0;}
   public double[][] getX() {return m_X;}
   public double     getPotential() {return m_potential;}
   public double[][] getGradient() {return m_gradient;}

   ////////////////////////////////////////////////////////////////////////////////
   // Constructor
   public Minimizer(APotentialField pfn, double[][] X0) {
      m_field = pfn;
      int dims = X0.length;
      int n = X0[0].length;
      m_X0 = new double[dims][n];
      for (int k = 0; k < dims; ++k) System.arraycopy(X0[k], 0, m_X0[k], 0, n);
      m_X = new double[dims][n];
      for (int k = 0; k < dims; ++k) System.arraycopy(X0[k], 0, m_X[k], 0, n);
      m_gradient = new double[dims][n];
      m_gradient0 = new double[dims][n];
      m_gradientNorms = new double[n];
      m_potential = 0;
      m_calculator = m_field.potentialCalculator(null, m_X, m_gradient, m_gradientNorms);
      m_stepWorkers = new StepWorker[Util.workers()];
      for (int i = 0; i < m_stepWorkers.length; ++i) {
         m_stepWorkers[i] = new StepWorker(i, Util.workers());
      }
      report = null;
   }

   ////////////////////////////////////////////////////////////////////////////////
   // Private Helper Functions
   public double gradNormCalc(int[] subset) throws Exception {
      int[] ss = m_calculator.subset;
      m_calculator.subset = subset;
      double tmp = m_calculator.gradNorms();
      m_calculator.subset = ss;
      m_gradientTotalNorm = m_calculator.gradientLength;
      return tmp;
   }
   private void takeStep(double dt, int[] subset) throws Exception {
      ExecutorService exc = Util.pool();
      Future[] fut = new Future[m_stepWorkers.length];
      m_stepSize = dt;
      for (int i = 0; i < m_stepWorkers.length; ++i) {
         m_stepWorkers[i].subset = subset;
         m_stepWorkers[i].direction = true;
         fut[i] = exc.submit(m_stepWorkers[i]);
      }
      for (int i = 0; i < fut.length; ++i) {
         fut[i].get();
      }
   }
   private void revertStep(int[] subset) throws Exception {
      ExecutorService exc = Util.pool();
      Future[] fut = new Future[m_stepWorkers.length];
      for (int i = 0; i < m_stepWorkers.length; ++i) {
         m_stepWorkers[i].subset = subset;
         m_stepWorkers[i].direction = false;
         fut[i] = exc.submit(m_stepWorkers[i]);
      }
      for (int i = 0; i < fut.length; ++i) {
         fut[i].get();
      }
   }

   ////////////////////////////////////////////////////////////////////////////////
   // The Step Function
   /** min.step(dt, z) follows the gradient of its potential-field and configuration until it has
    *  traveled for dt units of time and such that no vertex ever moves more than distance z in a
    *  single step. 
    */
   synchronized public Report step(double deltaT, int maxSteps, double z) throws Exception {
      double maxNorm, t, t0, dt, dx, pe0;
      int k = 0;
      if (deltaT <= 0) return null;
      if (z <= 0) throw new IllegalArgumentException("parameter z to step must be > 0");
      // first thing: calculate the total gradient!
      m_calculator.clear();
      m_calculator.calculate();
      m_potential = m_calculator.potential;
      if (Double.isNaN(m_potential))
         throw new IllegalArgumentException("Initial state has a NaN potential");
      else if (Double.isInfinite(m_potential))
         throw new IllegalArgumentException("Initial state has a non-finite potential");
      // also the gradient norms...
      maxNorm = gradNormCalc(null);
      // okay, iteratively take appropriately-sized steps...
      dx = 0;
      t = 0;
      pe0 = m_potential;
      Report re = new Report(pe0);
      try {
         while (t < deltaT && k < maxSteps) {
            if (Util.zeroish(maxNorm))
               throw new Exception("gradient is effectively 0");
            // pick our start step size; first would be z/maxNorm or timeLeft, whichever is smaller
            dt = z / maxNorm;
            if (dt + t > deltaT) dt = deltaT - t;
            t0 = t;
            // see if the current step-size works; if not we'll halve it and try again...
            while (t0 == t) {
               // make sure we aren't below a threshold...
               if (Util.zeroish(dt)) throw new Exception("Step-size decreased to effectively 0");
               // take a step; this copies the current coordinates (m_X) into m_X0; same for grad
               takeStep(dt, null);
               // calculate the new gradient/potential
               m_calculator.clear();
               m_calculator.calculate();
               m_potential = m_calculator.potential;
               if (Double.isNaN(m_potential)) {
                  throw new IllegalStateException("Potential function yielded NaN");
               } else if (Double.isInfinite(m_potential) || m_potential >= pe0) {
                  // we broke a triangle or we failed to reduce potential (perhaps due to a 
                  // too-large step-size); swap x0 back to x and grad0 back to grad and try with a
                  // smaller step
                  revertStep(null);
                  dt *= 0.5;
               } else {
                  ++k;
                  re.push(dt, m_gradientTotalNorm * dt, maxNorm, m_potential - pe0);
                  // the step was okay; we can cement it
                  t += dt;
                  pe0 = m_potential;
                  // we need to get the gradient norms too
                  maxNorm = gradNormCalc(null);
                  dx += m_gradientTotalNorm * dt;
               }
            }
         }
         // this was a success, so copy X over to X0
         for (int i = 0; i < m_X.length; ++i)
            System.arraycopy(m_X[i], 0, m_X0[i], 0, m_X0[i].length);
      } finally {
         // freeze the report and save it...
         re.freeze(m_potential);
         report = re;
      }
      return re;
   }


   /** The nimbleStep function selectively updates the vertices with the highest gradient more 
    *  frequently than vertices with lower gradients while performing a time-step in order to
    *  improve accuracy at a lower cost in terms of runtime.
    */
   synchronized public double nimbleStep(double deltaT, int maxSteps, double z) throws Exception {
      double maxNorm, t, t0, dt, dt0, dx, pe0, peStep;
      int partitions = 8; // how many partitions to make...
      int miniStepsPerStep = (1 << partitions);
      AInPlaceCalculator[] calcs;
      int k = 0;
      int miniStep, part;
      boolean cont;
      if (deltaT <= 0) return Double.NaN;
      if (z <= 0) throw new IllegalArgumentException("parameter z to step must be > 0");
      // first thing: calculate the total gradient!
      m_calculator.clear();
      m_calculator.calculate();
      pe0 = m_calculator.potential;
      maxNorm = gradNormCalc(null);
      if (Double.isNaN(m_potential))
         throw new IllegalArgumentException("Initial state has a NaN potential");
      else if (Double.isInfinite(m_potential))
         throw new IllegalArgumentException("Initial state has a non-finite potential");
      // okay, iteratively take appropriately-sized overall-steps...
      dx = 0;
      t = 0;
      while (t < deltaT && k++ < maxSteps) {
         if (Util.zeroish(maxNorm))
            throw new Exception("gradient is effectively 0");
         // pick our start step size; first would be z/maxNorm or timeLeft, whichever is smaller
         dt0 = z / maxNorm;
         if (dt0 + t > deltaT) dt0 = deltaT - t;
         t0 = t;
         // we want to make some substeps to run through...
         calcs = substeps(partitions);
         // okay, we make <ministepsPerStep> steps total...
         for (miniStep = 0; miniStep < miniStepsPerStep; ++miniStep) {
            // on this mini-step, we update the appropriate subsets using a max step-size of dt
            // scaled up to be appropriate for how often this subset is updated;
            for (part = 0; part < partitions; ++part) {
               // only do this part if it is divisible by the appropriate power
               if (miniStep % (1 << part) > 0) continue;
               // we need to recalculate potential etc, as it may have changed...
               calcs[part].clear();
               calcs[part].calculate();
               calcs[part].gradNorms();
               // also skip this part if the max norm is 0
               if (Util.zeroish(calcs[part].gradientLength)) continue;
               // we always start with this stepsize...
               dt = dt0;
               // see if the current step-size works; if not we'll halve it and try again...
               peStep = calcs[part].potential;
               cont = true;
               while (cont) {
                  // make sure we aren't below a threshold...
                  if (Util.zeroish(dt))
                     throw new Exception("Step-size decreased to effectively 0 -- " + 
                                         peStep + " ;; " + part + " ;; " + miniStep + " -- " + 
                                         calcs[part].maxGradientNorm + " ;; " + 
                                         calcs[part].gradientLength + " ;; " + 
                                         calcs[part].potential + " ;; " + dt0);
                  // save the start potential...
                  // take a step; this copies the current coordinates (m_X) into m_X0; 
                  // same for grad
                  takeStep(dt, calcs[part].subset);
                  // calculate the new gradient/potential
                  calcs[part].clear();
                  calcs[part].calculate();
                  if (Double.isNaN(calcs[part].potential)) {
                     throw new IllegalStateException("Potential function yielded NaN");
                  } else if (Double.isInfinite(calcs[part].potential) 
                             || peStep < calcs[part].potential) {
                     // we broke a triangle or we failed to reduce potential (perhaps due to a 
                     // too-large step-size); swap x0 back to x and grad0 back to grad and try 
                     // with a smaller step
                     revertStep(calcs[part].subset);
                     dt *= 0.5;
                  } else {
                     // the step was okay; we can cement it
                     calcs[part].gradNorms();
                     cont = false;
                  }
               }
            }
         }
      }
      // this was a success, so copy X over to X0
      for (int i = 0; i < m_X.length; ++i)
         System.arraycopy(m_X[i], 0, m_X0[i], 0, m_X0[i].length);
      return 0;
   }

   // nimbleStep uses this function to partition the gradNorm's into a plan of action:
   /** min.substeps() examines the current gradient norms (without recalculations) and partitions 
    *  them into a sequence of 8 calculators: R; for the next 128 steps, the recommended subsets to
    *  use on step k are the subsets represented in R[i] for every i &lt; 8 such that 
    *  mod(k, 2^i) == 0.
    */
   public AInPlaceCalculator[] substeps() {return substeps(8);}
   public AInPlaceCalculator[] substeps(int k) {
      if (k < 2 || k > 32) throw new IllegalArgumentException("steps must be in the range 2-32");
      int i, j, n = m_gradientNorms.length;
      int steps = (1 << k);
      AInPlaceCalculator[] calcs = new AInPlaceCalculator[k];
      // first, sort the gradient norms... this gives us the cutoffs
      double[] gnorms = new double[n];
      System.arraycopy(m_gradientNorms, 0, gnorms, 0, gnorms.length);
      Arrays.sort(gnorms);
      // now we can go ahead and get the bucket sizes...
      int[][] buckets = new int[k][];
      int[] ranks = new int[k];
      double[] cutoffs = new double[k];
      int left = n;
      for (i = 0; i < k-1; ++i) {
         buckets[i] = new int[n / (1 << (k - i))];
         left -= buckets[i].length;
         ranks[i] = left;
         cutoffs[i] = gnorms[left]; // must be >= than cutoff[i] to be in i
      }
      buckets[k-1] = new int[left];
      ranks[k-1] = 0;
      cutoffs[k-1] = gnorms[0];
      // now we put each index in its place
      int[] counts = new int[k];
      double tmp;
      for (i = 0; i < n; ++i) {
         tmp = m_gradientNorms[i];
         for (j = 0; j < k; ++j) {
            if (j == k-1 || (tmp >= cutoffs[j] && counts[j] < buckets[j].length)) {
               buckets[j][counts[j]++] = i;
               break;
            }
         }
      }
      // okay, got the buckets; just make subset calculators...
      for (i = 0; i < k; ++i)
         calcs[i] = m_field.potentialCalculator(buckets[i], m_X, m_gradient, m_gradientNorms);
      return calcs;
   }
}

