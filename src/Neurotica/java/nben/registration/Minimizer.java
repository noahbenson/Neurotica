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

/** Minimizer is a class that handles the registration minimization. Minimizer, by default, uses the
 *  Util class's executor and worker suggestions.
 */
public class Minimizer {

   ////////////////////////////////////////////////////////////////////////////////
   // Private Classes
   private class GradientWorker implements Runnable {
      public final int id;
      public final int workers;
      public int[] subset;
      public double max;
      public double totalNormSquared;
      public GradientWorker(int i, int nwork) {
         id = i;
         workers = nwork;
         subset = null;
      }
      public void run() {
         max = 0;
         totalNormSquared = 0;
         if (subset == null) {
            int j;
            for (int i = id; i < m_gradientNorms.length; i += workers) {
               m_gradientNorms[i] = 0;
               for (j = 0; j < m_gradient.length; ++j)
                  m_gradientNorms[i] += m_gradient[j][i]*m_gradient[j][i];
               totalNormSquared += m_gradientNorms[i];
               m_gradientNorms[i] = Math.sqrt(m_gradientNorms[i]);
               if (m_gradientNorms[i] > max) max = m_gradientNorms[i];
            }
         } else {
            int j, u;
            for (int i = id; i < subset.length; i += workers) {
               u = subset[i];
               m_gradientNorms[u] = 0;
               for (j = 0; j < m_gradient.length; ++j)
                  m_gradientNorms[u] += m_gradient[j][u]*m_gradient[j][u];
               totalNormSquared += m_gradientNorms[u];
               m_gradientNorms[u] = Math.sqrt(m_gradientNorms[u]);
               if (m_gradientNorms[u] > max) max = m_gradientNorms[u];
            }
         }
      }
   }
   private class StepWorker implements Runnable {
      public final int id;
      public final int workers;
      public int[] subset;
      public StepWorker(int i, int nwork) {
         id = i;
         workers = nwork;
         subset = null;
      }
      public void run() {
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
   private GradientWorker[]   m_gradientWorkers;
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
      m_calculator = m_field.potentialCalculator(null, m_X, m_gradient);
      m_gradientWorkers = new GradientWorker[Util.workers()];
      m_stepWorkers = new StepWorker[Util.workers()];
      for (int i = 0; i < m_gradientWorkers.length; ++i) {
         m_gradientWorkers[i] = new GradientWorker(i, Util.workers());
         m_stepWorkers[i] = new StepWorker(i, Util.workers());
      }
      report = null;
   }

   ////////////////////////////////////////////////////////////////////////////////
   // Private Helper Functions
   private double gradNormCalc(int[] subset) throws Exception {
      double max = 0;
      ExecutorService exc = Util.pool();
      Future[] fut = new Future[m_gradientWorkers.length];
      for (int i = 0; i < m_gradientWorkers.length; ++i) {
         m_gradientWorkers[i].subset = subset;
         fut[i] = exc.submit(m_gradientWorkers[i]);
      }
      m_gradientTotalNorm = 0;
      for (int i = 0; i < fut.length; ++i) {
         fut[i].get();
         m_gradientTotalNorm += m_gradientWorkers[i].totalNormSquared;
         if (m_gradientWorkers[i].max > max) max = m_gradientWorkers[i].max;
      }
      m_gradientTotalNorm = Math.sqrt(m_gradientTotalNorm);
      return max;
   }
   private void takeStep(double dt, int[] subset) throws Exception {
      ExecutorService exc = Util.pool();
      Future[] fut = new Future[m_stepWorkers.length];
      m_stepSize = dt;
      for (int i = 0; i < m_stepWorkers.length; ++i) {
         m_stepWorkers[i].subset = subset;
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
            if (maxNorm < 1e-30)
               throw new Exception("gradient is effectively 0");
            // pick our start step size; first would be z/maxNorm or timeLeft, whichever is smaller
            dt = z / maxNorm;
            if (dt + t > deltaT) dt = deltaT - t;
            t0 = t;
            // see if the current step-size works; if not we'll halve it and try again...
            while (t0 == t) {
               // make sure we aren't below a threshold...
               if (dt <= 1e-30) throw new Exception("Step-size decreased to effectively 0");
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
                  for (int i = 0; i < m_X.length; ++i) {
                     System.arraycopy(m_X0[i], 0, m_X[i], 0, m_X[i].length);
                     System.arraycopy(m_gradient0[i], 0, m_gradient[i], 0, m_gradient[i].length);
                  }
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
}

