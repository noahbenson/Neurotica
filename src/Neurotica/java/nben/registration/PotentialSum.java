////////////////////////////////////////////////////////////////////////////////////////////////////
// PotentialSum.java
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

import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;

/** The PotentialSum class aggregates a set of potential field objects and acts as a single field
 *  equal to the sum of these fields.
 *
 *  @author Noah C. Benson
 */
class PotentialSum extends APotentialField {

   private APotentialField[] m_terms;

   /** Constructs a PotentialSum object from the give list of fields */
   public PotentialSum(APotentialField... fields) {
      if (fields == null) {
         // copy the array just in case...
         m_terms = new APotentialField[fields.length];
         for (int i = 0; i < fields.length; ++i)
            m_terms[i] = fields[i];
      }
   }

   /** Adds a potential field to this object */
   public void addField(APotentialField f) {
      APotentialField[] fs;
      if (m_terms == null) {
         m_terms = new APotentialField[1];
         m_terms[0] = f;
      } else {
         fs = new APotentialField[m_terms.length + 1];
         System.arraycopy(m_terms, 0, fs, 0, m_terms.length);
         fs[m_terms.length] = f;
         m_terms = fs;
      }
   }

   /** Remove a potential field from this object */
   public void removeField(APotentialField f) {
      APotentialField[] fs;
      if (m_terms != null) {
         int i;
         for (i = 0; i < m_terms.length && m_terms[i] != f; ++i)
            ;
         if (i < m_terms.length) {            
            fs = new APotentialField[m_terms.length - 1];
            if (i > 0)
               System.arraycopy(m_terms, 0, fs, 0, i);
            if (i + 1 < m_terms.length)
               System.arraycopy(m_terms, i+1, fs, i, m_terms.length - i - 1);
            m_terms = fs;
         }
      }
   }

   /** Yields the list of potential fields tracked by this summation */
   public APotentialField[] terms() {return m_terms;}

   // Here we have the code/subclasses that handle the workers
   private final class SumCalculation extends AInPlaceCalculator {
      // all the data we need to store...
      private AInPlaceCalculator[] m_calcs;

      public final class SumConstructor implements Runnable {
         public final int id;
         public final int workers;
         public final int[] subset;
         public final double[][] X0;
         public final double[][] G;
         public SumConstructor(int i, int nworkers, int[] ss, double[][] X0, double[][] G) {
            id = i;
            workers = nworkers;
            subset = ss;
            this.X0 = X0;
            this.G = G;
         }
         public void run() {
            for (int i = id; i < m_calcs.length; i += workers)
               m_calcs[i] = m_terms[i].potentialCalculator(subset, X0, G);
         }
      }

      // constructor
      public SumCalculation(int[] ss, double[][] X0, double[][] G) {
         super(ss, X0, G);
         // if ss is null, we basically know that construction will be super-fast;
         // if not, we multi-thread it to make it a little faster
         if (m_terms == null) {
            m_calcs = null;
            return;
         }
         m_calcs = new AInPlaceCalculator[m_terms.length];
         if (ss == null) {
            for (int i = 0; i < m_terms.length; ++i)
               m_calcs[i] = m_terms[i].potentialCalculator(ss, X0, G);
         } else {
            ExecutorService exc = PotentialFields.pool();
            int nworkers = PotentialFields.workers();
            SumConstructor[] scs = new SumConstructor[nworkers];
            for (int i = 0; i < nworkers; ++i)
               scs[i] = new SumConstructor(i, nworkers, ss, X0, G);
            try {
               runThreads(scs);
            } catch (Exception e) {
               m_calcs = null;
            }
         }
      }

      // Actually run the calculations...
      public void calculate(int nworkers, ExecutorService exc) throws Exception {
         if (m_calcs == null) throw new Exception("SumCalculation failed to initialize!");
         potential = 0.0;
         // these multi-thread by themselves, so we can run them iteratively
         for (int i = 0; i < m_calcs.length; ++i) {
            m_calcs[i].calculate(nworkers, exc);
            // if we get an infinite or nan value, we can break; these are errors essentially
            if (Double.isInfinite(m_calcs[i].potential) || Double.isNaN(m_calcs[i].potential)) {
               potential = m_calcs[i].potential;
               break;
            } else
               potential += m_calcs[i].potential;
         }
      }
   }

   public final SumCalculation potentialCalculator(int[] subset, double[][] X, double[][] G) {
      return new SumCalculation(subset, X, G);
   }

}