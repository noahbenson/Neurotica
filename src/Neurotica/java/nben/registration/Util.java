////////////////////////////////////////////////////////////////////////////////////////////////////
// Util.java
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
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.LinkedList;
import java.util.Arrays;

/** The nben.registration.Util class is a static class with static functions that are useful for 
 *  various tasks involved in mesh registration.
 *
 *  @author Noah C. Benson
 */
public class Util {

   static final double ZERO_TOL = 1e-30;
   static final boolean zeroish(double d) {return Math.abs(d) < ZERO_TOL;}

   // Data about how to optimally multi-thread
   private static ExecutorService m_pool;
   private static int m_nthreads;
   /** Util.pool() yields an ExecutorService with Util.worker() threads */
   public synchronized static final ExecutorService pool() {return m_pool;}
   /** Util.workers() yields the optimal number of threads to use with the ExecutorService obtained
    *  by calling Util.worker().
    */
   public synchronized static final int workers() {return m_nthreads;}
   /** Sets the number of workers to n and modifies the executor service yielded by Util.pool() to
    *  be optimal for the number of workers.
    */
   public synchronized static final void setWorkers(int n) {
      m_nthreads = n;
      m_pool.shutdown();
      m_pool = Executors.newFixedThreadPool(n);
   }
   static {
      m_nthreads = Runtime.getRuntime().availableProcessors();
      m_pool = Executors.newFixedThreadPool(m_nthreads);
   }

   /** Util.buildSimplexIndex(n, data) yields an index idx for the given data such that idx[i] is an
    *  array of indices in the second dimension of data that you can find i. The assumption here is
    *  that data is an s x n array where s is the dimensionality of the simplices and n is the
    *  number of simplices. The parameter n must be the highest number in data.
    */
   public static int[][] buildSimplexIndex(int n, int[][] data) {
      int[][] index = new int[n][];
      int m = data[0].length;
      int[] u, v;
      for (int j = 0; j < data.length; ++j) {
         for (int i = 0; i < m; ++i) {
            u = index[data[j][i]];
            if (u == null) {
               u = new int[1];
               u[0] = i;
            } else {
               v = new int[u.length + 1];
               System.arraycopy(u, 0, v, 0, u.length);
               v[u.length] = i;
               u = v;
            }
            index[data[j][i]] = u;
         }
      }
      return index;
   }

   /** Util.subsampleIndex(subset, index) yields a set of all integers in the given
    *  index that are pointed to by any element of the given subset.
    */
   public static int[] subsampleIndex(int[] subset, int[][] index) {
      int i,j;
      int[] ss;
      HashSet<Integer> q = new HashSet<Integer>();
      for (i = 0; i < subset.length; ++i) {
         ss = index[subset[i]];
         if (ss != null) {
            for (j = 0; j < ss.length; ++j) 
               q.add(new Integer(ss[j]));
         }
      }
      int[] samp = new int[q.size()];
      i = 0;
      for (Iterator<Integer> it = q.iterator(); it.hasNext(); ++i)
         samp[i] = it.next().intValue();
      return samp;
   }

   /** Util.facesToEdges(faces) yields a 2 x m matrix of edges in the given 3 x p matrix of
    *  faces; m is the number of edges in the face list and p is the number of faces.
    */
   public static int[][] facesToEdges(int[][] faces) {
      if (faces == null || faces.length != 3 
          || faces[0].length != faces[2].length
          || faces[0].length != faces[1].length)
         return null;
      HashSet<List<Integer>> set = new HashSet<List<Integer>>(3 * faces[0].length);
      List<Integer> li;
      int u, v, i, j;
      for (j = 0; j < 3; ++j) {
         for (i = 0; i < faces[0].length; ++i) {
            u = faces[j][i];
            v = faces[(j + 1) % 3][i];
            if (u == v) throw new IllegalArgumentException("improper face found with form (u,u,v)");
            li = (u > v? Arrays.asList(new Integer(v), new Integer(u)) 
                       : Arrays.asList(new Integer(u), new Integer(v)));
            set.add(li);
         }
      }
      int[][] edges = new int[2][set.size()];
      i = 0;
      for (Iterator<List<Integer>> it = set.iterator(); it.hasNext(); ++i) {
         li = it.next();
         edges[0][i] = li.get(0).intValue();
         edges[1][i] = li.get(1).intValue();
      }
      return edges;
   }

   /** Util.facesToAngles(faces) yields a 3 x m matrix of angles in the given 3 x p matrix of faces
    */
   public static int[][] facesToAngles(int[][] faces) {
      if (faces == null || faces.length != 3 
          || faces[0].length != faces[2].length
          || faces[0].length != faces[1].length)
         return null;
      int n = faces[0].length;
      int i, j, k;
      int[][] angles = new int[3][3*n];
      for (j = 0; j < 3; ++j) { // angle number...
         for (k = 0; k < 3; ++k) { // angle part (A, B, C)
            for (i = 0; i < n; ++i) { // face number...
               angles[k][i + n*j] = faces[(j + k) % 3][i];
            }
         }
      }
      return angles;
   }

   /** Given a list of faces (3 x m), yields a list of the vertices in the perimeter in either a
    *  clockwise or counter-clockwise ordering. If the mesh is actually multiple meshes, or if the
    *  faces do not form a proper mesh, an exception is thrown.
    */
   public static int[] perimeter(int[][] faces) {
      if (faces == null || faces.length != 3 
          || faces[0].length != faces[2].length
          || faces[0].length != faces[1].length)
         return null;
      HashSet<List<Integer>> set = new HashSet<List<Integer>>(faces[0].length * 2);
      int u, v, i, j;
      List<Integer> li;
      for (j = 0; j < 3; ++j) {
         for (i = 0; i < faces[0].length; ++i) {
            u = faces[j][i];
            v = faces[(j + 1) % 3][i];
            if (u == v) throw new IllegalArgumentException("improper face found with form (u,u,v)");
            li = (u > v? Arrays.asList(new Integer(v), new Integer(u)) 
                       : Arrays.asList(new Integer(u), new Integer(v)));
            // we want these to stick around only if they occur an odd number of times
            if (set.contains(li)) set.remove(li);
            else                  set.add(li);
         }
      }
      // anything that is still in the set is part of the perimeter of the mesh, assuming that the
      // faces were correctly specified.
      // we can traverse them by making a hash then following them around...
      int n = set.size();
      Integer g, h, tmpint, end, front;
      LinkedList<Integer> lg, lh, lt = null;
      Iterator<Integer> iit;
      HashMap<Integer,LinkedList<Integer>> map = new HashMap<Integer,LinkedList<Integer>>(2*n);
      for (Iterator<List<Integer>> it = set.iterator(); it.hasNext(); ) {
         li = it.next();
         g = li.get(0);
         h = li.get(1);
         lg = map.get(g);
         lh = map.get(h);
         if (lg != null) map.remove(g);
         if (lh != null) map.remove(h);
         if (lg == null && lh == null) {
            // make a new pair of lists
            lt = new LinkedList<Integer>();
            lt.addFirst(g);
            lt.addLast(h);
         } else if (lh == null || lg == null) {
            tmpint = g;
            g = (lh == null? g : h);
            h = (lh == null? h : tmpint);
            lt = (lh == null? lg : lh);
            // now, regardless of ordering, g is the found one, h is the null one, and lt is
            // the list mapped to g; so we are adding h to lt
            tmpint = lt.getFirst();
            if (tmpint.equals(g))
               lt.addFirst(h);
            else
               lt.addLast(h);
         } else if (!lg.equals(lh)) {
            // we link two pieces together; first, we want to add the smaller to the larger
            if (lg.size() > lh.size()) {
               lt = lg;
               lg = lh;
               lh = lt;
               tmpint = g;
               g = h;
               h = tmpint;
            }
            // now, lg is smaller than lh and h/g match lh/lg; see if we are adding from the front
            // or the back of lg:
            iit = (lg.getFirst().equals(g)? lg.iterator() : lg.descendingIterator());
            // and to the front or back of lh...
            if (lh.getFirst().equals(h))
               while (iit.hasNext())
                  lh.addFirst(iit.next());
            else
               while (iit.hasNext())
                  lh.addLast(iit.next());
            lt = lh;
         // any clause below here assumes that lh == lg and g<->h completes the circle...
         } else if (it.hasNext()) {
            // there are more edges; so there's a bug in the face matrix
            throw new IllegalArgumentException("face matrix contained multiple meshes");
         } else {
            // everything seems hunky-dory! we can break out of this loop
            lt = lg;
            break;
         }
         map.put(lt.getFirst(), lt);
         map.put(lt.getLast(), lt);
      }
      if (lt == null) throw new IllegalStateException("lt == null?!");
      // now, we just need to make sure that the edges get written out in the right order; lt stores
      // the vertices...
      int[] result = new int[lt.size()];
      // we want to set iit to either a forward or reverse iterator, so see how the first<->last 
      // edge is ordered in the face matrix
      u = lt.getFirst().intValue();
      v = lt.getLast().intValue();
      for (i = 0; i < faces[0].length; ++i) {
         j = ((faces[0][i] == u || faces[0][i] == v? 1 : 0)
              + (faces[1][i] == u || faces[1][i] == v? 1 : 0)
              + (faces[2][i] == u || faces[2][i] == v? 1 : 0));
         if (j == 2) break;
      }
      if (i == faces[0].length)
         throw new IllegalStateException("initial edge not found in faces (" + u + ", " + v + ")");
      iit = (faces[0][i] == u? (faces[1][i] == v? lt.descendingIterator() : lt.iterator())
             : (faces[1][i] == u? (faces[2][i] == v? lt.descendingIterator() : lt.iterator())
                : (faces[0][i] == v? lt.descendingIterator() : lt.iterator())));
      for (i = 0; iit.hasNext(); ++i)
         result[i] = iit.next().intValue();
      return result;
   }
}