////////////////////////////////////////////////////////////////////////////////////////////////////
// InfiniteWellFunction.java
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
import nben.registration.IDifferentiatedFunction;

/** The InfiniteWellFunction class simply contains the code for calculating potential and derivative
 *  values of inescepable well-like functions. Infinite well functions are used to create 
 *  singularities at certain boundaries in a range; this is most commonly used for preventing
 *  triangles from inverting via boundaries on 0 and 180 degree angles.
 *
 *  @author Noah C. Benson
 */
public final class InfiniteWellFunction implements IDifferentiatedFunction {
   public final double scale;
   public final double min;
   public final double max;
   public final double order;

   public InfiniteWellFunction(double scale, double min, double max, double order) {
      this.scale = scale;
      this.min = min;
      this.max = max;
      this.order = order;
   }
   public InfiniteWellFunction(double scale, double min, double max) {
      this.scale = scale;
      this.min = min;
      this.max = max;
      this.order = 0.5;
   }
   public InfiniteWellFunction(double scale) {
      this.scale = scale;
      this.min = 0.0;
      this.max = Math.PI;
      this.order = 0.5;
   }
   public InfiniteWellFunction() {
      this.scale = 1.0;
      this.min = 0.0;
      this.max = Math.PI;
      this.order = 0.5;
   }

   /** InfiniteWellFunction.y(x, x0, s, mn, mx, q) yields the infinite-well potential y at the value
    *  x given the minimum point x0, the scale s, the min/max mn an mx, and the order q. The formula
    *  for the potential is 
    *  y = s*((((t0 - min)/(t - min))^q - 1)^2 + (((max - t0)/(max - t))^q - 1)^2)
    *
    *  @param x the value at which to calculate the potential
    *  @param x0 the reference or 0-point of the potential
    *  @param s the scale of the potential well
    *  @param mn the minimum of the potential well
    *  @param mx the maximum of the potential well
    *  @param q the shape (exponent) of the potential well
    */
   static public double y(double x, double x0, double s, double min, double max, double q) {
      double tl = (x0 - min)/(x - min), tr = (max - x0)/(max - x);
      if (x <= min || x >= max) {
         return Double.POSITIVE_INFINITY;
      } else if (q == 0.5) {
         tr = Math.sqrt(tr) - 1.0;
         tl = Math.sqrt(tl) - 1.0;
      } else {
         tr = Math.pow(tr, q) - 1.0;
         tl = Math.pow(tl, q) - 1.0;
      }
      return s*(tr*tr + tl*tl);
   }

   /** InfiniteWellFunction.dy(x, x0, s, mn, mx, q) yields the derivative of the infinite-well 
    *  potential y in terms of the value x given the zero-point x0, the scale s, the min and max 
    *  mn and mx, and the shape q.
    *
    *  @param x the value at which to calculate the potential
    *  @param x0 the reference or 0-point of the potential
    *  @param s the scale of the potential well
    *  @param mn the minimum of the potential well
    *  @param mx the maximum of the potential well
    *  @param q the shape (exponent) of the potential well
    */
   static public double dy(double x, double x0, double s, double min, double max, double q) {
      double tl = (x0 - min)/(x - min), tr = (max - x0)/(max - x);
      if (x <= min || x >= max) {
         return Double.POSITIVE_INFINITY;
      } if (q == 0.5) {
         tl = Math.sqrt(tl);
         tr = Math.sqrt(tr);
      } else {
         tl = Math.pow(tl, q);
         tr = Math.pow(tr, q);
      }
      return 2.0*s*q * (tl*(tl - 1.0)/(min - x) + tr*(tr - 1.0)/(max - x));
   }

   // The IDifferentiatedFunction's
   public double y(double x, double x0) {
      return InfiniteWellFunction.y(x, x0, scale, min, max, order);
   }
   public double dy(double x, double x0) {
      return InfiniteWellFunction.dy(x, x0, scale, min, max, order);
   }
}

