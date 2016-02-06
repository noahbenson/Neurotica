////////////////////////////////////////////////////////////////////////////////////////////////////
// LennardJonesFunction.java
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

/** LennardJonesFunction class simply contains the code for calculating potential and derivative 
 *  values of Lennard-Jones-like functions. Lennard-Jones functions are used in modeling van der
 *  Waals forces and have the feature that they are minimal at a small positive x and become
 *  infinite as x approaches 0, but asymptotically approach a constant as x approaches infinity.
 *
 *  @author Noah C. Benson
 */
public final class LennardJonesFunction implements IDifferentiatedFunction {
   public final double scale;
   public final double order;

   public LennardJonesFunction(double scale, double order) {
      this.scale = scale;
      this.order = order;
   }
   public LennardJonesFunction(double scale) {
      this.scale = scale;
      this.order = 2.0;
   }
   public LennardJonesFunction() {
      this.scale = 1.0;
      this.order = 2.0;
   }

   /** LennardJonesFunction.y(x, x0, s, q) yields the Lennard-Jones potential y at the value x
    *  given the minimum point x0, the scale s, and the shape q. The formula for the potential is
    *  y = s* (1 + (r0/r)^(q) - 2*(r0/r)^(q/2)).
    *
    *  @param x the value at which to calculate the potential
    *  @param x0 the reference or 0-point of the potential
    *  @param s the scale of the Lennard-Jones potential function
    *  @param q the shape (exponent) of the Lennard-Jones potential function
    */
   static public double y(double x, double x0, double s, double q) {
      x = x0 / x;
      if (x == 0) return Double.POSITIVE_INFINITY;
      return s * (1 + Math.pow(x, q) - 2.0*Math.pow(x, 0.5*q));
   }

   /** LennardJonesFunction.dy(x, x0, s, q) yields the derivative of the Lennard-Jones potential y
    *  in terms of the value x given the zero-point x0, the scale s, and the shape q.  The formula
    *  for the derivative is dy/dx = -s*q/r * ((r0/r)^q - (r0/r)^(q/2)).
    *
    *  @param x the value at which to calculate the potential
    *  @param x0 the reference or 0-point of the potential
    *  @param s the scale of the Lennard-Jones potential function
    *  @param q the shape (exponent) of the Lennard-Jones potential function
    */
   static public double dy(double x, double x0, double s, double q) {
      x0 /= x;
      if (x0 == 0) return Double.POSITIVE_INFINITY;
      return -s*q/x * (Math.pow(x0, q) - Math.pow(x0, 0.5*q));
   }

   // The IDifferentiatedFunction's
   public double y(double x, double x0) {return LennardJonesFunction.y(x, x0, scale, order);}
   public double dy(double x, double x0) {return LennardJonesFunction.dy(x, x0, scale, order);}
}

