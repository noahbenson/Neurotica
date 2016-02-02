////////////////////////////////////////////////////////////////////////////////////////////////////
// GaussianFunction.java
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

/** The GaussianFunction class simply contains the code for calculating potential and derivative 
 *  values of inverted Gaussian-like functions. Inverted Gaussian functions are used when one does
 *  not want a large value (such as a distance) to result in an enormous gradient.
 *
 *  @author Noah C. Benson
 */
public final class GaussianFunction implements IDifferentiatedFunction {
   public final double scale;
   public final double sigma;
   public final double shape;

   public GaussianFunction(double scale, double sigma, double shape) {
      this.scale = scale;
      this.sigma = sigma;
      this.shape = shape;
   }
   public GaussianFunction(double scale, double sigma) {
      this.scale = scale;
      this.sigma = sigma;
      this.shape = 2.0;
   }

   /** GaussianFunction.y(x, x0, k, s, q) yields the inverted Gaussian potential y at the value x
    *  given the minimum point x0, the scale k, the standard deviation s, and the shape q. The 
    *  formula for the potential is y = s * (1 - exp(-0.5((x - x0)/s)^2))
    *
    *  @param x the value at which to calculate the potential
    *  @param x0 the reference or 0-point of the potential
    *  @param k the scale of the potential function
    *  @param s the standard deviation/width of the Gaussian
    *  @param q the shape (exponent) of the Gaussian
    */
   static public double y(double x, double x0, double k, double s, double q) {
      x = (x - x0)/s;
      if (q == 2.0)
         return k * (1.0 - Math.exp(-0.5 * x*x));
      else
         return k * (1.0 - Math.exp(-0.5 * Math.pow(x, q)));
   }

   /** GaussianFunction.dy(x, x0, k, s, q) yields the derivative of the inverted Gaussian potential
    *  y in terms of the value x given the zero-point x0, the scale k, the standard deviation/width
    *  s, and the shape/exponent q. The formula for the derivative is 
    *  dy/dx = 0.5 q ((x - x0)/s)^(q-1) / s * exp(-0.5 * ((x-x0)/s)^q)
    *
    *  @param x the value at which to calculate the potential
    *  @param x0 the reference or 0-point of the potential
    *  @param k the scale of the potential function
    *  @param s the standard deviation/width of the Gaussian
    *  @param q the shape (exponent) of the Gaussian
    */
   static public double dy(double x, double x0, double k, double s, double q) {
      x = (x - x0)/s;
      if (q == 2.0)
         return 0.5*k*q*x * Math.exp(-0.5 * x*x);
      else
         return 0.5*k*q*Math.pow(x, q-1.0) * Math.exp(-0.5 * Math.pow(x, q));
   }

   // The IDifferentiatedFunction's
   public double y(double x, double x0) {return GaussianFunction.y(x, x0, scale, sigma, shape);}
   public double dy(double x, double x0) {return GaussianFunction.dy(x, x0, scale, sigma, shape);}
}

