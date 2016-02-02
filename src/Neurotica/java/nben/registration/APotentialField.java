////////////////////////////////////////////////////////////////////////////////////////////////////
// APotentialField.java
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
import nben.registration.AInPlaceCalculator;
import java.util.concurrent.ExecutorService;

/** The APotentialField abstract class defines the kinds of queries that can be made of a potential
 *  field object. Potentials define how vertices are compelled to move/minimize in registrations.
 *
 *  @author Noah C. Benson
 */
public abstract class APotentialField {
   /** field.potentialCalculator(S, X, G) yields a potential energy calculator for the system given
    *  the mesh coordinate matrix X (size dims x vertices) such that only the vertices indexed by 
    *  the array S will be considered. Whenever the calculator updates its values, it places the 
    *  resulting gradient in the matrix G (size dims x vertices). If G is null, then a matrix is 
    *  allocated automatically. If S is null, then all vertices are assumed to be included.
    *
    *  @returns the potential energy calculator of the coordinate configuration X
    *  @param subset either null (indicating the full set) or an array of vertex indices for which
    *                the potential and gradient should be calculated; note that this does not change
    *                the size of the gradient or coordinate matrices, but it does instruct the 
    *                function to ignore certain vertices
    *  @param X the (dims x vertices)-sized matrix of vertex coordinates at which to evaluate the
    *           potential
    *  @param G the (dims x vertices)-sized output matrix to which the gradient matrix should be
    *           added (note that overwriting this matrix would be an error)
    */
   public abstract AInPlaceCalculator potentialCalculator(int[] subset, double[][] X, double[][] G);

   /** potential(subset, X, G) is simply a wrapper that constucts an in-place-calculator using
    *  potentialCalculator(subset, X, G), runs it, and returns the resulting potential value.
    */
   public double potential(int[] subset, double[][] X, double[][] G) throws Exception {
      // we basically just create a calculation object and execute it...
      AInPlaceCalculator c = this.potentialCalculator(subset, X, G);
      c.calculate();
      return c.potential;
   }
}
