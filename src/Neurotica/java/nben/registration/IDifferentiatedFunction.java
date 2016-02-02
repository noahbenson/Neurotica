////////////////////////////////////////////////////////////////////////////////////////////////////
// IDifferentiatedFunction.java
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

/** An IDifferentiatedFunction is an object that accepts a value x and a reference x0, returns a
 *  single value y(x, x0), and can additionally return the derivative dy(x, x0).
 *
 *  @author Noah C. Benson
 */
public interface IDifferentiatedFunction {
   public double y(double x, double x0);
   public double dy(double x, double x0);
}



