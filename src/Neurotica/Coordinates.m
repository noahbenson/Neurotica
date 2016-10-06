(* Coordinates.m
 *
 * The Neurotica`Coordinates package includes code for switching between different coordinate
 * systems used for representing the cortical surface.
 *
 * Copyright (C) 2014-2015 by Noah C. Benson.
 * This file is part of the Neurotica library.
 *
 * This program is free software: you can redistribute it and/or modify it under the terms of the
 * GNU General Public License as published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See
 * the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with this program.
 * If not, see <http://www.gnu.org/licenses/>.
 *)

(**************************************************************************************************)
BeginPackage["Neurotica`Coordinates`", {"Neurotica`Global`", "Neurotica`Util`"}];
Unprotect["Neurotica`Coordinates`*", "Neurotica`Coordinates`Private`*"];
ClearAll["Neurotica`Coordinates`*", "Neurotica`Coordinates`Private`*"];

Longitude::usage = "Longitude is a keyword that represents the longitude (identical to azimuth angle) in spherical coordinates.";
Latitude::usage = "Latitude is a keyword that represents the latitude in spherical coordinates.";

CartesianToSpherical::usage = "CartesianToSpherical[data] yields an equivalent dataset to data except that the parts of data corresponding to (x,y,z) coordinates are replaced with (\[Theta],\[Phi]) coordinates. The data parameter must be a list of either lists or rules OR a single list or a single rule. In data is a list of lists, then the first three elements of each sublist is taken to be the (x,y,z) part; if data is a list of rules or a rule, then the first element of the rule is taken to be the (x,y,z) part. CartesianToSpherical yields data in the same format as it is given, but with the appropriate replacements."
CartesianToSpherical::badfmt = "Bad data format given to CartesianToSpherical";
CartesianToSphercial::badopt = "Bad option `1` given to CartesianToSpherical: `2`";

SphericalToCartesian::usage = "SphericalToCartesian is identical to CartesianToSpherical, but performs the inverse transform."
SphericalToCartesian::badfmt = "Bad data format given to SphericalToCartesian";
SphercialToCartesian::badopt = "Bad option `1` given to SphericalToCartesian: `2`";

SphericalCoordinateStyle::usage = "SphericalCoordinateStyle is an option for functions such as Surface that accept spherical coordinates. SphericalCoordinateStyle describes coordinates of the form {a, b} such that a SphericalCoordinateStyle value of {x, y} indicates that a represents measurement x and b represents measurement y. The possible measurements are Azimuth/Longitude (which are identical), Latitude, and PolarAngle. At least one of x and y must be Azimuth or Longitude and one of them must be either Latitude or PolarAngle. The ConvertCoordinates function can be used to convert between representations.";
SphericalPolarAngle::usage = "SphericalPolarAngle is a keyword that either represents polar angle data of the visual field or the angle from the pole (z-axis) in spherical coordinate data.";
SphericalAzimuth::usage = "SphericalAzimuth is a keyword that represents the azimuth angle (identical to longitude) in spherical coordinates.";

ConvertCoordinates::usage = "ConvertCoordinates[data, fromStyle -> toStyle] yields a transformation of data such that its coordinates are converted from fromStyle to toStyle; fromStyle and toStyle may be Cartesian or any argument appropriate for SphericalCoordinateStyle. The data parameter may be a list of points or a list of rules whose heads are points.";
ConvertCoordinates::badfmt = "Bad input to ConvertCoordinates: `1`";

Cartesian::usage = "Cartesian is a keyword used by the Neurotica`Coordinates` package to denote cartesian coordinates.";

Begin["`Private`"];

Protect[Latitude, Longitude, SphericalCoordinateStyle, SphericalPolarAngle, SphericalAzimuth, 
        Cartesian];

(* #CartesianToSpherical **************************************************************************)
Options[CartesianToSpherical] = {SphericalCoordinateStyle -> {Longitude, Latitude}};
CartesianToSpherical[dat_, OptionsPattern[]] := (Message[CartesianToSpherical::badfmt]; $Failed);
CartesianToSpherical[dat : {{_, _, _} ..}, OptionsPattern[]] := With[
  {scFn = Replace[
     OptionValue[SphericalCoordinateStyle],
     {{Latitude, Longitude|SphericalAzimuth} -> Function[
        {ArcSin[#[[3]]/Norm[#]],
         If[#[[1]] == 0 && #[[2]] == 0, 0, ArcTan[#[[1]], #[[2]]]]}],
      {Longitude|SphericalAzimuth, Latitude} -> Function[
        {If[#[[1]] == 0 && #[[2]] == 0, 0, ArcTan[#[[1]], #[[2]]]],
         ArcSin[#[[3]]/Norm[#]]}],
      {SphericalPolarAngle, Longitude|SphericalAzimuth} -> Function[
        {ArcCos[#[[3]]/Norm[#]],
         If[#[[1]] == 0 && #[[2]] == 0, 0, ArcTan[#[[1]], #[[2]]]]}],
      {Longitude|SphericalAzimuth, SphericalPolarAngle} -> Function[
        {If[#[[1]] == 0 && #[[2]] == 0, 0, ArcTan[#[[1]], #[[2]]]],
         ArcCos[#[[3]]/Norm[#]]}],
      _ :> (
        Message[CartesianToSpherical::badopt, SphericalCoordinateStyle, "Style not recognized"];
        $Failed)}]},
  If[scFn === $Failed, $Failed, Map[scFn, dat]]];
CartesianToSpherical[dat : {Rule[{_, _, _}, _] ..}, OptionsPattern[]] := With[
  {scFn = Replace[
     OptionValue[SphericalCoordinateStyle],
     {{Latitude, Longitude|SphericalAzimuth} -> Function[
        Rule[
          {ArcSin[#[[1,3]]/Norm[#[[1]]]],
           If[#[[1,1]] == 0 && #[[1,2]] == 0, 0, ArcTan[#[[1,1]], #[[1,2]]]]},
          #[[2]]]],
      {Longitude|SphericalAzimuth, Latitude} -> Function[
        Rule[
          {If[#[[1,1]] == 0 && #[[1,2]] == 0, 0, ArcTan[#[[1,1]], #[[1,2]]]],
           ArcSin[#[[1,3]]/Norm[#[[1]]]]},
          #[[2]]]],
      {SphericalPolarAngle, Longitude|SphericalAzimuth} -> Function[
        Rule[
          {ArcCos[#[[1,3]]/Norm[#[[1]]]],
           If[#[[1,1]] == 0 && #[[1,2]] == 0, 0, ArcTan[#[[1,1]], #[[1,2]]]]},
          #[[2]]]],
      {Longitude|SphericalAzimuth, SphericalPolarAngle} -> Function[
        Rule[
          {If[#[[1,1]] == 0 && #[[1,2]] == 0, 0, ArcTan[#[[1,1]], #[[1,2]]]],
           ArcCos[#[[1,3]]/Norm[#[[1]]]]},
          #[[2]]]],
      _ :> (
        Message[CartesianToSpherical::badopt, SphericalCoordinateStyle, "Style not recognized"];
        $Failed)}]},
  If[scFn === $Failed, $Failed, Map[scFn, dat]]];
CartesianToSpherical[dat : {_?AtomQ, _?AtomQ, _?AtomQ}, OptionsPattern[]] := First[
  CartesianToSpherical[
    {dat},
    SphericalCoordinateStyle -> OptionValue[SphericalCoordinateStyle]]];
CartesianToSpherical[dat : Rule[{_?AtomQ, _?AtomQ, _?AtomQ}, _], OptionsPattern[]] := First[
  CartesianToSpherical[
    {dat},
    SphericalCoordinateStyle -> OptionValue[SphericalCoordinateStyle]]];
Protect[CartesianToSpherical];

(* #SphericalToCartesian **************************************************************************)
Options[SphericalToCartesian] = {
  SphericalCoordinateStyle -> {Longitude, Latitude}};
SphericalToCartesian[dat_, OptionsPatterm[]] := (Message[SphericalToCartesian::badfmt]; $Failed);
SphericalToCartesian[dat : {{_, _} ..}, OptionsPattern[]] := With[
  {scFn = Replace[
     OptionValue[SphericalCoordinateStyle],
     {{Latitude, Longitude|SphericalAzimuth} -> Function[
        Append[Cos[#[[1]]] * {Cos[#[[2]]], Sin[#[[2]]]}, Sin[#[[1]]]]],
      {Longitude|SphericalAzimuth, Latitude} -> Function[
        Append[Cos[#[[2]]] * {Cos[#[[1]]], Sin[#[[1]]]}, Sin[#[[2]]]]],
      {SphericalPolarAngle, Longitude|SphericalAzimuth} -> Function[
        Append[Sin[#[[1]]] * {Cos[#[[2]]], Sin[#[[2]]]}, Cos[#[[1]]]]],
      {Longitude|SphericalAzimuth, SphericalPolarAngle} -> Function[
        Append[Sin[#[[2]]] * {Cos[#[[1]]], Sin[#[[1]]]}, Cos[#[[2]]]]],
      _ :> (
        Message[SphericalToCartesian::badopt, SphericalCoordinateStyle, "Style not recognized"];
        $Failed)}]},
  If[scFn === $Failed, $Failed, Map[scFn, dat]]];
SphericalToCartesian[dat : {Rule[{_, _}, _] ..}, OptionsPattern[]] := With[
  {scFn = Replace[
     OptionValue[SphericalCoordinateStyle],
     {{Latitude, Longitude|SphericalAzimuth} -> Function[
        Rule[
          Append[Cos[#[[1,1]]] * {Cos[#[[1,2]]], Sin[#[[1,2]]]}, Sin[#[[1,1]]]],
          #[[2]]]],
      {Longitude|SphericalAzimuth, Latitude} -> Function[
        Rule[
          Append[Cos[#[[1,2]]] * {Cos[#[[1,1]]], Sin[#[[1,1]]]}, Sin[#[[1,2]]]],
          #[[2]]]],
      {SphericalPolarAngle, Longitude|SphericalAzimuth} -> Function[
        Rule[
          Append[Sin[#[[1,1]]] * {Cos[#[[1,2]]], Sin[#[[1,2]]]}, Cos[#[[1,1]]]],
          #[[2]]]],
      {Longitude|SphericalAzimuth, SphericalPolarAngle} -> Function[
        Rule[
          Append[Sin[#[[1,2]]] * {Cos[#[[1,1]]], Sin[#[[1,1]]]}, Cos[#[[1,2]]]],
          #[[2]]]],
      _ :> (
        Message[SphericalToCartesian::badopt, SphericalCoordinateStyle, "Style not recognized"];
        $Failed)}]},
  If[scFn === $Failed, $Failed, Map[scFn, dat]]];
SphericalToCartesian[dat : {_?AtomQ, _?AtomQ}, OptionsPattern[]] := SphericalToCartesian[
  {dat}, 
  SphericalCoordinateStyle -> OptionValue[SphericalCoordinateStyle]];
SphericalToCartesian[dat : Rule[{_, _}, _], OptionsPattern[]] := SphericalToCartesian[
  {dat}, 
  SphericalCoordinateStyle -> OptionValue[SphericalCoordinateStyle]];
Protect[SphericalToCartesian];

(* #ConvertCoordinates ****************************************************************************)
ConvertCoordinates[data:{Rule[{_,_,_},_]..}, Rule[Cartesian, toStyle_]] := If[
  toStyle === Cartesian, data,
  CartesianToSpherical[data, SphericalCoordinateStyle -> toStyle]];
ConvertCoordinates[data:{{_,_,_}..}, Rule[Cartesian, toStyle_]] := If[
  toStyle === Cartesian, data,
  CartesianToSpherical[data, SphericalCoordinateStyle -> toStyle]];

ConvertCoordinates[data:{Rule[{_,_},_]..}, Rule[fromStyle_, Cartesian]] := SphericalToCartesian[
  data,
  SphericalCoordinateStyle -> fromStyle];
ConvertCoordinates[data:{{_,_}..}, Rule[fromStyle_, Cartesian]] := SphericalToCartesian[
  data,
  SphericalCoordinateStyle -> fromStyle];
ConvertCoordinates[data:{Rule[{_,_},_]..}, Rule[Cartesian, toStyle_]] := CartesianToSpherical[
  data,
  SphericalCoordinateStyle -> toStyle];
ConvertCoordinates[data:{{_,_}..}, Rule[Cartesian, toStyle_]] := CartesianToSpherical[
  data,
  SphericalCoordinateStyle -> toStyle];

ConvertCoordinates[data:{Rule[{_,_},_]..},
                   Rule[fromStyle:Except[Cartesian], toStyle:Except[Cartesian]]
                   ] := CartesianToSpherical[
  SphericalToCartesian[data, SphericalCoordinateStyle -> fromStyle],
  SphericalCoordinateStyle -> toStyle];
ConvertCoordinates[data:{{_,_}..},
                   Rule[fromStyle:Except[Cartesian], toStyle:Except[Cartesian]]
                   ] := CartesianToSpherical[
  SphericalToCartesian[data, SphericalCoordinateStyle -> fromStyle],
  SphericalCoordinateStyle -> toStyle];

ConvertCoordinates[data_List, _] := (
  Message[ConvertCoordinates::badfmt, "Cannot deduce conversion of data"];
  $Failed);
Protect[ConvertCoordinates];

End[];
EndPackage[];
