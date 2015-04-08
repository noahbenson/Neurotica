(* MRImage.m
 *
 * Utility functions for dealing with (mostly FreeSurfer) cortical volume image data in Mathematica.
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
BeginPackage["Neurotica`MRImage`", {"Neurotica`Global`", "Neurotica`Util`"}];
Unprotect["Neurotica`MRImage`*", "Neurotica`MRImage`Private`*"];
ClearAll["Neurotica`MRImage`*", "Neurotica`MRImage`Private`*"];

MeshVertexToVoxelIndex::usage = "MeshVertexToVoxelIndex[p, voldims] gives a translation of the surface point p ({x, y, z}) to an index ({i, j, k}) such that rounding the index values will give an approximate position of the surface point in the volume with dims given in voldims. If not provided, then voldims is taken to ba {256, 256, 256}. The first argument may also be a list of points, in which case the equivalent list of indices is returned.
Note that this method is intended to work with typically processed FreeSurfer volumes and surfaces and is not designed for other programs or for volumes that are used as input to FreeSurfer.";

VoxelToCoordinateMatrix::usage = "VoxelToCoordinateMatrix[vol, coordinate] yields the 3 x 4 matrix that can be used to translate voxel indices to coordinates.";

VoxelToCoordinate::usage = "VoxelToCoordinate[vol, coordinate] yields a translation of the given voxel index to an (x,y,z) coordinate, according to the volume vol's orientation matrix.";
VoxelToCoordinate::badarg = "Bad argument given to VoxelToCoordinate: `1`";

CoordinateToVoxelMatrix::usage = "CoordinateToVoxelMatrix[vol, coordinate] yields the 3 x 4 matrix that can be used to translate coordinates to voxel indices.";

CoordinateToVoxel::usage = "CoordinateToVoxel[vol, coordinate] yields a translation of the given (x,y,z) coordinate to a position in voxel space, according to the volume vol's orientation matrix.";
CoordinateToVoxel::badarg = "Bad argument given to CoordinateToVoxel: `1`";

OrientationMatrix::usage = "OrientationMatrix is an option to CortivalImage that specifies the orientation of the given volume. The matrix should be a 3 x 4 matrix or Automatic (in which case the default matrix is {{-1,0,0,0},{0,0,1,0},{0,-1,0,0}}) in which the final column is the displacement.
OrientationMatrix[vol] yields the orientation matrix for the given cortical volume object vol.";

InverseOrientationMatrix::usage = "InverseOrientationMatrix[vol] yields the inverse of the orientation matrix for the given cortical volume object vol.";

MRImage3D::usage = "MRImage3D[data] yields an Image3D-like form that can be used with Image3D functions but which stores additional relevant data regarding MR images.";
MRImage3D::badarg = "Bad argument given to MRImage3D: `1`";

MRImage::usage = "MRImage[data] yields an Image-like form that can be used with Image functions but which stores additional relevant data regarding MR images.
Note that MRImage is intended for MRI slices rather than for 3D images; for a 3D image representation, use MRImage3D.";
MRImage::badarg = "Bad argument given to MRImage: `1`";

MRImageMax::usage = "MRImageMax[mrimg] yields the maximum value of all voxels in the given MRImage object mrimg.";
MRImageMin::usage = "MRImageMin[mrimg] yields the minimum value of all voxels in the given MRImage object mrimg.";

MRImageQ::usage = "MRImageQ[img] yields True if img is a valid MRImage3D object and yields False otherwise.";
MRImageSliceQ::usage = "MRImageSliceQ[img] yields True if img is a valid MRImage object and yields False otherwise.";
MRImageObjectQ::usage = "MRImageObjectQ[img] is equivalent to Or[MRImageQ[img], MRImageSliceQ[img]].";

(**************************************************************************************************)
Begin["`Private`"];

(* #MRImage3D *************************************************************************************)

(* Because the MRImage3D and MRImage objects are so similar, we have a lot of overlapping functions,
 * which we store in this large Hold form for use in DefineImmutable[...].
 *)
$MRImageSharedMethods = Hold[
  (* Access the scale factors of min and max for the data *)
  MRImageMax[img] -> Max[ImageData[img]],
  MRImageMin[img] -> Min[ImageData[img]],
  
  (* Simple accessors for Image3D computed data *)
  ImageAspectRatio[img, opts___] := ImageAspectRatio[Image3D[img], opts],
  ImageChannels[img, opts___] := ImageChannels[Image3D[img], opts],
  ImageColorSpace[img, opts___] := ImageColorSpace[Image3D[img], opts],
  ImageDimensions[img, opts___] := ImageDimensions[Image3D[img], opts],
  ImageLevels[img, opts___] := ImageLevels[Image3D[img], opts],
  ImageMeasurements[img, opts___] := ImageMeasurements[Image3D[img], opts],

  (* Directly access values *)
  ImageValue[img, opts___] := ImageValue[Image3D[img], opts],
  PixelValue[img, opts___] := PixelValue[Image3D[img], opts],
  ImageValuePosition[img, opts___] := ImageValuePosition[Image3D[img], opts],
  PixelValuePosition[img, opts___] := PixelValuePosition[Image3D[img], opts],

  (* More complex query operations over the image *)
  ImageAccumulate[img, opts___] := ImageAccumulate[Image3D[img], opts],
  ImageHistogram[img, opts___] := ImageHistogram[Image3D[img], opts],
  DominantColors[img, opts___] := DominantColors[Image3D[img], opts],
  ImageCooccurance[img, opts___] := ImageCooccurance[Image3D[img], opts],
  BorderDimensions[img, opts___] := BorderDimensions[Image3D[img], opts],
  BinaryImageQ[img, opts___] := BinaryImageQ[Image3D[img], opts],
  AlphaChannel[img, opts___] := AlphaChannel[Image3D[img], opts],
  FindThreshold[img, opts___] := Plus[
    MRImageMin[img],
    (MRImageMax[img] - MRImageMin[img]) * FindThreshold[Image3D[img], opts]],

  (* Image Transformations *)
  Blur[img, opts___] := With[
    {res = Check[Blur[Image3D[img], opts], $Failed]},
    If[res === $Failed, 
      res,
      With[
        {newdat = ImageData[res] * (MRImageMax[img] - MRImageMin[img]) + MRImageMin[img],
         theOpts = Options[res]},
        MRImage3D[newdat, theOpts]]]],
  ImageAlign[img, opts___] := ImageAlign[Image3D[img], opts],
  ImageConvolve[img, opts___] := With[
    {res = Check[ImageConvolve[Image3D[img], opts], $Failed]},
    If[res === $Failed, 
      res,
      With[
        {newdat = ImageData[res] * (MRImageMax[img] - MRImageMin[img]) + MRImageMin[img],
         theOpts = Options[res]},
        MRImage3D[newdat, theOpts]]]],
  ImageCorrelate[img, opts___] := ImageCorrelate[Image3D[img], opts],
  ImageCorrespondingPoints[img, opts___] := ImageCorrespondingPoints[Image3D[img], opts],
  ImageCrop[img, opts___] := With[
    {res = Check[ImageCrop[Image3D[img], opts], $Failed]},
    If[res === $Failed, 
      res,
      With[
        {newdat = ImageData[res] * (MRImageMax[img] - MRImageMin[img]) + MRImageMin[img],
         theOpts = Options[res]},
        MRImage3D[newdat, theOpts]]]],
  ImageDeconvolve[img, opts___] := With[
    {res = Check[ImageDeconvolve[Image3D[img], opts], $Failed]},
    If[res === $Failed, 
      res,
      With[
        {newdat = ImageData[res] * (MRImageMax[img] - MRImageMin[img]) + MRImageMin[img],
         theOpts = Options[res]},
        MRImage3D[newdat, theOpts]]]],
  ImageDistance[img, opts___] := ImageDistance[Image3D[img], opts],
  ImageForwardTransformation[img, opts___] := With[
    {res = Check[ImageForwardTransformation[Image3D[img], opts], $Failed]},
    If[res === $Failed, 
      res,
      With[
        {newdat = ImageData[res] * (MRImageMax[img] - MRImageMin[img]) + MRImageMin[img],
         theOpts = Options[res]},
        MRImage3D[newdat, theOpts]]]],
  ImagePad[img, opts___] := With[
    {res = Check[ImagePad[Image3D[img], opts], $Failed]},
    If[res === $Failed, 
      res,
      With[
        {newdat = ImageData[res] * (MRImageMax[img] - MRImageMin[img]) + MRImageMin[img],
         theOpts = Options[res]},
        MRImage3D[newdat, theOpts]]]],
  ImagePartition[img, opts___] := With[
    {res = Check[ImagePartition[Image3D[img], opts], $Failed]},
    If[res === $Failed, 
      res,
      With[
        {newdat = ImageData[res] * (MRImageMax[img] - MRImageMin[img]) + MRImageMin[img],
         theOpts = Options[res]},
        MRImage3D[newdat, theOpts]]]],
  ImagePerspectiveTransformation[img, opts___] := With[
    {res = Check[ImagePerspectiveTransformation[Image3D[img], opts], $Failed]},
    If[res === $Failed, 
      res,
      With[
        {newdat = ImageData[res] * (MRImageMax[img] - MRImageMin[img]) + MRImageMin[img],
         theOpts = Options[res]},
        MRImage3D[newdat, theOpts]]]],
  ImageReflect[img, opts___] := With[
    {res = Check[ImageReflect[Image3D[img], opts], $Failed]},
    If[res === $Failed, 
      res,
      With[
        {newdat = ImageData[res] * (MRImageMax[img] - MRImageMin[img]) + MRImageMin[img],
         theOpts = Options[res]},
        MRImage3D[newdat, theOpts]]]],
  ImageResize[img, opts___] := With[
    {res = Check[ImageResize[Image3D[img], opts], $Failed]},
    If[res === $Failed, 
      res,
      With[
        {newdat = ImageData[res] * (MRImageMax[img] - MRImageMin[img]) + MRImageMin[img],
         theOpts = Options[res]},
        MRImage3D[newdat, theOpts]]]],
  ImageRotate[img, opts___] := With[
    {res = Check[ImageRotate[Image3D[img], opts], $Failed]},
    If[res === $Failed, 
      res,
      With[
        {newdat = ImageData[res] * (MRImageMax[img] - MRImageMin[img]) + MRImageMin[img],
         theOpts = Options[res]},
        MRImage3D[newdat, theOpts]]]],
  ImageTake[img, opts___] := With[
    {res = Check[ImageTake[Image3D[img], opts], $Failed]},
    If[res === $Failed, 
      res,
      With[
        {newdat = ImageData[res] * (MRImageMax[img] - MRImageMin[img]) + MRImageMin[img],
         theOpts = Options[res]},
        MRImage3D[newdat, theOpts]]]],
  ImageTransformation[img, opts___] := With[
    {res = Check[ImageTransformation[Image3D[img], opts], $Failed]},
    If[res === $Failed, 
      res,
      With[
        {newdat = ImageData[res] * (MRImageMax[img] - MRImageMin[img]) + MRImageMin[img],
         theOpts = Options[res]},
        MRImage3D[newdat, theOpts]]]],
  ImageTrim[img, opts___] := With[
    {res = Check[ImageTrim[Image3D[img], opts], $Failed]},
    If[res === $Failed, 
      res,
      With[
        {newdat = ImageData[res] * (MRImageMax[img] - MRImageMin[img]) + MRImageMin[img],
         theOpts = Options[res]},
        MRImage3D[newdat, theOpts]]]],
  Sharpen[img, opts___] := With[
    {res = Check[Sharpen[Image3D[img], opts], $Failed]},
    If[res === $Failed, 
      res,
      With[
        {newdat = ImageData[res] * (MRImageMax[img] - MRImageMin[img]) + MRImageMin[img],
         theOpts = Options[res]},
        MRImage3D[newdat, theOpts]]]],

  (* Image Filters *)
  BandpassFilter[img, opts___] := With[
    {res = Check[BandpassFilter[Image3D[img], opts], $Failed]},
    If[res === $Failed, 
      res,
      With[
        {newdat = ImageData[res] * (MRImageMax[img] - MRImageMin[img]) + MRImageMin[img],
         theOpts = Options[res]},
        MRImage3D[newdat, theOpts]]]],
  BandstopFilter[img, opts___] := With[
    {res = Check[BandstopFilter[Image3D[img], opts], $Failed]},
    If[res === $Failed, 
      res,
      With[
        {newdat = ImageData[res] * (MRImageMax[img] - MRImageMin[img]) + MRImageMin[img],
         theOpts = Options[res]},
        MRImage3D[newdat, theOpts]]]],
  BilateralFilter[img, opts___] := With[
    {res = Check[BilateralFilter[Image3D[img], opts], $Failed]},
    If[res === $Failed, 
      res,
      With[
        {newdat = ImageData[res] * (MRImageMax[img] - MRImageMin[img]) + MRImageMin[img],
         theOpts = Options[res]},
        MRImage3D[newdat, theOpts]]]],
  CommonestFilter[img, opts___] := With[
    {res = Check[CommonestFilter[Image3D[img], opts], $Failed]},
    If[res === $Failed, 
      res,
      With[
        {newdat = ImageData[res] * (MRImageMax[img] - MRImageMin[img]) + MRImageMin[img],
         theOpts = Options[res]},
        MRImage3D[newdat, theOpts]]]],
  CurvatureFlowFilter[img, opts___] := With[
    {res = Check[CurvatureFlowFilter[Image3D[img], opts], $Failed]},
    If[res === $Failed, 
      res,
      With[
        {newdat = ImageData[res] * (MRImageMax[img] - MRImageMin[img]) + MRImageMin[img],
         theOpts = Options[res]},
        MRImage3D[newdat, theOpts]]]],
  DerivativeFilter[img, opts___] := With[
    {res = Check[DerivativeFilter[Image3D[img], opts], $Failed]},
    If[res === $Failed, 
      res,
      With[
        {newdat = ImageData[res] * (MRImageMax[img] - MRImageMin[img]) + MRImageMin[img],
         theOpts = Options[res]},
        MRImage3D[newdat, theOpts]]]],
  DifferentiatorFilter[img, opts___] := With[
    {res = Check[DifferentiatorFilter[Image3D[img], opts], $Failed]},
    If[res === $Failed, 
      res,
      With[
        {newdat = ImageData[res] * (MRImageMax[img] - MRImageMin[img]) + MRImageMin[img],
         theOpts = Options[res]},
        MRImage3D[newdat, theOpts]]]],
  EntropyFilter[img, opts___] := With[
    {res = Check[EntropyFilter[Image3D[img], opts], $Failed]},
    If[res === $Failed, 
      res,
      With[
        {newdat = ImageData[res] * (MRImageMax[img] - MRImageMin[img]) + MRImageMin[img],
         theOpts = Options[res]},
        MRImage3D[newdat, theOpts]]]],
  GaborFilter[img, opts___] := With[
    {res = Check[GaborFilter[Image3D[img], opts], $Failed]},
    If[res === $Failed, 
      res,
      With[
        {newdat = ImageData[res] * (MRImageMax[img] - MRImageMin[img]) + MRImageMin[img],
         theOpts = Options[res]},
        MRImage3D[newdat, theOpts]]]],
  GaussianFilter[img, opts___] := With[
    {res = Check[GaussianFilter[Image3D[img], opts], $Failed]},
    If[res === $Failed, 
      res,
      With[
        {newdat = ImageData[res] * (MRImageMax[img] - MRImageMin[img]) + MRImageMin[img],
         theOpts = Options[res]},
        MRImage3D[newdat, theOpts]]]],
  GeometricMeanFilter[img, opts___] := With[
    {res = Check[GeometricMeanFilter[Image3D[img], opts], $Failed]},
    If[res === $Failed, 
      res,
      With[
        {newdat = ImageData[res] * (MRImageMax[img] - MRImageMin[img]) + MRImageMin[img],
         theOpts = Options[res]},
        MRImage3D[newdat, theOpts]]]],
  GradientFilter[img, opts___] := With[
    {res = Check[GradientFilter[Image3D[img], opts], $Failed]},
    If[res === $Failed, 
      res,
      With[
        {newdat = ImageData[res] * (MRImageMax[img] - MRImageMin[img]) + MRImageMin[img],
         theOpts = Options[res]},
        MRImage3D[newdat, theOpts]]]],
  GradientOrientedFilter[img, opts___] := With[
    {res = Check[GradientOrientedFilter[Image3D[img], opts], $Failed]},
    If[res === $Failed, 
      res,
      With[
        {newdat = ImageData[res] * (MRImageMax[img] - MRImageMin[img]) + MRImageMin[img],
         theOpts = Options[res]},
        MRImage3D[newdat, theOpts]]]],
  HarmonicMeanFilter[img, opts___] := With[
    {res = Check[HarmonicMeanFilter[Image3D[img], opts], $Failed]},
    If[res === $Failed, 
      res,
      With[
        {newdat = ImageData[res] * (MRImageMax[img] - MRImageMin[img]) + MRImageMin[img],
         theOpts = Options[res]},
        MRImage3D[newdat, theOpts]]]],
  HighpassFilter[img, opts___] := With[
    {res = Check[HighpassFilter[Image3D[img], opts], $Failed]},
    If[res === $Failed, 
      res,
      With[
        {newdat = ImageData[res] * (MRImageMax[img] - MRImageMin[img]) + MRImageMin[img],
         theOpts = Options[res]},
        MRImage3D[newdat, theOpts]]]],
  HilbertFilter[img, opts___] := With[
    {res = Check[HilbertFilter[Image3D[img], opts], $Failed]},
    If[res === $Failed, 
      res,
      With[
        {newdat = ImageData[res] * (MRImageMax[img] - MRImageMin[img]) + MRImageMin[img],
         theOpts = Options[res]},
        MRImage3D[newdat, theOpts]]]],
  KuwaharaFilter[img, opts___] := With[
    {res = Check[KuwaharaFilter[Image3D[img], opts], $Failed]},
    If[res === $Failed, 
      res,
      With[
        {newdat = ImageData[res] * (MRImageMax[img] - MRImageMin[img]) + MRImageMin[img],
         theOpts = Options[res]},
        MRImage3D[newdat, theOpts]]]],
  LaplacianFilter[img, opts___] := With[
    {res = Check[LaplacianFilter[Image3D[img], opts], $Failed]},
    If[res === $Failed, 
      res,
      With[
        {newdat = ImageData[res] * (MRImageMax[img] - MRImageMin[img]) + MRImageMin[img],
         theOpts = Options[res]},
        MRImage3D[newdat, theOpts]]]],
  LaplacianGaussianFilter[img, opts___] := With[
    {res = Check[LaplacianGaussianFilter[Image3D[img], opts], $Failed]},
    If[res === $Failed, 
      res,
      With[
        {newdat = ImageData[res] * (MRImageMax[img] - MRImageMin[img]) + MRImageMin[img],
         theOpts = Options[res]},
        MRImage3D[newdat, theOpts]]]],
  LowpassFilter[img, opts___] := With[
    {res = Check[LowpassFilter[Image3D[img], opts], $Failed]},
    If[res === $Failed, 
      res,
      With[
        {newdat = ImageData[res] * (MRImageMax[img] - MRImageMin[img]) + MRImageMin[img],
         theOpts = Options[res]},
        MRImage3D[newdat, theOpts]]]],
  MaxFilter[img, opts___] := With[
    {res = Check[MaxFilter[Image3D[img], opts], $Failed]},
    If[res === $Failed, 
      res,
      With[
        {newdat = ImageData[res] * (MRImageMax[img] - MRImageMin[img]) + MRImageMin[img],
         theOpts = Options[res]},
        MRImage3D[newdat, theOpts]]]],
  MeanFilter[img, opts___] := With[
    {res = Check[MeanFilter[Image3D[img], opts], $Failed]},
    If[res === $Failed, 
      res,
      With[
        {newdat = ImageData[res] * (MRImageMax[img] - MRImageMin[img]) + MRImageMin[img],
         theOpts = Options[res]},
        MRImage3D[newdat, theOpts]]]],
  MeanShiftFilter[img, opts___] := With[
    {res = Check[MeanShiftFilter[Image3D[img], opts], $Failed]},
    If[res === $Failed, 
      res,
      With[
        {newdat = ImageData[res] * (MRImageMax[img] - MRImageMin[img]) + MRImageMin[img],
         theOpts = Options[res]},
        MRImage3D[newdat, theOpts]]]],
  MedianFilter[img, opts___] := With[
    {res = Check[MedianFilter[Image3D[img], opts], $Failed]},
    If[res === $Failed, 
      res,
      With[
        {newdat = ImageData[res] * (MRImageMax[img] - MRImageMin[img]) + MRImageMin[img],
         theOpts = Options[res]},
        MRImage3D[newdat, theOpts]]]],
  MinFilter[img, opts___] := With[
    {res = Check[MinFilter[Image3D[img], opts], $Failed]},
    If[res === $Failed, 
      res,
      With[
        {newdat = ImageData[res] * (MRImageMax[img] - MRImageMin[img]) + MRImageMin[img],
         theOpts = Options[res]},
        MRImage3D[newdat, theOpts]]]],
  PeronaMalikFilter[img, opts___] := With[
    {res = Check[PeronaMalikFilter[Image3D[img], opts], $Failed]},
    If[res === $Failed, 
      res,
      With[
        {newdat = ImageData[res] * (MRImageMax[img] - MRImageMin[img]) + MRImageMin[img],
         theOpts = Options[res]},
        MRImage3D[newdat, theOpts]]]],
  RangeFilter[img, opts___] := With[
    {res = Check[RangeFilter[Image3D[img], opts], $Failed]},
    If[res === $Failed, 
      res,
      With[
        {newdat = ImageData[res] * (MRImageMax[img] - MRImageMin[img]) + MRImageMin[img],
         theOpts = Options[res]},
        MRImage3D[newdat, theOpts]]]],
  RidgeFilter[img, opts___] := With[
    {res = Check[RidgeFilter[Image3D[img], opts], $Failed]},
    If[res === $Failed, 
      res,
      With[
        {newdat = ImageData[res] * (MRImageMax[img] - MRImageMin[img]) + MRImageMin[img],
         theOpts = Options[res]},
        MRImage3D[newdat, theOpts]]]],
  StandardDeviationFilter[img, opts___] := With[
    {res = Check[StandardDeviationFilter[Image3D[img], opts], $Failed]},
    If[res === $Failed, 
      res,
      With[
        {newdat = ImageData[res] * (MRImageMax[img] - MRImageMin[img]) + MRImageMin[img],
         theOpts = Options[res]},
        MRImage3D[newdat, theOpts]]]],
  TotalVariationFilter[img, opts___] := With[
    {res = Check[TotalVariationFilter[Image3D[img], opts], $Failed]},
    If[res === $Failed, 
      res,
      With[
        {newdat = ImageData[res] * (MRImageMax[img] - MRImageMin[img]) + MRImageMin[img],
         theOpts = Options[res]},
        MRImage3D[newdat, theOpts]]]],
  WienerFilter[img, opts___] := With[
    {res = Check[WienerFilter[Image3D[img], opts], $Failed]},
    If[res === $Failed, 
      res,
      With[
        {newdat = ImageData[res] * (MRImageMax[img] - MRImageMin[img]) + MRImageMin[img],
         theOpts = Options[res]},
        MRImage3D[newdat, theOpts]]]],

  (* Morphological Operations *)
  Binarize[img, opts___] := Binarize[Image3D[img], opts],
  BottomHatTransform[img, opts___] := BottomHatTransform[Image3D[img], opts],
  Closing[img, opts___] := With[
    {res = Check[Closing[Image3D[img], opts], $Failed]},
    If[res === $Failed, 
      res,
      With[
        {newdat = ImageData[res] * (MRImageMax[img] - MRImageMin[img]) + MRImageMin[img],
         theOpts = Options[res]},
        MRImage3D[newdat, theOpts]]]],
  ColorNegate[img, opts___] := ColorNegate[Image3D[img], opts],
  ComponentMeasurements[img, opts___] := ComponentMeasurements[Image3D[img], opts],
  Colorize[img, opts___] := Colorize[Image3D[img], opts],
  DeleteBorderComponents[img, opts___] := With[
    {res = Check[DeleteBorderComponents[Image3D[img], opts], $Failed]},
    If[res === $Failed, 
      res,
      With[
        {newdat = ImageData[res] * (MRImageMax[img] - MRImageMin[img]) + MRImageMin[img],
         theOpts = Options[res]},
        MRImage3D[newdat, theOpts]]]],
  DeleteSmallComponents[img, opts___] := With[
    {res = Check[DeleteSmallComponents[Image3D[img], opts], $Failed]},
    If[res === $Failed, 
      res,
      With[
        {newdat = ImageData[res] * (MRImageMax[img] - MRImageMin[img]) + MRImageMin[img],
         theOpts = Options[res]},
        MRImage3D[newdat, theOpts]]]],
  Dilation[img, opts___] := With[
    {res = Check[Dilation[Image3D[img], opts], $Failed]},
    If[res === $Failed, 
      res,
      With[
        {newdat = ImageData[res] * (MRImageMax[img] - MRImageMin[img]) + MRImageMin[img],
         theOpts = Options[res]},
        MRImage3D[newdat, theOpts]]]],
  DistanceTranform[img, opts___] := DistanceTransgorm[Image3D[img], opts],
  Erosion[img, opts___] := With[
    {res = Check[Erosion[Image3D[img], opts], $Failed]},
    If[res === $Failed, 
      res,
      With[
        {newdat = ImageData[res] * (MRImageMax[img] - MRImageMin[img]) + MRImageMin[img],
         theOpts = Options[res]},
        MRImage3D[newdat, theOpts]]]],
  FillingTransform[img, opts___] := With[
    {res = Check[FillingTransform[Image3D[img], opts], $Failed]},
    If[res === $Failed, 
      res,
      With[
        {newdat = ImageData[res] * (MRImageMax[img] - MRImageMin[img]) + MRImageMin[img],
         theOpts = Options[res]},
        MRImage3D[newdat, theOpts]]]],
  HighlightImage[img, opts___] := HighlightImage[Image3D[img], opts],
  HitMissTransform[img, opts___] := With[
    {res = Check[HitMissTransform[Image3D[img], opts], $Failed]},
    If[res === $Failed, 
      res,
      With[
        {newdat = ImageData[res] * (MRImageMax[img] - MRImageMin[img]) + MRImageMin[img],
         theOpts = Options[res]},
        MRImage3D[newdat, theOpts]]]],
  InverseDistanceTransform[img, opts___] := With[
    {res = Check[InverseDistanceTransform[Image3D[img], opts], $Failed]},
    If[res === $Failed, 
      res,
      With[
        {newdat = ImageData[res] * (MRImageMax[img] - MRImageMin[img]) + MRImageMin[img],
         theOpts = Options[res]},
        MRImage3D[newdat, theOpts]]]],
  MaxDetect[img, opts___] := MaxDetect[Image3D[img], opts],
  MinDetect[img, opts___] := MinDetect[Image3D[img], opts],
  MorphologicalBinarize[img, opts___] := MorphologicalBinarize[Image3D[img], opts],
  MorphologicalBranchPoints[img, opts___] := MorphologicalBranchPoints[Image3D[img], opts],
  MorphologicalComponents[img, opts___] := MorphologicalComponents[Image3D[img], opts],
  MorphologicalEulerNumber[img, opts___] := MorphologicalEulerNumber[Image3D[img], opts],
  MorphologicalGraph[img, opts___] := MorphologicalGraph[Image3D[img], opts],
  MorphologicalPerimeter[img, opts___] := MorphologicalPerimeter[Image3D[img], opts],
  MorphologicalTransform[img, opts___] := MorphologicalTransform[Image3D[img], opts],
  Opening[img, opts___] := With[
    {res = Check[Opening[Image3D[img], opts], $Failed]},
    If[res === $Failed, 
      res,
      With[
        {newdat = ImageData[res] * (MRImageMax[img] - MRImageMin[img]) + MRImageMin[img],
         theOpts = Options[res]},
        MRImage3D[newdat, theOpts]]]],
  Pruning[img, opts___] := With[
    {res = Check[Pruning[Image3D[img], opts], $Failed]},
    If[res === $Failed, 
      res,
      With[
        {newdat = ImageData[res] * (MRImageMax[img] - MRImageMin[img]) + MRImageMin[img],
         theOpts = Options[res]},
        MRImage3D[newdat, theOpts]]]],
  SelectComponents[img, opts___] := SelectComponents[Image3D[img], opts],
  SkeletonTranform[img, opts___] := SkeletonTransgorm[Image3D[img], opts],
  Thinning[img, opts___] := Thinning[Image3D[img], opts],
  TopHatTransform[img, opts___] := TopHatTransform[Image3D[img], opts]];

Options[MRImage3D] = Join[
  {CoordinateToVoxelMatrix -> Automatic,
   VoxelToCoordinateMatrix -> Automatic},
  Replace[
    Options[Image3D],
    {(ColorFunction -> _) -> (ColorFunction -> "XRay"),
     (Boxed -> _) -> (Boxed -> False),
     (Background -> _) -> (Background -> Black)},
    {1}]];
DefineImmutable[
  MRImage3D[data_, OptionsPattern[]] :> img,
  Evaluate[
    Join[
      $MRImageSharedMethods,
      Hold[
        (* Retreive (or edit) the raw image data *)
        ImageData[img] = data,
        (* Retreive (or edit) the raw image options *)
        Options[img] = Map[(# -> OptionValue[#])&, Options[MRImage3D][[All, 1]]],
        
        (* The image as it is represented internally *)
        Image3D[img] -> With[
          {dat = ImageData[img],
           opts = Options[img]},
          If[!ArrayQ[dat, 3|4, NumericQ], 
            Message[MRImage3D::badarg, "ImageData must be a 3 or 4D numeric array"],
            Image3D[
              (dat - MRImageMin[img]) / (MRImageMax[img] - MRImageMin[img]),
              (*"Real32",*)
              Sequence@@FilterRules[opts, Options[Image3D][[All,1]]]]]],
        
        (* This is an MRImage... Checks of data quality should go here as well. *)
        MRImageQ[img] -> True,

        (* #todo this should eventually be fixed to give back MRImage's instead of Image's*)
        Image3DSlices[img, opts___] := With[
          {max = MRImageMax[img],
           min = MRImageMin[img]},
          Map[
            Function[
              MRImage[
                ImageData[#] * (max - min) + min,
                With[
                  {o = If[Count[#, ("MetaInformation"|MetaInformation) -> _, {1}] == 0,
                     Append[#, "MetaInformation" -> {}], #]&[Options[#]]},
                  Sequence @@ Replace[
                    o,
                    Rule[("MetaInformation"|MetaInformation), l_] :> Rule[
                      "MetaInformation",
                      Append[l, "SourceImage" -> img]],
                    {1}]]]],
          Image3DSlices[Image3D[img], opts]]],
        
        (* We want to support the coordinate-to-voxel matrices, if provided *)
        CoordinateToVoxelMatrix[img] := With[
          {opt = Replace[
             Replace[CoordinateToVoxelMatrix, Options[img]],
             Automatic :> With[
               {meta = Replace[MetaInformation, Options[img]]},
               With[
                 {mtx = If[ListQ[meta], "VOXToRASMatrix" /. meta, "None"]},
                 If[StringQ[mtx],
                   {{1,0,0,0},{0,1,0,0},{0,0,1,0}},
                   Plus[
                     ReplacePart[
                       If[Length[mtx] == 4, Most[mtx], mtx],
                       {_,4} -> 0],
                     Transpose @ Append[
                       ConstantArray[0, {3,3}],
                       0.5 * ImageDimensions[img]]]]]]]},
          Which[
            opt === CoordinateToVoxelMatrix, None,
            !ListQ[opt], $Failed,
            !MatchQ[Dimensions[opt], {3|4, 4}], $Failed,
            Length[opt] == 4, Most[opt],
            True, opt]],
        VoxelToCoordinateMatrix[img] := With[
          {opt = Replace[VoxelToCoordinateMatrix, Options[img]]},
          Which[
            opt === VoxelToCoordinateMatrix, None,
            opt === Automatic, With[
              {c2v = CoordinateToVoxelMatrix[img]},
              Which[
                c2v === None, None,
                c2v === $Failed, $Failed,
                True, Most[Inverse[Append[c2v, {0,0,0,1}]]]]],
            !ListQ[opt], $Failed,
            !MatchQ[Dimensions[opt], {3|4, 4}], $Failed,
            True, opt]]]]],
  SetSafe ->True,
  Symbol -> MRImage3D];

(* #MRImage ***************************************************************************************)
Options[MRImage] = Join[
  {},
  Replace[
    Options[Image],
    {(ColorFunction -> _) -> (ColorFunction -> "XRay"),
     (Boxed -> _) -> (Boxed -> False),
     (Background -> _) -> (Background -> Black)},
    {1}]];
DefineImmutable[
  MRImage[data_, OptionsPattern[]] :> img,
  Evaluate[
    Join[
      ReplaceAll[$MRImageSharedMethods, {Image3D -> Image, MRImage3D -> MRImage}],
      Hold[
        (* Retreive (or edit) the raw image data *)
        ImageData[img] = data,
        (* Retreive (or edit) the raw image options *)
        Options[img] = Map[(# -> OptionValue[#])&, Options[MRImage][[All, 1]]],
        
        (* The image as it is represented internally *)
        Image[img] -> With[
          {dat = ImageData[img],
           opts = Options[img]},
          If[!ArrayQ[dat, 2|3, NumericQ], 
            Message[MRImage::badarg, "ImageData must be a 2 or 3D numeric array"],
            Image[
              (dat - MRImageMin[img]) / (MRImageMax[img] - MRImageMin[img]),
              (*"Real32",*)
              Sequence@@FilterRules[opts, Options[Image][[All,1]]]]]],
        
        (* This is an MRImage... Checks of data quality should go here as well. *)
        MRImageSliceQ[img] -> True]]],
  SetSafe ->True,
  Symbol -> MRImage];

(* #MakeBoxes allows us to display the image as a 3d image *)
MakeBoxes[img:MRImage3D[_,Except[_Rule],__], form_] := MakeBoxes[#, form]&[Image3D[img]];
MakeBoxes[img:MRImage[_,Except[_Rule],__], form_] := MakeBoxes[#, form]&[Image[img]];

(* We want to define the #MRImageObjectQ and related functions now... *)
MRImageQ[x_] := False;
MRImageSliceQ[x_] := False;
MRImageObjectQ[img_] := Or[MRImageQ[img], MRImageSliceQ[img]];
Protect[MRImageObjectQ, MRImageQ, MRImageSliceQ];

(* For all of the functions below, we want to write them only once, for the 3D case, but to define
 * them twice, also for the 2D case. Their code is almost identical, however, so a hack is employed
 * to prevent code duplication.
 *)
(ReleaseHold[#]; ReleaseHold[# //. {Image3D -> Image, MRImage3D -> MRImage}];)& @ Hold[

   (* #ImageAlign is a slightly special case *)
   Unprotect[ImageAlign];
   ImageAlign[a_, b_MRImage3D, opts___] := With[
     {res = Check[ImageAlign[a, Image3D[b], opts], $Failed]},
     If[res === $Failed, 
       res,
       With[
         {res = Check[ImageCrop[Image3D[img], opts], $Failed]},
         If[res === $Failed, 
           res,
           With[
             {newdat = ImageData[res] * (MRImageMax[b] - MRImageMin[b]) + MRImageMin[b],
              theOpts = Options[res]},
             MRImage3D[newdat, theOpts]]]]]];
   Protect[ImageAlign];

   (* #ImageCorrespondingPoints also is a special case *)
   Unprotect[ImageCorrespondingPoints];
   ImageCorrespondingPoints[a_, b_MRImage3D, opts___] := ImageCorrespondingPoints[
     a,
     Image3D[b],
     opts];
   Protect[ImageCorrespondingPoints];
   
   (* #ImageDistance is also a special case *)
   Unprotect[ImageDistance];
   ImageDistance[a_, b_MRImage3D, opts___] := ImageDistance[a, Image3D[b], opts];
   Protect[ImageDistance];

   (* As is #ImageFilter *)
   Unprotect[ImageFilter];
   (* I have no idea why this is necessary, but it seems ImageFilter reprotects itself *)
   ImageFilter;
   Unprotect[ImageFilter];
   ImageFilter[f_, img_MRImage3D, opts___] := With[
     {res = Check[ImageFilter[f, Image3D[img], opts], $Failed]},
     If[res === $Failed, 
       res,
       With[
         {newdat = ImageData[res] * (MRImageMax[img] - MRImageMin[img]) + MRImageMin[img],
          theOpts = Options[res]},
         MRImage3D[newdat, theOpts]]]];
   Protect[ImageFilter];];

(* #MeshVertexToVoxelIndex **********************************************************************)
MeshVertexToVoxelIndex[p_List, voxdim:{_NumericQ, _NumericQ, _NumericQ}] := Which[
  Length[p] == 0, p,
  ListQ[First[p]], Map[MeshVertexToVoxelIndex[#, voxdim]&, p],
  True, List[
    0.5*voxdim[[1]] - p[[1]] + 1.0, 
    0.5*voxdim[[2]] - p[[3]] + 1.0, 
    0.5*voxdim[[3]] + p[[2]] + 1.0]];
MeshVertexToVoxelIndex[p_List, vox_List /; ArrayQ[vox, 3]] := MeshVertexToVoxelIndex[
  p,
  Dimensions[vox]];
MeshVertexToVoxelIndex[p_List] := MeshVertexToVoxelIndex[p, {256.0, 256.0, 256.0}];
Protect[MeshVertexToVoxelIndex];

(* #VoxelToCoordinate *****************************************************************************)
VoxelToCoordinate[vol_?MRImageQ, idx:{_?NumericQ, _?NumericQ, _?NumericQ}] := With[
  {dims = ImageDimensions[vol],
   matrix = VoxelToCoordinateMatrix[vol]},
  If[And@@MapThread[TrueQ[0 < #1 <= #2]&, {idx, dims}],
    (*Dot[matrix, Append[(idx - 1.0) * (dims / (dims - 1.0)), 1.0]],*)
    Dot[matrix, Append[{idx[[1]], dims[[2]] - idx[[2]] + 1, dims[[3]] - idx[[3]] + 1} - 0.5, 1.0]],
    (Message[VoxelToCoordinate::badarg, "index is not a valid voxel"];
     $Failed)]];
VoxelToCoordinate[vol_?MRImageQ, idcs:{{_?NumericQ, _?NumericQ, _?NumericQ}..}] := With[
  {dims = ImageDimensions[vol],
   matrix = VoxelToCoordinateMatrix[vol]},
  Map[
    Function[
      If[And@@MapThread[TrueQ[0 < #1 <= #2]&, {#, dims}],
        (*Dot[matrix, Append[(# - 1.0) * (dims / (dims - 1.0)), 1.0]],*)
        Dot[matrix, Append[{#[[1]], dims[[2]] - #[[2]] + 1, dims[[3]] - #[[3]] + 1} - 0.5, 1.0]],
        (Message[VoxelToCoordinate::badarg, "index is not a valid voxel"];
         $Failed)]],
    idcs]];
Protect[VoxelToCoordinate];

(* #CoordinateToVoxel *****************************************************************************)
CoordinateToVoxel[vol_?MRImageQ, coord:{_?NumericQ, _?NumericQ, _?NumericQ}] := With[
  {dims = ImageDimensions[vol]},
  With[
    {idx0 = Dot[CoordinateToVoxelMatrix[vol], Append[coord, 1.0]] + 0.5},
    With[
      {idx = {idx0[[1]], dims[[2]] - idx0[[2]] + 1, dims[[3]] - idx0[[3]] + 1}},
      If[!And@@MapThread[TrueQ[0 < #1 <= #2]&, {idx, dims}],
        Message[CoordinateToVoxel::badarg, "coordinate is outside voxel range"]];
      idx]]];
CoordinateToVoxel[vol_?MRImageQ, coords:{{_?NumericQ, _?NumericQ, _?NumericQ}..}] := With[
  {dims = ImageDimensions[vol],
   matrix = CoordinateToVoxelMatrix[vol]},
  With[
    {cTr0 = Dot[matrix, Append[Transpose@coords, Table[1.0, {Length@coords}]]] + 0.5},
    With[
      {cTr = {cTr0[[1]], dims[[2]] - cTr0[[2]] + 1, dims[[3]] - cTr0[[3]] + 1}},
      If[0 < Total[MapThread[Count[#1, x_ /; x <= 0 || x > #2] &, {cTr, dims}]],
        (Message[CoordinateToVoxel::badarg, "coordinate is outside voxel range"]; $Failed),
        Transpose @ cTr]]]];
Protect[CoordinateToVoxel];


End[];
EndPackage[];

