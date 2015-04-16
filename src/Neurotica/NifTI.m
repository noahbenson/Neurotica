(* NifTI.m
 *
 * The Neurotica`NifTI namespace contains functions for reading and interpreting the NIH standard
 * neuroimaging file formats, including NifTI, GifTI, and SifTI.
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
BeginPackage[
  "Neurotica`NifTI`",
  {"Neurotica`Global`","Neurotica`Util`","Neurotica`Mesh`", "Neurotica`MRImage`", "JLink`"}];
Unprotect["Neurotica`NifTI`*", "Neurotica`NifTI`Private`*"];
ClearAll[ "Neurotica`NifTI`*", "Neurotica`NifTI`Private`*"];

ImportNifTI::usage = "ImportNifTI[source, options...] is equivalent to Import[source, \"NIFTI\", optios...].";
ImportNifTI::badarg = "Bad argument given to ImportNifTI: `1`";
ImportNifTI::badfmt = "Bad NifTI file format: `1`";

ImportGifTI::usage = "ImportGifTI[source, options...] is equivalent to Import[source, \"GIFTI\", optios...].";
ImportGifTI::badarg = "Bad argument given to ImportGifTI: `1`";
ImportGifTI::badfmt = "Bad GifTI file format: `1`";

Base64ZLibDecode::usage = "Base64ZLibDecode[string] yields the sequence of bytes (suitable for calls to FromCharacterCode and via that ImportString) corresponding to the base64-encoded zlib-zipped string. This function is intended to fill a hole in the GifTI file format, which specifies that GzipBase64 data should be in gzip format but which is actually stored in zlib stream format, which Mathematica does not support.";
Base64ZLibEncode::usage = "Base64ZLibEncode[data] yields the Base64 string encoding of the zlib-compressed sequence of bytes given in data. See also Base64ZLibDecode.";

Begin["`Private`"];

(* Base64ZLipDecode and Encode require Java to be installed *)
InstallJava[];

(* #Base64ZLibDecode ******************************************************************************)
Base64ZLibDecode[string_String] := With[
  {x = Apply[
     Join,
     First @ Last @ Reap @ JavaBlock[
       With[
         {inflater = JavaNew[
            "java.util.zip.InflaterInputStream",
            JavaNew[
              "java.io.ByteArrayInputStream",
              ImportString[string, {"Base64","Binary"}]]],
          ar = JavaNew["[B", 1024]},
         While[
           inflater@available[] != 0,
           With[
             {k = inflater@read[ar, 0, 1024]},
             If[k > 0, Sow[JavaObjectToExpression[ar][[1 ;; k]]]]]]]]]},
  (* Java bytes can be negative, but we need only positive for FromCharacterCode to work *)
  Abs[Sign[x]]*((-Sign[x] + 1)/2 * (x + 256) + (Sign[x] + 1)/2 * x)];
Protect[Base64ZLibDecode];

(* #Base64ZLibEncode ******************************************************************************)
Base64ZLibEncode[bytes_List] := ExportString[
  With[
    {x = Apply[
       Join,
       First@Last@Reap@JavaBlock[
         With[
           {deflater = JavaNew[
              "java.util.zip.DeflaterInputStream",
              JavaNew[
                "java.io.ByteArrayInputStream",
                q]],
            ar = JavaNew["[B", 1024]},
           While[
             deflater@available[] != 0,
             With[
               {k = deflater@read[ar, 0, 1024]},
               If[k > 0,
                 Sow[JavaObjectToExpression[ar][[1 ;; k]]]]]]]]]},
    (* Java bytes can be negative, but we need only positive for FromCharacterCode *)
    FromCharacterCode[ -(Sign[x] - 1)/2*(x + 256) + (Sign[x] + 1)/2*d ]],
  "Base64"];
Protect[Base64ZLibDecode];


(* #NifTI File Format *****************************************************************************)

(* ============================================================================================== *)
(* NifTI specifies several possible datatypes; here we paste the NifTI header file information
   about datatypes:

   #define DT_NONE                    0
   #define DT_UNKNOWN                 0     /* what it says, dude           */
   #define DT_BINARY                  1     /* binary (1 bit/voxel)         */
   #define DT_UNSIGNED_CHAR           2     /* unsigned char (8 bits/voxel) */
   #define DT_SIGNED_SHORT            4     /* signed short (16 bits/voxel) */
   #define DT_SIGNED_INT              8     /* signed int (32 bits/voxel)   */
   #define DT_FLOAT                  16     /* float (32 bits/voxel)        */
   #define DT_COMPLEX                32     /* complex (64 bits/voxel)      */
   #define DT_DOUBLE                 64     /* double (64 bits/voxel)       */
   #define DT_RGB                   128     /* RGB triple (24 bits/voxel)   */
   #define DT_ALL                   255     /* not very useful (?)          */

   #define DT_UINT8                   2
   #define DT_INT16                   4
   #define DT_INT32                   8
   #define DT_FLOAT32                16
   #define DT_COMPLEX64              32
   #define DT_FLOAT64                64
   #define DT_RGB24                 128
   
   #define DT_INT8                  256     /* signed char (8 bits)         */
   #define DT_UINT16                512     /* unsigned short (16 bits)     */
   #define DT_UINT32                768     /* unsigned int (32 bits)       */
   #define DT_INT64                1024     /* long long (64 bits)          */
   #define DT_UINT64               1280     /* unsigned long long (64 bits) */
   #define DT_FLOAT128             1536     /* long double (128 bits)       */
   #define DT_COMPLEX128           1792     /* double pair (128 bits)       */
   #define DT_COMPLEX256           2048     /* long double pair (256 bits)  */
   #define DT_RGBA32               2304     /* 4 byte RGBA (32 bits/voxel)  */                    *)
(* ============================================================================================== *)
$NifTIDatatypes = {
  {1,    "Bit"},
  {2,    "UnsignedInteger8"},
  {4,    "Integer16"},
  {8,    "Integer32"},
  {16,   "Real32"},
  {32,   "Complex64"},
  {64,   "Real64"},
  {128,  "RGB24"},
  {256,  "Integer8"},
  {512,  "UnsignedInteger16"},
  {768,  "UnsignedInteger32"},
  {1024, "Integer64"},
  {1280, "UnsignedInteger64"},
  {1536, "Real128"},
  {1792, "Complex128"},
  {2038, "Complex256"},
  {2304, "RGB32"}};
NifTIDatatypeTranslate[type_Integer] := Replace[
  type,
  Append[
    (Rule@@#)& /@ $NifTIDatatypes,
    _ :> Message[ImportNifTI::badfmt, "Could not recognize datatype id: "<>ToString[type]]]];
NifTIDatatypeUntranslate[type_String] := Replace[
  type,
  Append[
    (Rule@@Reverse[#])& /@ $NifTIDatatypes,
    _ :> Message[ImportNifTI::badfmt, "Could not recognize datatype name: "<>type]]];
NifTIDataToBinaryType[type_String] := (type /. {"RGB24" -> "Integer24", "RGB32" -> "Integer32"});
NifTIColorTranslate[type_String, data_List] := If[
  StringLength[type] < 3 || StringTake[type, 3] != "RGB",
  data,
  With[
    {bits = ToExpression @ StringTake[type, {4, -1}]},
    IntegerDigits[BitAnd[data, 2^bits - 1], 256]]];
NifTIColorTranslate[type_Integer, data_List] := NifTIColorTranslate[
  NifTIDatatypeTranslate[type],
  data];
NifTIColorUntranslate[type_String, data_List] := If[
  StringLength[type] < 3 || StringTake[type, 3] != "RGB",
  data,
  With[
    {bits = ToExpression @ StringTake[type, {4, -1}]},
    With[
      {u = If[bits == 24, {256^2, 256, 1}, {256^3, 256^2, 256, 1}]},
      Map[Dot[#,u]&, data, {-1}]]]];
NifTIColorUntranslate[type_Integer, data_List] := NifTIColorUntranslate[
  NifTIDatatypeTranslate[type],
  data];
Protect[NifTIDatatypeTranslate, NifTIDatatypeUntranslate,
        NifTIDataToBinaryType,
        NifTIColorTranslate, NifTIColorUntranslate];

(* ============================================================================================== *)
(* Here we include the NifTI specs for unit codes
                                  /*! NIFTI code for unspecified units. */
   #define NIFTI_UNITS_UNKNOWN 0

                                  /** Space codes are multiples of 1. **/
                                  /*! NIFTI code for meters. */
   #define NIFTI_UNITS_METER   1
                                  /*! NIFTI code for millimeters. */
   #define NIFTI_UNITS_MM      2
                                  /*! NIFTI code for micrometers. */
   #define NIFTI_UNITS_MICRON  3
   
                                  /** Time codes are multiples of 8. **/
                                  /*! NIFTI code for seconds. */
   #define NIFTI_UNITS_SEC     8
                                  /*! NIFTI code for milliseconds. */
   #define NIFTI_UNITS_MSEC   16
                                  /*! NIFTI code for microseconds. */
   #define NIFTI_UNITS_USEC   24
   
                                  /*** These units are for spectral data: ***/
                                  /*! NIFTI code for Hertz. */
   #define NIFTI_UNITS_HZ     32
                                  /*! NIFTI code for ppm. */
   #define NIFTI_UNITS_PPM    40
                                  /*! NIFTI code for radians per second. */
   #define NIFTI_UNITS_RADS   48                                                                  *)
(* ============================================================================================== *)
NifTIUnitTranslate[unit_Integer] := {
  Switch[
    BitAnd[unit, 7],
    1, "Meters",
    2, "Millimeters",
    3, "Micrometers",
    _, None],
  Switch[
    BitAnd[unit, BitNot[7]],
    8, "Seconds",
    16, "Milliseconds",
    24, "Microseconds",
    32, "Hertz",
    40, "PartsPerMillion",
    48, "RadiansPerSecond",
    _, None]};
NifTIUnitTranslate[unit_] := BitOr[
  Switch[
    unit[[1]]
    None, 0,
    "Meters", 1,
    "Millimeters", 2,
    "Micrometers", 3,
    _, 0],
  Switch[
    unit[[2]],
    "Seconds", 8,
    "Milliseconds", 16,
    "Microseconds", 24,
    "Hertz", 32,
    "PartsPerMillion", 40,
    "RadiansPerSecond", 48]];
Protect[NifTIUnitTranslate, NifTIUnitUntranslate];

(* NifTI Header Extension codes *******************************************************************)
$NifTIExtensionCodeCifTI = 32;
Protect[$NifTIExtensionCodeCifTI];


(* ============================================================================================== *)
(* Here we include the NifTI file spec header:
   int   sizeof_hdr;    /*!< MUST be 348           */  /* int sizeof_hdr;      */    0-3,
   char  data_type[10]; /*!< ++UNUSED++            */  /* char data_type[10];  */    4-13,
   char  db_name[18];   /*!< ++UNUSED++            */  /* char db_name[18];    */    14-31,
   int   extents;       /*!< ++UNUSED++            */  /* int extents;         */    32-35,
   short session_error; /*!< ++UNUSED++            */  /* short session_error; */    36-37,
   char  regular;       /*!< ++UNUSED++            */  /* char regular;        */    38,
   char  dim_info;      /*!< MRI slice ordering.   */  /* char hkey_un0;       */    39,
 
   short dim[8];        /*!< Data array dimensions.*/  /* short dim[8];        */    40
   float intent_p1 ;    /*!< 1st intent parameter. */  /* short unused8;       */
                                                       /* short unused9;       */
   float intent_p2 ;    /*!< 2nd intent parameter. */  /* short unused10;      */
                                                       /* short unused11;      */
   float intent_p3 ;    /*!< 3rd intent parameter. */  /* short unused12;      */
                                                       /* short unused13;      */
   short intent_code ;  /*!< NIFTI_INTENT_* code.  */  /* short unused14;      */
   short datatype;      /*!< Defines data type!    */  /* short datatype;      */
   short bitpix;        /*!< Number bits/voxel.    */  /* short bitpix;        */
   short slice_start;   /*!< First slice index.    */  /* short dim_un0;       */
   float pixdim[8];     /*!< Grid spacings.        */  /* float pixdim[8];     */
   float vox_offset;    /*!< Offset into .nii file */  /* float vox_offset;    */
   float scl_slope ;    /*!< Data scaling: slope.  */  /* float funused1;      */
   float scl_inter ;    /*!< Data scaling: offset. */  /* float funused2;      */
   short slice_end;     /*!< Last slice index.     */  /* float funused3;      */
   char  slice_code ;   /*!< Slice timing order.   */
   char  xyzt_units ;   /*!< Units of pixdim[1..4] */
   float cal_max;       /*!< Max display intensity */  /* float cal_max;       */
   float cal_min;       /*!< Min display intensity */  /* float cal_min;       */
   float slice_duration;/*!< Time for 1 slice.     */  /* float compressed;    */
   float toffset;       /*!< Time axis shift.      */  /* float verified;      */
   int   glmax;         /*!< ++UNUSED++            */  /* int glmax;           */
   int   glmin;         /*!< ++UNUSED++            */  /* int glmin;           */
  
   char  descrip[80];   /*!< any text you like.    */  /* char descrip[80];    */
   char  aux_file[24];  /*!< auxiliary filename.   */  /* char aux_file[24];   */
  
   short qform_code ;   /*!< NIFTI_XFORM_* code.   */  /*-- all ANALYZE 7.5 ---*/
   short sform_code ;   /*!< NIFTI_XFORM_* code.   */  /*   fields below here  */
                                                       /*   are replaced       */
   float quatern_b ;    /*!< Quaternion b param.   */
   float quatern_c ;    /*!< Quaternion c param.   */
   float quatern_d ;    /*!< Quaternion d param.   */
   float qoffset_x ;    /*!< Quaternion x shift.   */
   float qoffset_y ;    /*!< Quaternion y shift.   */
   float qoffset_z ;    /*!< Quaternion z shift.   */
  
   float srow_x[4] ;    /*!< 1st row affine transform.   */
   float srow_y[4] ;    /*!< 2nd row affine transform.   */
   float srow_z[4] ;    /*!< 3rd row affine transform.   */

   char intent_name[16];/*!< 'name' or meaning of data.  */

   char magic[4] ;      /*!< MUST be "ni1\0" or "n+1\0". */                                       *)
(* ============================================================================================== *)

$NifTI1HeaderSize = 348;
$NifTI2HeaderSize = 540;
$MinNifTI1Offset = 352;
$MinNifTI2Offset = 544;
Protect[$NifTI1HeaderSize, $NifTI2HeaderSize, $MinNifTI1Offset, $MinNifTI2Offset];

(* #ImportNifTIHeader *****************************************************************************)

ImportNifTI1Header[stream_, opts___Rule] := Check[
  SetStreamPosition[stream, 40];
  Block[
    {$ByteOrdering = If[0 < BinaryRead[stream, "Integer16"] < 8, 1, -1] * $ByteOrdering},
    SetStreamPosition[stream, 0];
    With[
      {raw = ReadBinaryStructure[
         stream,
         {"HeaderSize" -> {"Integer32", 1, Function[
            If[#[[1]] != $NifTI1HeaderSize, 
              Message[ImportNifTI::badfmt, "Header does not start with 348"],
              #[[1]]]]},
          "DataType" -> {"Integer8", 10},
          "DBName" -> {"Character8", 18, BinaryStringFix},
          "Extents" -> "Integer32",
          "SessionError" -> "Integer16",
          "Regular" -> "Integer8",
          "DimensionsInformation" -> "Integer8",
          "Dimensions" -> {"Integer16", 8, #[[2 ;; (#[[1]] + 1)]]&},
          "IntentParameters" -> {"Real32", 3},
          "IntentCode" -> "Integer16",
          "Datatype" -> {"Integer16", 1, NifTIDatatypeTranslate[#[[1]]]&},
          "BitsPerVoxel" -> "Integer16",
          "SliceStart" -> "Integer16",
          "VoxelDimensions" -> {"Real32", 8},
          "VoxelOffset" -> {"Real32", 1, Max[{#[[1]], $MinNifTI1Offset}]&},
          "ScaleSlope" -> "Real32",
          "ScaleIntercept" -> "Real32",
          "SliceEnd" -> "Integer16",
          "SliceCode" -> "Integer8",
          "DimensionsUnits" -> {"Integer8", 1, NifTIUnitTranslate[#[[1]]]&},
          "DisplayIntensityRange" -> {"Real32", 2},
          "SliceDuration" -> "Real32",
          "TimeAxisShift" -> "Real32",
          "GLRange" -> {"Integer32", 2, Reverse},
          "Description" -> {"Character8", 80, BinaryStringFix},
          "AuxiliaryFilename" -> {"Character8", 24, BinaryStringFix},
          "QFormCode" -> "Integer16",
          "SFormCode" -> "Integer16",
          "Quaternions" -> {"Real32", 6},
          "AffineTransform" -> {"Real32", 12, Partition[#, 4]&},
          "IntentName" -> {"Character8", 16, BinaryStringFix},
          "MagicTerminus" -> {"Character8", 4, Function[
            With[
              {str = BinaryStringFix[#]},
              If[str == "ni1" || str == "n+1", 
                str,
                Message[ImportNifTI::badfmt, "header does not end in ni1 or n+1 string"]]]]}}]},
      If[raw === $Failed,
        raw,
        "Header" -> Join[
          With[
            {dims = "Dimensions" /. raw,
             pixdim = "VoxelDimensions" /. raw},
            Append[
              Replace[
                raw,
                ("VoxelDimensions" -> d_) :> ("VoxelDimensions"->d[[2 ;; (Length[dims] + 1)]]),
                {1}],
              "Pixdim0" -> pixdim[[1]]]],
          {"ByteOrdering" -> $ByteOrdering,
           "Extension" -> With[
             {ext = NestWhileList[
                Function @ With[
                  {esize = BinaryRead[stream, "Integer32"],
                   ecode = BinaryRead[stream, "Integer32"]},
                  {BinaryReadList[stream, "UnsignedInteger8", esize - 12],
                   BinaryReadList[stream, "Integer8", 4]}],
                {{}, BinaryReadList[stream, "Integer8", 4]},
                #[[2,1]] != 0 &]},
             If[Length[ext] == 1, None, ext[[2;;All, 1]]]]}]]]],
  $Failed];

(* Here we have the NifTI-2 header:
                         /***************************/ /**********************/ /************/
struct nifti_2_header {  /* NIFTI-2 usage           */ /* NIFTI-1 usage      */ /*  offset  */
                         /***************************/ /**********************/ /************/
   int   sizeof_hdr;     /*!< MUST be 540           */ /* int sizeof_hdr; (348) */   /*   0 */
   char  magic[8] ;      /*!< MUST be valid signature. */  /* char magic[4];    */   /*   4 */
   int16_t datatype;     /*!< Defines data type!    */ /* short datatype;       */   /*  12 */
   int16_t bitpix;       /*!< Number bits/voxel.    */ /* short bitpix;         */   /*  14 */
   int64_t dim[8];       /*!< Data array dimensions.*/ /* short dim[8];         */   /*  16 */
   double intent_p1 ;    /*!< 1st intent parameter. */ /* float intent_p1;      */   /*  80 */
   double intent_p2 ;    /*!< 2nd intent parameter. */ /* float intent_p2;      */   /*  88 */
   double intent_p3 ;    /*!< 3rd intent parameter. */ /* float intent_p3;      */   /*  96 */
   double pixdim[8];     /*!< Grid spacings.        */ /* float pixdim[8];      */   /* 104 */
   int64_t vox_offset;   /*!< Offset into .nii file */ /* float vox_offset;     */   /* 168 */
   double scl_slope ;    /*!< Data scaling: slope.  */ /* float scl_slope;      */   /* 176 */
   double scl_inter ;    /*!< Data scaling: offset. */ /* float scl_inter;      */   /* 184 */
   double cal_max;       /*!< Max display intensity */ /* float cal_max;        */   /* 192 */
   double cal_min;       /*!< Min display intensity */ /* float cal_min;        */   /* 200 */
   double slice_duration;/*!< Time for 1 slice.     */ /* float slice_duration; */   /* 208 */
   double toffset;       /*!< Time axis shift.      */ /* float toffset;        */   /* 216 */
   int64_t slice_start;  /*!< First slice index.    */ /* short slice_start;    */   /* 224 */
   int64_t slice_end;    /*!< Last slice index.     */ /* short slice_end;      */   /* 232 */
   char  descrip[80];    /*!< any text you like.    */ /* char descrip[80];     */   /* 240 */
   char  aux_file[24];   /*!< auxiliary filename.   */ /* char aux_file[24];    */   /* 320 */
   int qform_code ;      /*!< NIFTI_XFORM_* code.   */ /* short qform_code;     */   /* 344 */
   int sform_code ;      /*!< NIFTI_XFORM_* code.   */ /* short sform_code;     */   /* 348 */
   double quatern_b ;    /*!< Quaternion b param.   */ /* float quatern_b;      */   /* 352 */
   double quatern_c ;    /*!< Quaternion c param.   */ /* float quatern_c;      */   /* 360 */
   double quatern_d ;    /*!< Quaternion d param.   */ /* float quatern_d;      */   /* 368 */
   double qoffset_x ;    /*!< Quaternion x shift.   */ /* float qoffset_x;      */   /* 376 */
   double qoffset_y ;    /*!< Quaternion y shift.   */ /* float qoffset_y;      */   /* 384 */
   double qoffset_z ;    /*!< Quaternion z shift.   */ /* float qoffset_z;      */   /* 392 */
   double srow_x[4] ;    /*!< 1st row affine transform. */  /* float srow_x[4]; */   /* 400 */
   double srow_y[4] ;    /*!< 2nd row affine transform. */  /* float srow_y[4]; */   /* 432 */
   double srow_z[4] ;    /*!< 3rd row affine transform. */  /* float srow_z[4]; */   /* 464 */
   int slice_code ;      /*!< Slice timing order.   */ /* char slice_code;      */   /* 496 */
   int xyzt_units ;      /*!< Units of pixdim[1..4] */ /* char xyzt_units;      */   /* 500 */
   int intent_code ;     /*!< NIFTI_INTENT_* code.  */ /* short intent_code;    */   /* 504 */
   char intent_name[16]; /*!< 'name' or meaning of data. */ /* char intent_name[16]; */  /* 508 */
   char dim_info;        /*!< MRI slice ordering.   */      /* char dim_info;        */  /* 524 */
   char unused_str[15];  /*!< unused, filled with \0 */                                  /* 525 */
   } ;
*)
ImportNifTI2Header[stream_, opts___Rule] := Check[
  SetStreamPosition[stream, 40];
  Block[
    {$ByteOrdering = If[0 < BinaryRead[stream, "Integer64"] < 8, 1, -1] * $ByteOrdering},
    SetStreamPosition[stream, 0];
    With[
      {raw = ReadBinaryStructure[
         stream,
         {"HeaderSize" -> {"Integer32", 1, Function[
            If[#[[1]] != $NifTI2HeaderSize, 
              Message[ImportNifTI::badfmt, "Header does not start with 540"],
              #[[1]]]]},
          "MagicString" -> {"Character8", 4, Function[
            With[
              {str = BinaryStringFix[#[[1;;4]]]},
              If[str == "ni2" || str == "n+2", 
                str,
                Message[ImportNifTI::badfmt, "header does not end in nil or n+1 string"]]]]},
          "MagicSignature" -> {"Integer8", 4, Function[
            If[# == {13, 10, 26, 10},
              #,
              Message[ImportNifTI::badfmt, "nifti2 magic signature invalid"]]]},
          "Datatype" -> "Integer16",
          "BitsPerVoxel" -> "Integer16",
          "Dimensions" -> {"Integer64", 8, #[[2 ;; (#[[1]] + 1)]]&},
          "IntentParameters" -> {"Real64", 3},
          "VoxelDimensions" -> {"Real64", 8},
          "VoxelOffset" -> {"Integer64", 1, Max[{#[[1]], $MinNifTI2Offset}]&},          
          "ScaleSlope" -> "Real64",
          "ScaleIntercept" -> "Real64",
          "DisplayIntensityRange" -> {"Real64", 2, Reverse},
          "SliceDuration" -> "Real64",
          "TimeAxisShift" -> "Real64",
          "SliceStart" -> "Integer64",
          "SliceEnd" -> "Integer64",
          "Description" -> {"Character8", 80, BinaryStringFix},
          "AuxiliaryFilename" -> {"Character8", 24, BinaryStringFix},
          "QFormCode" -> "Integer32",
          "SFormCode" -> "Integer32",
          "Quaternions" -> {"Real64", 6},
          "AffineTransform" -> {"Real64", 12, Partition[#, 4]&},
          "SliceCode" -> "Integer32",
          "DimensionsUnits" -> {"Integer32", 1, NifTIUnitTranslate[#[[1]]]&},
          "IntentCode" -> "Integer32",
          "IntentName" -> {"Character8", 16, BinaryStringFix},
          "DimensionsInformation" -> "Integer8",
          "UnusedTerminus" -> {"Character8", 15}}]},
      If[raw === $Failed,
        raw,
        "Header" -> Join[
          With[
            {dims = "Dimensions" /. raw,
             pixdims = "VoxelDimensions" /. raw},
            Append[
              Replace[
                raw, 
                ("VoxelDimensions" -> d_) :> ("VoxelDimensions" -> d[[2 ;; (Length[dims] + 1)]]),
                {1}],
              "Pixdim0" -> pixdims[[1]]]],
          {"ByteOrdering" -> $ByteOrdering,
           "Extension" -> With[
             {ext = NestWhileList[
                Function @ With[
                  {esize = BinaryRead[stream, "Integer32"],
                   ecode = BinaryRead[stream, "Integer32"]},
                  {Which[
                    ecode == $NifTIExtensionCodeCifTI, ImportString[
                      StringJoin @ BinaryReadList[stream, "Character8", esize - 12],
                      "XML"],
                    True, BinaryReadList[stream, "UnsignedInteger8", esize - 12]],
                   BinaryReadList[stream, "Integer8", 4]}],
                {{}, BinaryReadList[stream, "Integer8", 4]},
                #[[2,1]] != 0 &]},
             If[Length[ext] == 1, None, ext[[2;;All, 1]]]]}]]]],
  $Failed];

ImportNifTIHeader[stream_, opts___Rule] := Check[
  SetStreamPosition[stream, 0];
  With[
    {sz = BinaryRead[stream, "Integer32"]},
    Which[
      sz == $NifTI1HeaderSize, ImportNifTI1Header[stream, opts],
      sz == $NifTI2HeaderSize, ImportNifTI2Header[stream, opts],
      True, Message[ImportNifTI::badfmt, "Header size is invalid!"]]],
  $Failed];

ImportNifTIMetaInformation[stream_, opts___Rules] := Check[
  With[
    {header = Replace[
       "Header" /. {opts, "Header" :> ("Header" /. ImportNifTIHeader[stream, opts])},
       "Header" :> Message[ImportNifTI::badfmt, "Invalid header"]]},
    "MetaInformation" -> {"Header" -> header}],
  $Failed];

ImportNifTIVoxels[stream_, opts___Rule] := Check[
  With[
    {header = Replace[
       "Header" /. {opts, "Header" :> ("Header" /. ImportNifTIHeader[stream, opts])},
       "Header" :> Message[ImportNifTI::badfmt, "Invalid header"]]},
    Block[
      {$ByteOrdering = ("ByteOrdering" /. header)},
      With[
        {offset = "VoxelOffset" /. header,
         bitpix = "BitsPerVoxel" /. header,
         datatype = "Datatype" /. header,
         dims = "Dimensions" /. header},
        SetStreamPosition[stream, Round[offset]];
        With[
          {raw = BinaryReadList[stream, datatype, Times @@ dims]},
          "Voxels" -> NifTIColorTranslate[
            datatype,
            Transpose[
              Fold[
                Partition,
                raw,
                Most @ Reverse[
                  If[Length[dims] >= 4 && dims[[4]] == 1, Delete[dims, 4], dims]]],
              {3,2,1}]]]]]],
  $Failed];

ImportNifTIData[stream_, opts___Rule] := Check[
  With[
    {header = Replace[
       "Header" /. {opts, "Header" :> ("Header" /. ImportNifTIHeader[stream, opts])},
       "Header" :> Message[ImportNifTI::badfmt, "Invalid header"]]},
    With[
      {voxels = ImportNifTIVoxels[stream, "Header" -> header, opts]},
      "Data" -> {"Header" -> header, voxels}]],
  $Failed];

InterpretNifTI[data_List] := With[
  {header = "Header" /. data,
   voxels = "Voxels" /. data},
  With[
    {dims = Dimensions[voxels]},
    If[Count[dims, 1, {1}] == Length[dims] - 1,
      Flatten[voxels],
      With[
        {qcode = "QFormCode" /. header,
         scode = "SFormCode" /. header},
      MRImage3D[
        Map[Reverse, voxels, {0,1}],
        VoxelDimensions -> ("VoxelDimensions" /. header),
        Sequence @@ If[qcode > 0, 
          With[
            {mtx = QuaternionToRotationMatrix[("Quaternions" /. header)[[1;;3]]]},
            {RightDirectionVector -> mtx[[1]],
             AnteriorDirectionVector -> mtx[[2]],
             SuperiorDirectionVector -> mtx[[3]]}],
          {}],
        MetaInformation -> header]]]]];

ImportNifTI[stream_, opts___Rule] := Check[
  With[
    {data = "Data" /. ImportNifTIData[stream, opts]},
    InterpretNifTI[data]],
  $Failed];

Protect[ImportNifTIHeader, ImportNifTI1Header, ImportNifTI2Header,
        ImportNifTIVoxels, ImportNifTIMetaInformation, 
        ImportNifTI, InterpretNifTI];

ImportExport`RegisterImport[
  "NifTI",
  {"MetaInformation" :> ImportNifTIMetaInformation, 
   "Header" :> ImportNifTIHeader,
   "Voxels" :> ImportNifTIVoxels,
   "Data" :> ImportNifTIData,
   ImportNifTI},
  "FunctionChannels" -> {"Streams"},
  "BinaryFormat" -> True];


(* #GifTI File Format *****************************************************************************)

GifTIExtractMetaData[xml_] := With[
  {meta = FirstCase[xml, XMLElement["MetaData", attr_, data_] :> {data, attr}, None, Infinity]},
  If[meta === None,
    None,
    If[meta[[2]] == {}, #, Append[#, MetaInformation -> meta[[2]]]]& @ Cases[
      meta[[1]],
      RuleDelayed[
        XMLElement["MD", _, {XMLElement["Name", _, {name_}], XMLElement["Value", _, {val_}]}],
        name -> val],
      Infinity]]];
Protect[GifTIExtractMetaData];

GifTIExtractData[xml_] := FirstCase[
  xml,
  XMLElement["Data", _, data_] :> If[ListQ[data], data[[1]], data],
  None,
  Infinity];
Protect[GifTIExtractData];

GifTIExtractDataSpace[xml_] := FirstCase[
  xml,
  XMLElement["DataSpace", _, data_] :> If[ListQ[data], data[[1]], data],
  None,
  Infinity];
Protect[GifTIExtractDataSpace];

GifTIExtractTransformedSpace[xml_] := FirstCase[
  xml,
  XMLElement["TransformedSpace", _, data_] :> If[ListQ[data], data[[1]], data],
  None,
  Infinity];
Protect[GifTIExtractTransformedSpace];

GifTIExtractMatrixData[xml_] := FirstCase[
  xml,
  XMLElement["MatrixData", _, data_] :> If[ListQ[data], data[[1]], data],
  None,
  Infinity];
Protect[GifTIExtractMatrixData];

GifTIExtractCoordinateSystemTransformMatrix[xml_] := Replace[
  Cases[
    xml,
    XMLElement["CoordinateSystemTransformMatrix", _, data_] :> {
      "DataSpace" -> GifTIExtractDataSpace[data],
      "TransformedSpace" -> GifTIExtractTransformedSpace[data],
      "MatrixData" -> GifTIExtractMatrixData[data]},
    Infinity],
  {} -> None];
Protect[GifTIExtractCoordinateSystemTransformMatrix];

GifTIExtractDataArray[xml_] := Map[
  #[[1,1]] -> #[[All, 2]] &,
  GatherBy[
    Cases[
      xml,
      XMLElement["DataArray", attr_, data_] :> With[
        {arrayIndexOrd = Replace[
           "ArrayIndexingOrder" /. attr,
           {"RowMajorOrder" -> (Fold[Partition, #1, Rest @ #2]&),
            "ColumnMajorOrder" -> Transpose[
              Fold[Partition, #1, Reverse @ Rest @ #2],
              Reverse[Range[Length@#2]]],
            _ :> Message[ImportGifTI::badfmt, "Unrecognized ArrayIndexingOrder in DataArray"]}],
         datatype = Replace[
           "DataType" /. attr,
           {"NIFTI_TYPE_UINT8" -> "UnsignedInteger8",
            "NIFTI_TYPE_INT32" -> "Integer32",
            "NIFTI_TYPE_FLOAT32" -> "Real32",
            _ :> Message[ImportGifTI::badfmt, "Unrecognized DataType in DataArray"]}],
         dimensionality = ToExpression["Dimensionality" /. attr],
         encoding = "Encoding" /. attr,
         byteOrder = ("Endian" /. attr) /. {"BigEndian" -> 1, "LittleEndian" -> -1},
         intent = Replace[
           "Intent" /. attr,
           {"NIFTI_INTENT_GENMATRIX" -> {"Tensor", Identity},
            "NIFTI_INTENT_LABEL" -> {"Labels", Identity},
            "NIFTI_INTENT_NODE_INDEX" -> {"Mask", (# + 1)&},
            "NIFTI_INTENT_POINTSET" -> {"Points", Identity},
            "NIFTI_INTENT_RGB_VECTOR" -> {"Overlay", Map[Apply[RGBColor, #]&, #, {-2}]&},
            "NIFTI_INTENT_RGBA_VECTOR" -> {"Overlay", Map[Apply[RGBColor, #]&, #, {-2}]&},
            "NIFTI_INTENT_SHAPE" -> {"Shape", Identity},
            "NIFTI_INTENT_TIME_SERIES" -> {"TimeSeries", Identity},
            "NIFTI_INTENT_TRIANGLE" -> {"Faces", (# + 1)&},
            "NIFTI_INTENT_VECTOR" -> {"Vectors", Identity},
            (s_String /; StringTake[s, 12] == "NIFTI_INTENT") -> {"Overlay", Identity},
            _ -> {"Other", Identity}}]},
        With[
          {dims = Table[ToExpression[("Dim"<>ToString[k]) /. attr], {k, 0, dimensionality - 1}]},
          Block[
            {$ByteOrdering = byteOrder},
            intent[[1]] -> With[
              {subdata = GifTIExtractData[data]},
              With[
                {decoded = arrayIndexOrd[
                   Switch[
                     encoding,
                     "ASCII", Flatten @ ImportString[subdata, "Table"],
                     "Base64Binary", ImportString[data, {"Base64", datatype}],
                     "GZipBase64Binary", ImportString[
                       FromCharacterCode[Base64ZLibDecode[subdata]],
                       datatype]],
                   dims]},
                {"Data" -> (intent[[2]])[decoded],
                 "CoordinateSystemTransformMatrix" -> GifTIExtractCoordinateSystemTransformMatrix[data],
                 "MetaData" -> GifTIExtractMetaData[xml],
                 MetaInformation -> attr}]]]]],
      Infinity],
    First]];
Protect[GifTIExtractDataArray];
  
GifTIExtractLabels[xml_] := SortBy[
  Cases[
    xml,
    XMLElement["Label", attr_, ___] :> attr,
    Infinity],
  ToExpression["Key" /. #]&];
Protect[GifTIExtractLabels];

GifTIExtractLabelTable[xml_] := FirstCase[
  xml,
  XMLElement["LabelTable", _, data_] :> GifTIExtractLabels[data],
  None,
  Infinity];
Protect[GifTIExtractLabelTable];

GifTIParseXML[xml_] := FirstCase[
  xml,
  XMLElement["GIFTI", attrs_, data_] :> With[
    {extr = {
       "LabelTable" -> GifTIExtractLabelTable[data],
       "DataArray" -> GifTIExtractDataArray[data],
       "MetaData" -> GifTIExtractMetaData[data]}},
    If[attrs == {}, extr, Append[extr, MetaInformation -> attrs]]],
  $Failed,
  Infinity];
Protect[GifTIParseXML];

ImportGifTIData[xml_] := Check[
  With[
    {data = GifTIParseXML[xml]},
    data],
  $Failed];
Protect[ImportGifTIData];

GifTILabelTableToAssociation[table_List] := MimicAssociation[
  Map[
    Function[("Key" /. #) -> RGBColor["Red" /. #, "Green" /. #, "Blue" /. #, "Alpha" /. #]],
    table]];
Protect[GifTILabelTableToAssociation];

GifTIConstructSurface[points_, faces_, overlays_, args___] := Module[
  {k = 1},
  Fold[
    Function @ Switch[
      Length[#2], 
      Length[points], SetProperty[{#1, VertexList}, ("Data" <> ToString[k++]) -> #2],
      Length[faces], SetProperty[{#1, FaceList}, ("Data" <> ToString[k++]) -> #2],
      _, #1],
    CorticalMesh[points, faces, args],
    overlays]];
Protect[GifTIConstructSurface];

InterpretGifTIData[xmlData_] := Check[
  With[
    {label = If[# === None, Missing["KeyAbsent",#]&, GifTILabelTableToAssociation[#]]&[
       FirstCase[xmlData, ("LabelTable" -> lbl_) :> lbl, None, Infinity]]},
    With[
      {dataArray = FirstCase[xmlData, ("DataArray" -> q_) :> q, $Failed, Infinity]},
      If[Length[dataArray] == 0, 
        Message[ImportGifTI::badfmt, "No Data element found in GifTI"]];
      With[
        {points = FirstCase[dataArray, ("Points" -> q_) :> q, {}, {1}],
         faces = FirstCase[dataArray, ("Faces" -> q_) :> q, {}, {1}],
         overlays = FirstCase[dataArray, ("Overlay" -> q_) :> q, {}, {1}],
         masks = FirstCase[dataArray, ("Mask" -> q_) :> q, {}, {1}],
         labels = FirstCase[dataArray, ("Labels" -> q_) :> q, {}, {1}]},
        Which[
          (* one surface... *)
          Length[points] == 1 && Length[faces] == 1, GifTIConstructSurface[
            "Data" /. points[[1]], 
            "Data" /. faces[[1]], 
            Join[
              ("Data" /. #)& /@ overlays,
              (Normal @ SparseArray["Data" /. #, {Length[points[[1]]]}])& /@ masks,
              label[[#]]& /@ labels]
            MetaInformation -> {"GifTIData" -> q}],
          (* exactly one other data type; we want to yield just the data *)
          Lengh[Join[points, faces, overlays, masks, labels]] == 1, First @ Which[
            Length[points] > 0, points,
            Length[faces] > 0, faces,
            Length[overlays] > 0, overlays,
            Length[masks] > 0, masks,
            True, labels],
          (* Otherwuse, we don't have a clear interpretation, so just yield the data *)
          True, xmlData]]]],
  $Failed];
Protect[InterpretGifTIData];

ImportGifTI[filename_, opts___Rule] := Check[
  InterpretGifTIData[ImportGifTIData[Import[filename, "XML", opts]]],
  $Failed];
Protect[ImportGifTI];

(* Register the GifTI importer... *)
ImportExport`RegisterImport[
  "GifTI",
  {"Data" :> ImportGifTIData,
   ImportGifTI}];


End[];
EndPackage[];

