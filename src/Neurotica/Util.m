(* Util.m
 *
 * Basic utility functions, e.g. for automatically caching data in quick-to-load mx files.
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
BeginPackage["Neurotica`Util`", {"Neurotica`Global`"}];
Unprotect["Neurotica`Util`*", "Neurotica`Util`Private`*"];
ClearAll ["Neurotica`Util`*", "Neurotica`Util`Private`*"];

$CacheDirectory::usage = "$CacheDirectory is the directory in which cached data for this user is placed. If the directory does not exist at cache time and $AutoCreateCacheDirectory is True, then this directory is automatically created. By default it is the directory FileNameJoin[{$UserBaseDirectory, \"AutoCache\"}]. If set to Temporary, then the next use of AutoCache will create a temporary directory and set the $CacheDirectory to this path.";

$AutoCreateCacheDirectory::usage = "$AutoCreateCacheDirectory is True if and only if the $CacheDirectory should be automatically created when AutoCache requires it.";

AutoCache::usage = "AutoCache[name, body] checks first to see if a cache file exist for the file named by the given name and yields its contents if so. Otherwise, yields the result of evaluating body and caches it in the given filename. If the cache file is out of date relative to the evaluating cell, the contents are erased and recalculated. If the given filename is an absolute path, it is used as such, otherwise it is localized to the cache directory. Note that the mx extension is added automatically if not included.
The following options may be used:
 * Quiet (default: False) prevents the function from producing messages when updating or overwriting the cache file if False.
 * Check (default: True) if True instructs AutoCache to yield $Failed and not generate the cache file whenever messages are produced during the execution of body.
 * Directory (default: Automatic) names the directory in which the cache file should go; if Automatic, then uses $CacheDirectory.
 * CreateDirectory (default: Automatic) determines whether AutoCache can automatically create the $CacheDirectory if it is needed. If True, then AutoCache will always create the directory; if False it will never create the directory; if Automatically, it will defer to $AutoCreateCacheDirectory.";
AutoCache::expired = "AutoCache file `1` has expired; recalculating";
AutoCache::nodir = "AutoCache directory (`1`) does not exist and cannot be created";
AutoCache::nomkdir = "AutoCache directory (`1`) does not exist and AutoCache is not allowed to create it";
AutoCache::badopt = "Bad option to AutoCache: `1`"

AutoCacheFilename::usage = "AutoCacheFilename[name] yields the filename that would be used to store the auto-cache with the given name.
AutoCacheFilename accepts the same arguments, Directory and CreateDirectory, as AutoCache and will create directories is specified.";

AutoCachedQ::usage = "AutoCachedQ[name] yields True if and only if the given cache exists in the AutoCache directories.
The Directory option may be passed as with AutoCache.";

SetSafe::usage = "SetSafe[a, b] is equivalent t0 the expression (a = b) except that if any messages are generated during the evaluation of b, $Failed is yielded and a is not set.";

Memoize::usage = "Memoize[f := expr] constructs a version of the pattern f that will auto-memoize its values and remember them.
Memoize[\!\(\*SubscriptBox[\(f\),\(1\)]\) := \!\(\*SubscriptBox[\(expr\),\(1\)]\), \!\(\*SubscriptBox[\(f\),\(2\)]\) := \!\(\*SubscriptBox[\(expr\),\(2\)]\), ...] memoizes each of the given functions.
The function Forget may be used to forget the memoized values of functions defined within a Memoize form.";
MemoizeSafe::usage = "MemoizeSafe[f := expr, ...] is identical to calling Memoize except that it uses SetSafe instead of Set to memoize values, so of the evaluation of the expression fails or generates a message, the result is not memoized and $Failed is yielded.";
Forget::usage = "Forget[f] forces all memoized values associated with the symbol f to be forgotten.
Forget[\!\(\*SubscriptBox[\(f\),\(1\)]\), \!\(\*SubscriptBox[\(f\),\(2\)]\), ...] forgets each of the symbols.";

Indices::usage = "Indices[list, patt] yields a list of the indices of elements that match patt in list; this is equivalent to Flatten[Position[list, patt, {1}, Heads->False]], but is considerably optimized.";
Index::usage = "Index[list] yields an association containing, as keys, each unique element that occurs in the given list and, as values, the indices at which the corresponding element occurs in the list. For example, Index[{a,b,c,b}] == <|a -> {1}, b -> {2,4}, c -> {3}|>.";

TemporarySymbol::usage = "TemporarySymbol[] yields a unique symbol that is marked as temporary.
TemporarySymbol[string] yields a unique symbol that begins with the given string and is marked temporary.";

DefineImmutable::usage = "DefineImmutable[constructor[args...] :> symbol, {members...}] defines an immutable data type with the given constructor and the accessor functions given by the list of members. Each member must be of the form (f[symbol] ** form) where f is the function that is used to access the data, the form yields the value and may refer to the arguments of the constructor or the other member functions, and ** must be either =, ->, :>, or :=, each of which is interpreted differently:
 * any value assigned with a := will become a generic function whose value is calculated each time the function is called;
 * any value assigned with a :> rule becomes a function that is memoized for each object once it has been calculated; note, however, that memoized values in an immutable object cannot be forgotted by the Forget[] function;
 * any value assigned with a -> rule will be calculated and saved to the object when it is constructed;
 * any value assigned with a = will also be assigned when the object is constructed, just as with the -> rule forms; however, values assigned with an = form may be modified in the Clone[] function.
In addition, the functions assigned using the := form may have additional arguments, as long as the first argument is the symbol.
DefineImmutable[constructor[args...] :> symbol, Hold[members...]] is identical to DefineImmutable[constructor[args...] :> symbol, {members...}].
The following options may be given:
 * Symbol (default: Automatic) specifies the box form used to hold the type; if no box type is given, then a unique symbol is generated to encapsulate the type.
 * SetSafe (default: True) specifies whether memoization performed on :> rule forms should be done using SetSafe (True) or Set (False).
See also Clone.";
DefineImmutable::badarg = "Bad argument given to DefineImmutable: `1`";
Clone::usage = "Clone[immutable, edits...] yields an immutable object (created via a constructor defined by the DefineImmutable interface) identical to the given object immutable with the changes given by the sequence of edits. Each edit must take the form accessor -> value where accessor is the function name for a member assigned using an = form in the DefineImmutable function and value is the new value it should take.";
Clone::badarg = "Bad argument given to Clone: `1`";

Let::usage = "Let[{bindings...}, body] is similar to With, Module, and Block, but includes the following features:
  * Each binding may refer to any binding that precedes it in the list of bindings.
  * Bindings of the form symbol = form may be reassigned in the body (like a Module binding).
  * Bindings of the form symbol -> form may not be reassigned in the body (like a With binding).
  * Bindings of the form symbol := form may be reassigned but are initially assigned using SetDelayed.
  * Bindings of the form symbol :> form may not be reassigned, but are not assigned until requested, then are memoized.";

NormalizeRows::usage = "NormalizeRows[X] yields a transformation of the matrix X in which each row of X has been normalized; this is the equivalent of (Normalize /@ X) but is significantly optimized.";
NormalizeColumns::usage = "NormalizeColumns[X] yields a transformation of the matrix X in which each column of X has been normalized. This is equivalent to Transpose[Normalze /@ Transpose[X]], but has been significantly optimized.";
RowNorms::usage = "RowNorms[X] yields the equivalent of Norm /@ X.";
ColumnNorms::usage = "ColumnNorms[X] yields the equivalent of Norm /@ Transpose[X] but has been optimized for speed.";

QuaternionToRotationMatrix::usage = "QuaternionToRotationMatrix[{x,y,z,w}] yields the rotation matrix associated with the quaternion {x,y,z,w}.
QuaternionToRotationMatrix[{y, z, q}] is equivalent to QuaternionToRotationMatrix[{x,y,z,q}] where x is Sqrt[1 - (x^2 + y^2 + z^2)].";

ReadBinaryStructure::usage = "ReadBinaryStructure[stream, instructionsList] yields the data structure result of importing the given instructions via BinaryReadList. Instructions may take the following format:
  * type (e.g., \"Integer32\" or \"Real64\") are read as is;
  * {type, n} reads in a list of n of the given type (if n is 1 then a 1-element list is returned);
  * {type, n, postprocessFn} yields postprocessFn[data] where data is the list of n of the given type;
  * name -> spec yields the rule name -> result where result is the result of reading in the given instruction spec;
  * {spec1, spec2...} yields a list of the given instruction specs, each read in order.";
ReadBinaryStructure::badinstr = "Bad instruction given to ReadBinaryStructure: `1`";

BinaryStringFix::usage = "BinaryStringFix[{c1, c2, ...}] yields the string formed from the given list of characters, truncating the string upon encountering the first null character (character code 0).";

NeuroticaPermanentDatum::usage = "NeuroticaPermanentDatum[name, value] yields Null but adds the given name/value pair to the current user's Neurotica auto-initialization file. The data contained in this file are loaded at startup and must consist only of string data. The file itself is stored as a \"Table\" style file at FileNameJoin[{$UserBasePath, \".NeuroticaPermanentData.dat\"}].
If the value given in the second argument is None, then the datum is removed from the permanent cache.
NeuroticaPermanentDatum[name] yields the value of the permanent datum with the given name; if no such datum is found, then None is returned.";

FlipIntegerBytes::usage = "FlipIntegerBytes[i] yields a version of the integer i such that the bytes are interpreted in reverse order; i.e., switches between Big and Little Endian.
FlipIntegerBytes[i, k] treats i as a k-byte integer; k may be a number (e.g., 4 for a 32-bit integer) or a data type (e.g., \"Integer32\" for a 32-bit integer).
FlipIntegerBytes[i] is equivalent to FlipIntegerBytes[i, \"Integer32\"].";
FlipIntegerBytes::badarg = "Bad argument given to FlipIntegerBytes: `1`";

Unboole::usage = "Unboole[u] yields a value or list the same length as the list u in which all zero elements have been converted to False and all non-zero elements have been converted to True; effectively this is the opposite of the Boole function. Unboole automatically threads over lists.
Note: If the argument to Unboole is not a real number or an integer, then Unboole yields a Piecewise function that checks for the argument == 0.";

MapNamed::usage = "MapNamed[f, {name1 -> list1, name2 -> list2, ...}] is identical to MapThread[f, {list1, list2...}] except that instead of passing the members of the various lists as the first, second, etc. argument to f, MapNamed passes a single Association argument to f with the appropriate values named according to the names in the given list.
MapNamed[f, name -> list] is equivalent to MapNamed[f, {name -> list}].
MapNamed[f, {data...}, level] applies f to the given level of the lists in data (equivalent to MapThread's use of level).
MapNamed[f, name -> list, level] also applies f to the given level, but expects level to be equivalent to the level passed to Map.";

DivideCheck::usage = "DivideCheck[a,b] yields a/b if Chop[b] is not equal to 0 and yields 0 if Chop[b] is equal to 0.
DivideCheck[a,b,c] yields c if Chop[b] is equal to 0.
Note that DivideCheck works with arrays as well as single values.";

FlatOuter::usage = "FlatOuter[args...] is identical to Outer[args...] except that it yields a 1D map that is equivalent to a flattened version of the latter form.";
FlatTable::usage = "FlatTable[args...] is identical to Table[args...] except that it yields a 1D map that is equivalent to a flattened version of the latter form.";

ForwardOptions::usage = "ForwardOptions[f[args...], opts...] yields the result of calling the function f with the given arguments in args... followed by the options in the sequence of opts... where the options are filtered: ForwardOptions[f[args...], opts...] is equivalent to f[args..., Sequence@@FilterRules[{opts}, Options[f]].";

$NeuroticaPermanentData::usage = "$NeuroticaPermanentData is an Association of the permanent Neurotica data, as saved using the NeuroticaPermanentDatum function.";

Begin["`Private`"];

(* #FlipIntegerBytes ******************************************************************************)
(* This code was adapted (by Noah C. Benson, <nben@nyu.edu>) from a pull request given by 
   Szabolcs Horvát <szhorvat@gmail.com> in order to match the Neurotica library aesthetic *)
FlipIntegerBytes[i_Integer, byteCount_Integer /; 0 < byteCount] := If[2^(8*byteCount) <= i,
  Message[FlipIntegerBytes::badarg, "given integer is out of range for number of bytes"],
  FromDigits[Reverse@IntegerDigits[i, 256, byteCount], 256]];
FlipIntegerBytes[i_Integer, type_String] := Replace[
  type,
  {"Integer8"   :> FlipIntegerBytes[i, 1],
   "Integer16"  :> FlipIntegerBytes[i, 2],
   "Integer24"  :> FlipIntegerBytes[i, 3],
   "Integer32"  :> FlipIntegerBytes[i, 4],
   "Integer64"  :> FlipIntegerBytes[i, 8],
   "Integer128" :> FlipIntegerBytes[i, 16],
   _ :> Message[
     FlipIntegerBytes::badarg, 
     "type must be an integer or \"IntegerN\" where N is 8, 16, 32, ... 128."]}];
FlipIntegerBytes[i_Integer] := FlipIntegerBytes[i, 4];
Protect[FlipIntegerBytes];

(* #$CacheDirectory *******************************************************************************)
$CacheDirectory = FileNameJoin[{$UserBaseDirectory, "AutoCache"}];

(* #$AutoCreateCacheDirectory *********************************************************************)
$AutoCreateCacheDirectory = True;
$AutoCreateCacheDirectory /: Set[$AutoCreateCacheDirectory, b:(True|False)] := (
  Unprotect[$AutoCreateCacheDirectory];
  $AutoCreateCacheDirectory = b;
  Protect[$AutoCreateCacheDirectory];
  b);
Protect[$AutoCreateCacheDirectory];

(* #AutoCacheFilename *****************************************************************************)
Options[AutoCacheFilename] = {Directory -> Automatic, CreateDirectory -> Automatic};
AutoCacheFilename[name_String, OptionsPattern[]] := With[
  {splitName = FileNameSplit[name],
   canCreate = Replace[
       OptionValue[CreateDirectory],
       {Automatic :> $AutoCreateCacheDirectory,
        Except[True|False] :> (
          Message[AutoCache::badopt, "CreateDirectory must be True, False, or Automatic"];
          Throw[$Failed])}]},
  With[
    {dir = Replace[
       OptionValue[Directory],
       {Automatic :> Which[
          First[splitName] == "", FileNameJoin[Most[splitName]],
          First[splitName] == "~", FileNameJoin[Most[splitName]],
          $CacheDirectory === Temporary, ($CacheDirectory = CreateDirectory[]),
          True, $CacheDirectory],
        Except[_String] :> (
          Message[AutoCache::badopt, "Directory must be Automatic or a string"];
          Throw[$Failed])}],
     file = StringReplace[Last[splitName], {s:(__ ~~ ".mx") :> s, s__ :> (s <> ".mx")}, 1]},
    Catch@FileNameJoin[
      {Which[
         DirectoryQ[dir], dir,
         canCreate, Replace[
           CreateDirectory[dir],
           $Failed :> (Message[AutoCache::nodir, dir]; Throw[$Failed])],
         True, (Message[AutoCache::nomkdir, dir]; Throw[$Failed])],
       file}]]];
Protect[AutoCacheFilename];

(* #AutoCache *************************************************************************************)
Options[AutoCache] = {
  Quiet -> False,
  Check -> True,
  Directory -> Automatic,
  CreateDirectory -> Automatic};
Attributes[AutoCache] = {HoldRest};
AutoCache[name_String, body_, OptionsPattern[]] := Catch@With[
  {quiet = TrueQ[OptionValue[Quiet]],
   check = TrueQ[OptionValue[Check]],
   fileName = AutoCacheFilename[
     name,
     Directory -> OptionValue[Directory],
     CreateDirectory -> OptionValue[CreateDirectory]]},
  With[
    {cachedRes = If[StringQ[fileName] && FileExistsQ[fileName],
       Block[{Global`data = None}, Get[fileName]; Global`data],
       None]},
    With[
      {theRes = If[cachedRes === None || cachedRes === $Failed,
         If[check, Check[body, $Failed], body],
         cachedRes]},
      If[cachedRes === None || cachedRes == $Failed,
        Block[{Global`data = theRes}, If[theRes =!= $Failed, DumpSave[fileName, Global`data]]]];
      theRes]]];
Protect[AutoCache];

(* #AutoCachedQ ***********************************************************************************)
Options[AutoCachedQ] = {Directory -> Automatic};
AutoCachedQ[name_String, opts:OptionsPattern[]] := Quiet@Check[
  FileExistsQ@AutoCacheFilename[name, CreateDirectory -> False, opts],
  False];
Protect[AutoCachedQ];

(* #SetSafe ***************************************************************************************)
SetAttributes[SetSafe, HoldAll];
SetSafe[f_, expr_] := With[
  {res = Check[expr, $Failed]},
  If[res === $Failed, $Failed, Set[f, res]]];
Protect[SetSafe];

(* #Memoize ***************************************************************************************)
$MemoizePatternName = $MemoizePatternName;
Protect[$MemoizePatternName];
SetAttributes[Memoize, HoldAll];
Memoize[forms:(SetDelayed[_,_]..)] := CompoundExpression@@Replace[
  Hold[forms],
  HoldPattern[SetDelayed[f_, expr_]] :> ReplacePart[
    Hold[
      Evaluate[Pattern[Evaluate[$MemoizePatternName], f]],
      Evaluate[Hold[Evaluate[$MemoizePatternName], expr]]],
    {{2,0} -> Set, 0 -> SetDelayed}],
  {1}];
Protect[Memoize];

(* #MemoizeSafe ***********************************************************************************)
SetAttributes[MemoizeSafe, HoldAll];
MemoizeSafe[forms:(SetDelayed[_,_]..)] := CompoundExpression@@Replace[
  Hold[forms],
  HoldPattern[SetDelayed[f_, expr_]] :> ReplacePart[
    Hold[
      Evaluate[Pattern[Evaluate[$MemoizePatternName], f]],
      Evaluate[Hold[Evaluate[$MemoizePatternName], expr]]],
    {{2,0} -> SetSafe, 0 -> SetDelayed}],
  {1}];
Protect[MemoizeSafe];

(* #Forget ****************************************************************************************)
SetAttributes[Forget, HoldAll];
Forget[symbol_] := With[
  {sym1 = Unique[],
   sym2 = Unique[],
   down = DownValues[symbol]},
  With[
    {patterns = Cases[
       down,
       ReplacePart[
         Hold[
           Evaluate[
             Pattern @@ {
               sym1,
               Hold[
                 Evaluate[HoldPattern@@{Pattern[Evaluate[$MemoizePatternName], _]}],
                 Evaluate[Hold[Evaluate[$MemoizePatternName], _]]]}],
           Evaluate[Hold[Evaluate[{sym1}], Evaluate[{1,1} -> sym2]]]],
         {{1,2,2,0} -> (Set|SetSafe),
          {1,2,0} -> RuleDelayed,
          {2,1,0} -> First,
          {2,0} -> ReplacePart,
          0 -> RuleDelayed}],
       {1}]},
    With[
      {memoized = Cases[
         down[[All, 1]], 
         Evaluate[(Alternatives @@ patterns) :> Evaluate[Hold[Evaluate[sym2]]]],
         {2}]},
      Scan[Apply[Unset, #]&, memoized]]]];
Protect[Forget];      

(* #Indices ***************************************************************************************)
Indices[list_List, patt_] := Pick[Range[Length@list], list, patt];
Protect[Indices];

(* #TemporarySymbol *******************************************************************************)
TemporarySymbol[] := With[{sym = Unique[]}, SetAttributes[Evaluate[sym], Temporary]; sym];
TemporarySymbol[s_String] := With[{sym = Unique[s]}, SetAttributes[Evaluate[sym], Temporary]; sym];
Protect[TemporarySymbol];

(* #DefineImmutable *******************************************************************************)
SetAttributes[DefineImmutable, HoldAll];
Options[DefineImmutable] = {
  SetSafe -> True,
  Symbol -> Automatic};
DefineImmutable[RuleDelayed[pattern_, sym_Symbol], args_, OptionsPattern[]] := Check[
  Block[
    {sym},
    With[
      {instructions0 = If[#[[1,0]] === Hold, ReplacePart[#, {1,0} -> List], #]&[Hold[args]],
       type = Unique[SymbolName[Head[pattern]]],
       box = Replace[
         OptionValue[Symbol],
         {Automatic :> Unique[ToString[sym]],
          Except[_Symbol] :> Message[DefineImmutable::badarg, "Symbol option must be a symbol"]}],
       setFn = Replace[OptionValue[SetSafe], {True -> SetSafe, Except[True] -> Set}]},
      With[
        {instructions = Thread @ NestWhile[#[[{1},2]]&, instructions0, #[[1,0]] === With&],
         outerWithPos = Last @ NestWhile[
           {(#[[1]])[[{1}, 2]], Append[#[[2]], 2]}&,
           {instructions0, {}},
           (#[[1]])[[1, 0]] === With&]},
        With[
          {patterns = Apply[HoldPattern, #]& /@ instructions[[All, {1}, 1]],
           heads = Replace[
             instructions[[All, {1}, 0]],
             {Hold[Set] -> Set,
              Hold[Rule] -> Rule,
              Hold[RuleDelayed] -> RuleDelayed,
              Hold[SetDelayed] -> SetDelayed,
              x_ :> Message[
                DefineImmutable::badarg, 
                "Instructions bust be ->, =, :=, or :> forms: " <> ToString[x]]},
             {1}],
           bodies = instructions[[All, {1}, 2]]},
          If[Count[patterns, form_ /; form[[1,1]] === sym, {1}] < Length[patterns],
            Message[
              DefineImmutable::badarg, 
              "All field patterns must use the object as their first parameter"]];
          MapThread[
            Function[
              If[#1 =!= SetDelayed && #1 !MatchQ[Hold @@ #2, Hold[_[sym]]],
                Message[DefineImmutable::badarg, "rule and set functions must be unary"]]],
            {heads, patterns}];
          (* Setup the graph of dependencies *)
          With[
            {deps = Flatten@Last@Reap[
               (* We sort instructions into 4 groups:
                *   Init: instructions that must be calculated at startup (a -> b or a = b)
                *   Delay: instructions that needn't be calculated/saved (a := b)
                *   Settable: instructions that can be changed via Clone (a = b)
                *   Lazy: instructions that are delayed then memoized (a :> b)
                * We also collect Edges:
                *   An edge occurs when the RHS of an instruction references another instruction
                *)
               Sow[None, "Edge"];
               Sow[None, "Init"];
               Sow[None, "Delay"];
               Sow[None, "Settable"];
               Sow[None, "Lazy"];
               MapThread[
                 Function[{head, body, id},
                   If[head === SetDelayed,
                     Sow[id, "Delay"],
                     With[
                       {matches = Select[
                          Range[Length@patterns],
                          (heads[[#]] =!= SetDelayed && Count[body,patterns[[#]],Infinity] > 0)&]},
                       If[head =!= RuleDelayed, Sow[id, "Init"], Sow[id, "Lazy"]];
                       If[head === Set, Sow[id, "Settable"]];
                       Scan[
                         Function[
                           If[id == #,
                             Message[DefineImmutable::badarg, "Self-loop detected in dependencies"],
                             Sow[DirectedEdge[#, id], "Edge"]]],
                         matches]]]],
                 {heads, bodies, Range[Length@heads]}],
               {"Edge", "Init", "Delay", "Settable", "Lazy"},
               Rule[#1, Rest[#2]]&]},
            With[
              {delay = "Delay" /. deps,
               members = Complement[Range[Length@heads], "Delay" /. deps],
               init = "Init" /. deps,
               lazy = "Lazy" /. deps,
               lazyIdx = Part[
                 SortBy[
                   Join[
                     MapIndexed[#1 -> {2, #2[[1]]}&, "Lazy" /. deps],
                     MapIndexed[#1 -> {1, #2[[1]]}&, "Init" /. deps]],
                   First],
                 All, 2]},
              (* The easiest part to setup is the delayed items; do them now *)
              Do[
                TagSetDelayed @@ Join[
                  Hold[box],
                  ReplacePart[Hold @@ {patterns[[id]]}, {1,1,1} -> Pattern@@{sym, Blank[box]}],
                  bodies[[id]]],
                {id, "Delay" /. deps}];
              (* The direct accessors are also pretty easy to create... 
                 We create these slightly differently for lazy and non-lazy members. *)
              Do[
                TagSetDelayed @@ Join[
                  Hold[box],
                  ReplacePart[
                    Hold @@ {patterns[[members[[id]]]]}, 
                    {1,1,1} -> Pattern@@{sym, Blank[box]}],
                  ReplacePart[
                    (* For lazy: sym[[2, k]]; for init: sym[[1, k]] *)
                    Hold @@ {Hold[sym, Evaluate[lazyIdx[[id, 1]]], Evaluate[lazyIdx[[id, 2]]]]},
                    {1,0} -> Part]],
                {id, Range[Length@members]}];
              (* The constructor requires that we evaluate things in certain orders... 
                 we start by making some symols to use for each value *)
              With[
                {syms = Table[Unique[], {Length@members}],
                 (* We also make a graph of the dependencies *)
                 depGraph = With[
                   {tmp = Graph[members, "Edge" /. deps]},
                   If[!AcyclicGraphQ[tmp],
                     Message[DefineImmutable::badarg, "Cyclical dependency detected"],
                     tmp]]},
                (* with the transitive closure of the dependency graph, we know all the downstream
                   dependencies of any given element *)
                With[
                  {depStream = With[
                     {tc = TransitiveClosureGraph[depGraph]},
                     Map[
                       Function[
                         {EdgeList[tc, DirectedEdge[_,#]][[All, 1]],
                          EdgeList[tc, DirectedEdge[#,_]][[All, 2]]}],
                       members]],
                   (* Given the dependency graph, we can also determine the order non-delayed values
                      need to be evaluated in to prevent dependency errors *)
                   initOrder = Reverse@Part[
                     NestWhileList[
                       Function@With[
                         {goal = #[[1]], known = #[[3]]},
                         With[
                           {knowable = Select[
                              goal,
                              Function[
                                0 == Length@Complement[
                                  EdgeList[depGraph, DirectedEdge[_,#]][[All,1]], 
                                  known]]]},
                           {Complement[goal, knowable], knowable, Union[Join[known, knowable]]}]],
                       With[
                         {easy = Select[
                            "Init" /. deps,
                            EdgeCount[depGraph, DirectedEdge[_,#]] == 0&]},
                         {Complement["Init" /. deps, easy], easy, easy}],
                       Length[#[[1]]] > 0&],
                     All, 2]},
                  (* okay, now we build the constructor;
                     this is a big ugly code chunk due to the macro-nature of it... *)
                  SetDelayed @@ Join[
                    Hold[pattern],
                    With[
                      {constructorBody = ReplaceAll[
                         (* To construct the body of the custructor:
                          * (1) make a core: a body that sets the lazy symbols and makes the form
                          * (2) wrap the core in layers of With[{init...}, core]
                          *)
                         Fold[
                           (* (2): wrap the core-so-far in the next outer layer of init vars *)
                           Function[
                             ReplacePart[
                               Hold@Evaluate@Join[
                                 Hold@Evaluate@Map[
                                   Join[Hold@@{syms[[#]]}, bodies[[members[[#]]]]]&,
                                   #2],
                                 #1],
                               {{1,1,_,0} -> Set,
                                {1,0} -> With}]],
                           (* (1): make a core that sets the lazy symbols and yields the form *)
                           With[
                             {idcs = Map[
                                Position[members,#][[1,1]]&,
                                Complement[members, Flatten[initOrder]]],
                              formSym = Unique[],
                              lazyRevIdx = Part[
                                Map[
                                  Rest[SortBy[#[[All, 2 ;; 3]], Last]]&,
                                  SortBy[
                                    GatherBy[
                                      Join[
                                        MapIndexed[{#1[[1]], #2[[1]], #1[[2]]}&, lazyIdx],
                                        {{1, 0, 0}, {2, 0, 0}}],
                                      First],
                                    #[[1,1]]&]],
                                All, All, 1]},
                             ReplacePart[
                               Hold@Evaluate@Join[
                                 ReplacePart[
                                   Hold@Evaluate[Hold[#1,Unique[]]& /@ syms[[idcs]]],
                                   {1,_,0} -> Set],
                                 (* Here we make the actual form, prior to setting the symbols so
                                    that they are easier to hold: *)
                                 With[
                                   {form = Hold@Evaluate[
                                      box @@ ReplacePart[
                                        syms[[#]]& /@ lazyRevIdx,
                                        {2,0} -> Hold]]},
                                   ReplacePart[
                                     Hold@Evaluate@Join[
                                       ReplacePart[
                                         Hold @@ {{Flatten[Hold[Hold@Evaluate@formSym, form]]}},
                                         {1,1,0} -> Set],
                                       ReplacePart[
                                         Hold@Evaluate@Join[
                                           Flatten[
                                             Hold @@ MapThread[
                                               Function[
                                                 ReplacePart[
                                                   ReplaceAll[
                                                     Hold@Evaluate@Hold[
                                                       {#1}, 
                                                       Evaluate[Join[Hold[#1], #2]]],
                                                     {sym -> formSym}],
                                                   {{1,2,0} -> setFn,
                                                    {1,0} -> SetDelayed,
                                                    {1,1,0} -> Evaluate}]],
                                               {syms[[idcs]], bodies[[members[[idcs]]]]}]],
                                           Hold@Evaluate@formSym],
                                         {1,0} -> CompoundExpression]],
                                     {1,0} -> With]]],
                               {1,0} -> With]],
                           Map[Position[members, #][[1,1]]&, initOrder, {2}]],
                         MapThread[Rule, {patterns[[members]], syms}]]},
                      (* At this point, we need to return a Hold[] pattern, and constructorBody is
                         all that we really need; however, we may need to wrap some extra data
                         around our definitions... *)
                      If[outerWithPos == {}, 
                        constructorBody,
                        ReplacePart[
                          ReplacePart[instructions0, Prepend[outerWithPos, 1] -> constructorBody],
                          Join[{1}, outerWithPos, {0}] -> Identity]]]];
                  (* we also want to build a clone function *)
                  With[
                    {settables = "Settable" /. deps,
                     patternNames = #[[1,0]]& /@ patterns},
                    Evaluate[box] /: Clone[b_box, changes___Rule] := Check[
                      With[
                        {rules = {changes}},
                        Which[
                          (* No changes to be made; just return the current object *)
                          rules == {}, b,
                          (* multiple changes: do them one at a time *)
                          Length[rules] > 1, Fold[Clone, b, rules],
                          (* otherwise, we have one change; just handle it *)
                          True, With[
                            {id = Replace[
                               Select[
                                 Indices[patternNames, rules[[1,1]]],
                                 MemberQ[settables, #]&],
                               {{i_} :> i,
                                _ :> Message[
                                  Clone::badarg,
                                  "unrecognized or unsettable rule: " <> ToString[rules[[1,1]]]]}],
                             val = rules[[1,2]]},
                            (* Iterate through the values that depend on this one *)
                            With[
                              {harvest = Reap[
                                 Fold[
                                   Function@With[
                                     {bCur = #1, iids = #2},
                                     ReplacePart[
                                       bCur,
                                       Map[
                                         Function[
                                           lazyIdx[[ Position[members, #1][[1,1]] ]] -> If[
                                             heads[[#1]] === RuleDelayed,
                                             With[
                                               {tmp = TemporarySymbol[]},
                                               Sow[tmp -> bodies[[#1]]];
                                               tmp],
                                             ReleaseHold@ReplaceAll[
                                               bodies[[#1]],
                                               sym -> bCur]]],
                                         iids]]],
                                   ReplacePart[b, lazyIdx[[ Position[members, id][[1,1]] ]] -> val],
                                   Rest@First@Last@Reap[
                                     NestWhile[
                                       Function[{G},
                                         With[
                                           {vs = Select[VertexList[G], VertexInDegree[G,#] == 0&]},
                                           Sow[vs];
                                           VertexDelete[G, vs]]],
                                       Subgraph[
                                         depGraph,
                                         VertexOutComponent[depGraph, id, Infinity]],
                                       VertexCount[#] > 0&]]]]},
                              (* harvest[[1]] is the new box'ed object; we need to setup the 
                                 temporary symbols, though. *)
                              With[
                                {finalSym = Unique[]},
                                Scan[
                                  Function@With[
                                    {memSym = #[[1]], memBody = #[[2]]},
                                    ReplacePart[
                                      {memSym, Join[
                                         Hold@Evaluate@memSym,
                                         ReplaceAll[memBody, sym -> finalSym]]},
                                      {{2,0} -> Set,
                                       0 -> SetDelayed}]],
                                  Flatten[harvest[[2]]]];
                                Set @@ Join[Hold[Evaluate@finalSym], Hold[harvest][[{1},1]]];
                                finalSym]]]]],
                      $Failed]];
                SetAttributes[Evaluate[box], Protected];
                True]]]]]]]],
  $Failed];

(* #Let *******************************************************************************************)
SetAttributes[Let, HoldAll];
With[
  {tagCompoundExpr = Unique["ComplexExpression"],
   tagWith = Unique["With"]},
  Let[
    {variables:((_Symbol|((Set|SetDelayed|Rule|RuleDelayed)[_,_]))...)},
    body_
   ] := ReplaceRepeated[
    Fold[
     Function@Replace[
       #2,
       {Hold[s_Symbol] :> Hold[
          tagWith,
          {Set[s, Unique@ToString[s]]},
          SetAttributes[s, Temporary],
          #1],
        Hold[Rule[s_Symbol, expr_]] :> Hold[tagWith, {Set[s, expr]}, #1],
        Hold[Set[s_Symbol, expr_]] :> Hold[
          tagWith,
          {Set[s, Unique@ToString[s]]},
          SetAttributes[s, Temporary],
          s = expr,
          #1],
        Hold[RuleDelayed[s_Symbol, expr_]] :> Hold[
          tagWith,
          {Set[s, Unique@ToString[s]]},
          SetAttributes[s, Temporary],
          SetDelayed[s, Set[s, expr]],
          #1],
        Hold[SetDelayed[s_Symbol, expr_]] :> Hold[
          tagWith,
          {Set[s, Unique@ToString[s]]},
          SetAttributes[s, Temporary],
          SetDelayed[s, expr],
          #1],
        Hold[Rule[s_Symbol[args___], expr_]] :> Hold[
          tagWith,
          {Set[s, Unique@ToString@s]},
          Set[s[args], expr],
          SetAttributes[s, {Temporary, Protected}],
          #1],
        Hold[Set[s_Symbol[args___], expr_]] :> Hold[
          tagWith,
          {Set[s, Unique@ToString@Head[s]]},
          SetAttributes[s, Temporary],
          s[args] = expr,
          #1],
        Hold[RuleDelayed[s_Symbol[args___], expr_]] :> With[
          {pattSym = Unique["patt"]},
          Hold[
           tagWith,
           {Set[s, Unique@ToString[s]]},
           SetAttributes[s, Temporary],
           SetDelayed[pattSym : s[args], Set[pattSym, expr]],
           #1]],
        Hold[SetDelayed[s_Symbol[args___], expr_]] :> Hold[
          tagWith,
          {Set[s, Unique@ToString[s]]},
          SetAttributes[s, Temporary],
          SetDelayed[s[args], expr],
          #1]}],
     Hold[tagCompoundExpr, body],
     Hold /@ Reverse@Hold[variables]],
    {Hold[tagCompoundExpr, b__] :> CompoundExpression[b],
     Hold[tagWith, v_, b_] :> With[v, b],
     Hold[tagWith, v_, b_, bs__] :> With[v, CompoundExpression[b, bs]]}]];
SyntaxInformation[Let] = {"ArgumentsPattern" -> {{__},_}};
Protect[Let];

(* #NormalizeRows *********************************************************************************)
NormalizeRows[X_] := With[
  {tr = Transpose[X]},
  Transpose[tr / Table[#, {Length[tr]}]&@Sqrt@Total[tr^2]]];
NormalizeRows[{}] = {};
Protect[NormalizeRows];

(* #NormalizeColumns ******************************************************************************)
NormalizeColumns[X_] := X / ConstantArray[# + (1 - Unitize[#]), {Length[X]}]&@Sqrt@Total[X^2];
NormalizeColumns[{}] = {};
Protect[NormalizeColumns];

(* #RowNorms **************************************************************************************)
RowNorms[X_] := Norm /@ X;
RowNorms[{}] = {};
Protext[RowNorms];

(* #ColumnNorms ***********************************************************************************)
ColumnNorms[Xt_] := Sqrt @ Total[Xt^2];
ColumnNorms[{}] = {};
Protext[ColumnNorms];

(* #QuaternionToRotationMatrix *******************************************************************)
QuaternionToRotationMatrix[{a_, b_, c_, d_}] := {
  {a^2 + b^2 - c^2 - d^2, 2*b*c - 2*a*d, 2*b*d + 2*a*c},
  {2*b*c + 2*a*d, a^2 + c^2 - b^2 - d^2,   2*c*d - 2*a*b},
  {2*b*d - 2*a*c, 2*c*d + 2*a*b, a^2 + d^2 - c^2 - b^2}};
QuaternionToRotationMatrix[{b_, c_, d_}] := QuaternionToRotationMatrix[
  {Sqrt[1.0 - (b^2 + c^2 + d^2)], b, c, d}];
Protect[QuaternionToRotationMatrix];

(* #BinaryStringFix *******************************************************************************)
BinaryStringFix[clist:{_String..}] := StringJoin[TakeWhile[clist, ToCharacterCode[#] != {0}&]];
BinaryStringFix[clist:{_Integer..}] := FromCharacterCode @ TakeWhile[clist, # != 0&];
Protect[BinaryStringFix];

(* #ReadBinaryStructure ***************************************************************************)
$ReadBinaryStructureInstructionKey = {
  type_String :> Function[BinaryRead[#, type]],
  {type_String} :> Function[BinaryReadList[#, type, 1]],
  {type_String, n_Integer} :> Function[BinaryReadList[#, type, n]],
  {type_String, n_Integer, f_} :> Function[f[BinaryReadList[#, type, n]]],
  l_List :> With[
    {els = Replace[l, $ReadBinaryStructureInstructionKey, 1]},
    Function[With[{stream = #}, #[stream]& /@ els]]],
  Rule[name_, instr_] :> With[
    {r = (instr /. $ReadBinaryStructureInstructionKey)}, 
    Function[name -> r[#]]],
  q_ :> Message[ReadBinaryStructure::badinstr, ToString[q]]};
Protect[$ReadBinaryStructureInstructionKey];
  
ReadBinaryStructure[stream_, instructions_List] := Check[
  With[
    {instFns = Replace[instructions, $ReadBinaryStructureInstructionKey, {1}]},
    Map[
      #[stream]&,
      instFns]],
  $Failed];
Protect[ReadBinaryStructure];

(* #NeuroticaPermanentDatum ***********************************************************************)
$NeuroticaPermanentDataPath = FileNameJoin[{$UserBaseDirectory, ".NeuroticaPermanentData.dat"}];
$NeuroticaPermanentData = If[FileExistsQ[$NeuroticaPermanentDataPath],
  Association @ Map[
    Function[#[[1]] -> StringJoin[Riffle[Rest[#], "\t"]]],
    Import[$NeuroticaPermanentDataPath, "Table"]],
  Association[]];
Protect[$NeuroticaPermamentDataPath, $NeuroticaPermamentData];

NeuroticaPermanentDataSave[] := Check[
  (If[Length[$NeuroticaPermanentData] == 0,
     If[FileExistsQ[$NeuroticaPermanentDataPath], DeleteFile[$NeuroticaPermanentDataPath]],
     Export[
       $NeuroticaPermanentDataPath,
       Map[(List @@ #)&, Normal[$NeuroticaPermanentData]],
       "Table"]];
   True),
  False];
Protect[NeuroticaPermanentDataSave];

NeuroticaPermanentDatum[name_String, value_String] := With[
  {cur = Lookup[$NeuroticaPermanentData, name]},
  Unprotect[$NeuroticaPermanentData];
  $NeuroticaPermanentData[name] = value;
  Protect[$NeuroticaPermanentData];
  NeuroticaPermanentDataSave[]];
NeuroticaPermanentDatum[name_String, None] := With[
  {cur = Lookup[$NeuroticaPermanentData, name]},
  Unprotect[$NeuroticaPermanentData];
  KeyDropFrom[$NeuroticaPermanentData, name];
  Protect[$NeuroticaPermanentData];
  NeuroticaPermanentDataSave[]];
NeuroticaPermanentDatum[name_String] := Replace[
  $NeuroticaPermanentData[name],
  Missing[___] -> None];
Protect[NeuroticaPermanentDatum];


(* #Unboole ***************************************************************************************)
With[
  {FT = {False,True}},
  Unboole[x_Real] := FT[[Unitize[x] + 1]];
  Unboole[x_Integer] := FT[[Unitize[x] + 1]];
  Unboole[x_?NumericQ] := FT[[Unitize[x] + 1]];
  Unboole[x_] := Piecewise[{{False, x == 0}}, True];
  Unboole[x_List /; ArrayQ[x, 1, NumericQ]] := FT[[Unitize[x] + 1]];
  Unboole[x_List] := Map[Unboole, x]];
Protect[Unboole];

(* #MapNamed **************************************************************************************)
MapNamed[f_, args:{Rule[_,_List]..}, level_:1] := With[
  {names = args[[All, 1]]},
  MapThread[
    Function[f@Association@Thread[names -> {##}]],
    args[[All,2]],
    level]];
MapNamed[f_, Rule[name_, list_List], level_:{1}] := Map[
  Function[f[<|name -> #|>]],
  list,
  level];
Protect[MapNamed];

(* #DivideCheck ***********************************************************************************)
DivideCheck[a_, b_, c_:0] := With[
  {bne0 = Unitize@Chop[b]},
  With[
    {beq0 = 1 - unit},
    c * beq0 + bneq0 * a / (b + beq0)]];
SetAttributes[DivideCheck, NumericFunction];
Protect[DivideCheck];

(* #Index *****************************************************************************************)
Index[data_List] := Association@Last@Reap[
  MapThread[Sow, {Range@Length[data], data}],
  _,
  Rule];
Protect[Index];

(* #FlatOuter *************************************************************************************)
FlatOuter[args__] := With[
  {listArgs = TakeWhile[Rest[{args}], ListQ]},
  Flatten[Outer[args], Length[listArgs] - 2]];
Protect[FlatOuter];

(* #FlatTable *************************************************************************************)
FlatTable[args__] := Flatten[
  Table[args],
  Length@Hold[args] - 2];
SetAttributes[FlatTable, HoldAll];
Protect[FlatTable];

(* #ForwardOptions ********************************************************************************)
Attributes[ForwardOptions] = {HoldFirst};
ForwardOptions[f_[args___], opts___] := With[
  {head = f},
  head[args, Sequence@@FilterRules[{opts}, Options[head]]]];
Protect[ForwardOptions];

End[];
EndPackage[];
