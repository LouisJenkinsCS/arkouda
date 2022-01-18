/*
 * Copyright 2020-2021 Hewlett Packard Enterprise Development LP
 * Copyright 2004-2019 Cray Inc.
 * Other additional copyright holders may be indicated within.
 *
 * The entirety of this work is licensed under the Apache License,
 * Version 2.0 (the "License"); you may not use this file except
 * in compliance with the License.
 *
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

// Note: Need `--ccflags="_GNU_SOURCE"` to get the definition of `fallocate`

extern proc fopen(name: c_string, mode: c_string): c_void_ptr;
extern proc fread(data:c_void_ptr, size: int, n: int, f: c_void_ptr): int;
extern proc fwrite(data:c_void_ptr, size: int, n: int, f: c_void_ptr): int;
extern proc fclose(f: c_void_ptr) : int;
extern proc sizeof(x): int;
extern proc fflush(f: c_void_ptr) : int;
extern proc fsync(fd: int) : int;
extern proc fileno(f: c_void_ptr) : int;
extern proc fseek(f: c_void_ptr, offset : int, origin : int) : int;
extern proc ftell(f: c_void_ptr) : int;
extern proc ftruncate(fd: int, length: int): int;
extern proc fallocate(fd: int, mode: int, offset: int, length: int): int;
extern proc mmap(addr : c_void_ptr, length : int, prot : int, flags : int, fd : int, offset : int) : c_void_ptr;
extern proc munmap(addr : c_void_ptr, length : int) : int;
extern proc unlink(fname : c_string) : int;
extern const SEEK_SET : int;
extern const SEEK_CUR : int;
extern const SEEK_END : int;
extern const MAP_SHARED : int;
extern const PROT_READ : int;
extern const PROT_WRITE : int;
extern const MAP_FAILED : int;
extern const FALLOC_FL_ZERO_RANGE : int;
extern const FALLOC_FL_PUNCH_HOLE : int;
extern const FALLOC_FL_KEEP_SIZE : int;
extern const FALLOC_FL_COLLAPSE_RANGE : int;
extern const errno : int;


private use DSIUtil;
private use ChapelLocks;
private use CPtr;
private use SysCTypes;
private use Sys;

proc _determineRankFromStartIdx(startIdx) param {
  return if isTuple(startIdx) then startIdx.size else 1;
}

proc _determineIdxTypeFromStartIdx(startIdx) type {
  return if isTuple(startIdx) then startIdx(0).type else startIdx.type;
}

config param debugCyclicMMAPDist = false;
config param verboseCyclicMMAPDistWriters = false;
config param debugCyclicMMAPDistBulkTransfer = false;

//
// If the testFastFollowerOptimization flag is set to true, the
// follower will write output to indicate whether the fast follower is
// used or not.  This is used in regression testing to ensure that the
// 'fast follower' optimization is working.
//
config param testFastFollowerOptimization = false;

//
// This flag is used to disable lazy initialization of the RAD cache.
//
config param disableCyclicMMAPLazyRAD = defaultDisableLazyRADOpt;

// chpldoc TODO:
//   a good reference to
//     dataParTasksPerLocale, dataParIgnoreRunningTasks, dataParMinGranularity
//   supports RAD opt, Bulk Transfer optimization, localSubdomain
//   disableCyclicMMAPLazyRAD
//
/*
This CyclicMMAP distribution maps indices to locales in a round-robin pattern
starting at a given index.

Formally, consider a CyclicMMAP distribution with:

  =============  ====================================================
  rank           ``d``
  start index    ``(s_1, ...., s_d)``
  over locales   ``targetLocales: [0..N_1-1, ...., 0..N_d-1] locale``
  =============  ====================================================

It maps an index ``(i_1, ...., i_d)``
to the locale ``targetLocales[j_1, ...., j_d]``
where, for each ``k`` in ``1..d``,
we have:

  ``j_k = (i_k - s_k) (mod N_k)``


**Example**

The following code declares a domain ``D`` distributed
using a CyclicMMAP distribution with a start index of ``(1,1)``,
and declares an array ``A`` over that domain.
The `forall` loop sets each array element
to the ID of the locale to which it is mapped.

  .. code-block:: chapel

    use CyclicMMAPDist;

    const Space = {1..8, 1..8};
    const D: domain(2) dmapped CyclicMMAP(startIdx=Space.low) = Space;
    var A: [D] int;

    forall a in A do
      a = a.locale.id;

    writeln(A);

When run on 6 locales, the output is:

  ::

    0 1 0 1 0 1 0 1
    2 3 2 3 2 3 2 3
    4 5 4 5 4 5 4 5
    0 1 0 1 0 1 0 1
    2 3 2 3 2 3 2 3
    4 5 4 5 4 5 4 5
    0 1 0 1 0 1 0 1
    2 3 2 3 2 3 2 3


**Initializer Arguments**

The ``CyclicMMAP`` class initializer is defined as follows:

  .. code-block:: chapel

    proc CyclicMMAP.init(
      startIdx,
      targetLocales: [] locale = Locales,
      dataParTasksPerLocale     = // value of  dataParTasksPerLocale      config const,
      dataParIgnoreRunningTasks = // value of  dataParIgnoreRunningTasks  config const,
      dataParMinGranularity     = // value of  dataParMinGranularity      config const,
      param rank: int  = // inferred from startIdx argument,
      type idxType     = // inferred from startIdx argument )

The argument ``startIdx`` is a tuple of integers defining an index that
will be distributed to the first locale in ``targetLocales``.
In the 1-dimensional case, ``startIdx`` can be an integer
or a tuple with a single element.

The argument ``targetLocales`` is an array containing the target
locales to which this distribution maps indices and data.
The rank of ``targetLocales`` must match the rank of the distribution,
or be ``1``.  If the rank of ``targetLocales`` is ``1``, a greedy
heuristic is used to reshape the array of target locales so that it
matches the rank of the distribution and each dimension contains an
approximately equal number of indices.

The arguments ``dataParTasksPerLocale``, ``dataParIgnoreRunningTasks``,
and ``dataParMinGranularity`` set the knobs that are used to
control intra-locale data parallelism for CyclicMMAP-distributed domains
and arrays in the same way that the like-named config constants
control data parallelism for ranges and default-distributed domains
and arrays.

The ``rank`` and ``idxType`` arguments are inferred from the
``startIdx`` argument unless explicitly set.
They must match the rank and index type of the domains
"dmapped" using that CyclicMMAP instance.


**Convenience Initializer Functions**

It is common for a ``CyclicMMAP`` distribution to distribute its indices
across all locales. In this case, a convenience function can be used to
declare variables of CyclicMMAP-distributed domain or array type.  These functions
take a domain or list of ranges as arguments and return a CyclicMMAP-distributed
domain or array.

  .. code-block:: chapel

    use CyclicMMAPDist;

    var CyclicMMAPDom1 = newCyclicMMAPDom({1..5, 1..5});
    var CyclicMMAPArr1 = newCyclicMMAPArr({1..5, 1..5}, real);
    var CyclicMMAPDom2 = newCyclicMMAPDom(1..5, 1..5);
    var CyclicMMAPArr2 = newCyclicMMAPArr(1..5, 1..5, real);


**Data-Parallel Iteration**

A `forall` loop over a CyclicMMAP-distributed domain or array
executes each iteration on the locale where that iteration's index
is mapped to.

Parallelism within each locale is guided by the values of
``dataParTasksPerLocale``, ``dataParIgnoreRunningTasks``, and
``dataParMinGranularity`` of the respective CyclicMMAP instance.
Updates to these values, if any, take effect only on the locale
where the updates are made.


**Limitations**

This distribution has not been tuned for performance.
*/

var uniqueFileSuffix : atomic int;
config const mmapFilePrefix = "";
class CyclicMMAP: BaseDist {
  param rank: int;
  type idxType = int;

  const targetLocDom: domain(rank);
  var targetLocs: [targetLocDom] locale;
  var startIdx: rank*idxType;

  var locDist: [targetLocDom] unmanaged LocCyclicMMAP(rank, idxType);

  var dataParTasksPerLocale: int;
  var dataParIgnoreRunningTasks: bool;
  var dataParMinGranularity: int;

  proc init(startIdx,
            targetLocales: [] locale = Locales,
            dataParTasksPerLocale=getDataParTasksPerLocale(),
            dataParIgnoreRunningTasks=getDataParIgnoreRunningTasks(),
            dataParMinGranularity=getDataParMinGranularity(),
            param rank: int = _determineRankFromStartIdx(startIdx),
            type idxType = _determineIdxTypeFromStartIdx(startIdx))
    where isTuple(startIdx) || isIntegral(startIdx)
  {
    this.rank = rank;
    this.idxType = idxType;

    const ranges = setupTargetLocRanges(rank, targetLocales);
    this.targetLocDom = {(...ranges)};
    this.targetLocs = reshape(targetLocales, this.targetLocDom);

    var startIdxTemp: rank*idxType;
    for param i in 0..rank-1 {
      const startIdxI = if isTuple(startIdx) then startIdx(i) else startIdx;
      startIdxTemp(i) = chpl__mod(startIdxI, targetLocDom.dim(i).sizeAs(int));
    }
    this.startIdx = startIdxTemp;

    // Instead of 'dummyLC', we could give 'locDistTemp' a nilable element type.
    const dummyLC = new unmanaged LocCyclicMMAP(rank, idxType, dummy=true);
    var locDistTemp: [targetLocDom] unmanaged LocCyclicMMAP(rank, idxType)
          = dummyLC;
    coforall locid in targetLocDom do
      on targetLocs(locid) do
       locDistTemp(locid) =
         new unmanaged LocCyclicMMAP(rank, idxType, locid, startIdxTemp, ranges);

    delete dummyLC;
    this.locDist = locDistTemp;

    // NOTE: When these knobs stop using the global defaults, we will need
    // to add checks to make sure dataParTasksPerLocale<0 and
    // dataParMinGranularity<0
    this.dataParTasksPerLocale = if dataParTasksPerLocale==0
                                 then here.maxTaskPar
                                 else dataParTasksPerLocale;
    this.dataParIgnoreRunningTasks = dataParIgnoreRunningTasks;
    this.dataParMinGranularity = dataParMinGranularity;

    this.complete();

    if debugCyclicMMAPDist then
      for loc in locDist do writeln(loc);
  }

  proc dsiAssign(other: this.type) {
    coforall locid in targetLocDom do
      on targetLocs(locid) do
        delete locDist(locid);
    startIdx = other.startIdx;
    targetLocDom = other.targetLocDom;
    targetLocs = other.targetLocs;
    dataParTasksPerLocale = other.dataParTasksPerLocale;
    dataParIgnoreRunningTasks = other.dataParIgnoreRunningTasks;
    dataParMinGranularity = other.dataParMinGranularity;
    coforall locid in targetLocDom do
      on targetLocs(locid) do
        locDist(locid) = new unmanaged LocCyclicMMAP(rank, idxType, locid, this);
  }

  proc dsiEqualDMaps(that: CyclicMMAP(?)) {
    return (this.startIdx == that.startIdx &&
            this.targetLocs.equals(that.targetLocs));
  }

  proc dsiEqualDMaps(that) param {
    return false;
  }

  proc dsiClone() {
    return new unmanaged CyclicMMAP(startIdx, targetLocs,
                      dataParTasksPerLocale,
                      dataParIgnoreRunningTasks,
                      dataParMinGranularity);
  }

  override proc dsiDestroyDist() {
    coforall ld in locDist do {
      on ld do
        delete ld;
    }
  }

}

proc CyclicMMAP.chpl__locToLocIdx(loc: locale) {
  for locIdx in targetLocDom do
    if (targetLocs[locIdx] == loc) then
      return (true, locIdx);
  return (false, targetLocDom.first);
}

proc CyclicMMAP.getChunk(inds, locid) {
  const chunk = locDist(locid).myChunk((...inds.getIndices()));
  return chunk;
}

override proc CyclicMMAP.dsiDisplayRepresentation() {
  writeln("startIdx = ", startIdx);
  writeln("targetLocDom = ", targetLocDom);
  writeln("targetLocs = ", for tl in targetLocs do tl.id);
  writeln("dataParTasksPerLocale = ", dataParTasksPerLocale);
  writeln("dataParIgnoreRunningTasks = ", dataParIgnoreRunningTasks);
  writeln("dataParMinGranularity = ", dataParMinGranularity);
  for tli in targetLocDom do
    writeln("locDist[", tli, "].myChunk = ", locDist[tli].myChunk);
}

override proc CyclicMMAPDom.dsiSupportsAutoLocalAccess() param { return true; }

proc CyclicMMAP.init(other: CyclicMMAP, privateData,
                 param rank = other.rank,
                 type idxType = other.idxType) {
  this.rank = rank;
  this.idxType = idxType;
  targetLocDom = {(...privateData[1])};
  targetLocs = other.targetLocs;
  startIdx = privateData[0];
  locDist = other.locDist;
  dataParTasksPerLocale = privateData[2];
  dataParIgnoreRunningTasks = privateData[3];
  dataParMinGranularity = privateData[4];
}
                 
override proc CyclicMMAP.dsiSupportsPrivatization() param return true;

proc CyclicMMAP.dsiGetPrivatizeData() return (startIdx,
                                          targetLocDom.dims(),
                                          dataParTasksPerLocale,
                                          dataParIgnoreRunningTasks,
                                          dataParMinGranularity);

proc CyclicMMAP.dsiPrivatize(privatizeData) {
  return new unmanaged CyclicMMAP(_to_unmanaged(this), privatizeData);
}

proc CyclicMMAP.dsiGetReprivatizeData() return 0;

proc CyclicMMAP.dsiReprivatize(other, reprivatizeData) {
  targetLocDom = other.targetLocDom;
  targetLocs = other.targetLocs;
  locDist = other.locDist;
  startIdx = other.startIdx;
  dataParTasksPerLocale = other.dataParTasksPerLocale;
  dataParIgnoreRunningTasks = other.dataParIgnoreRunningTasks;
  dataParMinGranularity = other.dataParMinGranularity;
}

override proc CyclicMMAP.dsiNewRectangularDom(param rank: int, type idxType, param stridable: bool, inds) {
  if idxType != this.idxType then
    compilerError("CyclicMMAP domain index type does not match distribution's");
  if rank != this.rank then
    compilerError("CyclicMMAP domain rank does not match distribution's");
  const whole = createWholeDomainForInds(rank, idxType, stridable, inds);

  const dummyLCD = new unmanaged LocCyclicMMAPDom(rank, idxType);
  var locDomsTemp: [this.targetLocDom] unmanaged LocCyclicMMAPDom(rank, idxType)
        = dummyLCD;
  coforall localeIdx in this.targetLocDom do
    on this.targetLocs(localeIdx) do
      locDomsTemp(localeIdx) = new unmanaged LocCyclicMMAPDom(rank, idxType,
                                              this.getChunk(whole, localeIdx));
  delete dummyLCD;

  var dom = new unmanaged CyclicMMAPDom(rank, idxType, stridable,
                                    this: unmanaged, locDomsTemp, whole);
  return dom;
}

//
// Given a tuple of scalars of type t or range(t) match the shape but
// using types rangeType and scalarType e.g. the call:
// _matchArgsShape(range(int(32)), int(32), (1:int(64), 1:int(64)..5, 1:int(64)..5))
// returns the type: (int(32), range(int(32)), range(int(32)))
//
proc _CyclicMMAP_matchArgsShape(type rangeType, type scalarType, args) type {
  proc helper(param i: int) type {
    if i == args.size-1 {
      if isCollapsedDimension(args(i)) then
        return (scalarType,);
      else
        return (rangeType,);
    } else {
      if isCollapsedDimension(args(i)) then
        return (scalarType, (... helper(i+1)));
      else
        return (rangeType, (... helper(i+1)));
    }
  }
  return helper(0);
}

proc CyclicMMAP.writeThis(x) throws {
  x <~> this.type:string <~> "\n";
  x <~> "------\n";
  for locid in targetLocDom do
    x <~> " [" <~> locid <~> "=" <~> targetLocs(locid) <~> "] owns chunk: " <~>
      locDist(locid).myChunk <~> "\n";
}

proc CyclicMMAP.targetLocsIdx(i: idxType) {
  const numLocs:idxType = targetLocDom.sizeAs(idxType);
  // this is wrong if i is less than startIdx
  //return ((i - startIdx(0)) % numLocs):int;
  // this works even if i is less than startIdx
  return chpl__diffMod(i, startIdx(0), numLocs):idxType;
}

proc CyclicMMAP.targetLocsIdx(ind: rank*idxType) {
  var x: rank*int;
  for param i in 0..rank-1 {
    var dimLen = targetLocDom.dim(i).sizeAs(int);
    //x(i) = ((ind(i) - startIdx(i)) % dimLen):int;
    x(i) = chpl__diffMod(ind(i), startIdx(i), dimLen):int;
  }
  if rank == 1 then
    return x(0);
  else
    return x;
}

proc CyclicMMAP.dsiIndexToLocale(i: idxType) where rank == 1 {
  return targetLocs(targetLocsIdx(i));
}

proc CyclicMMAP.dsiIndexToLocale(i: rank*idxType) {
  return targetLocs(targetLocsIdx(i));
}


  proc chpl__computeCyclicMMAPDim(type idxType, lo, myloc, numlocs) {
    const lower = min(idxType)..(lo+myloc) by -numlocs;
    const upper = lo+myloc..max(idxType) by numlocs;
    return lower.last..upper.last by numlocs;
  }

proc chpl__computeCyclicMMAP(type idxType, locid, targetLocBox, startIdx) {
    type strType = chpl__signedType(idxType);
    param rank = targetLocBox.size;
    var inds: rank*range(idxType, stridable=true);
    for param i in 0..rank-1 {
      // NOTE: Not bothering to check to see if these can fit into idxType
      const lo = chpl__tuplify(startIdx)(i): idxType;
      const myloc = chpl__tuplify(locid)(i): idxType;
      // NOTE: Not checking for overflow here when casting to strType
      const numlocs = targetLocBox(i).sizeAs(strType);
      inds(i) = chpl__computeCyclicMMAPDim(idxType, lo, myloc, numlocs);
    }
    return inds;
  }

class LocCyclicMMAP {
  param rank: int;
  type idxType;

  const myChunk: domain(rank, idxType, true);

  proc init(param rank, type idxType, locid,
            distStartIdx: rank*idxType, distLocDims) {
    this.rank = rank;
    this.idxType = idxType;

    var locidx: rank*idxType;
    var startIdx = distStartIdx;

    // NOTE: Not bothering to check to see if these can fit into idxType
    if rank == 1 then
      locidx(0) = locid:idxType;
    else
      for param i in 0..rank-1 do locidx(i) = locid(i):idxType;

    var inds: rank*range(idxType, stridable=true);

    inds = chpl__computeCyclicMMAP(idxType, locid, distLocDims, startIdx);
    myChunk = {(...inds)};
  }

  // Used to create a dummy instance.
  proc init(param rank, type idxType, param dummy: bool) where dummy {
    this.rank = rank;
    this.idxType = idxType;
  }
}


class CyclicMMAPDom : BaseRectangularDom {
  const dist: unmanaged CyclicMMAP(rank, idxType);

  var locDoms: [dist.targetLocDom] unmanaged LocCyclicMMAPDom(rank, idxType);

  var whole: domain(rank, idxType, stridable);
}

proc CyclicMMAPDom.setup() {
    coforall localeIdx in dist.targetLocDom {
      on dist.targetLocs(localeIdx) {
        var chunk = dist.getChunk(whole, localeIdx);
        locDoms(localeIdx).myBlock = chunk;
      }
    }
}

override proc CyclicMMAPDom.dsiDestroyDom() {
    coforall localeIdx in dist.targetLocDom {
      on dist.targetLocs(localeIdx) do
        delete locDoms(localeIdx);
    }
}

proc CyclicMMAPDom.dsiBuildArray(type eltType, param initElts:bool) {
  const dom = this;
  const creationLocale = here.id;
  const dummyLCD = new unmanaged LocCyclicMMAPDom(rank, idxType);
  const dummyLCA = new unmanaged LocCyclicMMAPArr(eltType, rank, idxType,
                                              dummyLCD, false);
  var locArrTemp: [dom.dist.targetLocDom]
                    unmanaged LocCyclicMMAPArr(eltType, rank, idxType) = dummyLCA;
  var myLocArrTemp: unmanaged LocCyclicMMAPArr(eltType, rank, idxType)?;

  // formerly in CyclicMMAPArr.setup()
  coforall localeIdx in dom.dist.targetLocDom with (ref myLocArrTemp) {
    on dom.dist.targetLocs(localeIdx) {
      const LCA = new unmanaged LocCyclicMMAPArr(eltType, rank, idxType,
                                             dom.locDoms(localeIdx),
                                             initElts=initElts);
      locArrTemp(localeIdx) = LCA;
      if here.id == creationLocale then
        myLocArrTemp = LCA;
    }
  }
  delete dummyLCA, dummyLCD;

  var arr = new unmanaged CyclicMMAPArr(eltType=eltType, rank=rank,
                                    idxType=idxType, stridable=stridable,
         dom=_to_unmanaged(dom), locArr=locArrTemp, myLocArr=myLocArrTemp);


  return arr;
}

override proc CyclicMMAPDom.dsiDisplayRepresentation() {
  writeln("whole = ", whole);
  for tli in dist.targetLocDom do
    writeln("locDoms[", tli, "].myBlock = ", locDoms[tli].myBlock);
  dist.dsiDisplayRepresentation();
}

proc CyclicMMAPDom.dsiLow return whole.low;

proc CyclicMMAPDom.dsiHigh return whole.high;

proc CyclicMMAPDom.dsiAlignedLow return whole.alignedLow;

proc CyclicMMAPDom.dsiAlignedHigh return whole.alignedHigh;

proc CyclicMMAPDom.dsiAlignment return whole.alignment;

proc CyclicMMAPDom.dsiStride return whole.stride;

proc CyclicMMAPDom.dsiMember(i) return whole.contains(i);

proc CyclicMMAPDom.dsiIndexOrder(i) return whole.indexOrder(i);

proc CyclicMMAPDom.dsiDims() return whole.dims();

proc CyclicMMAPDom.dsiDim(d: int) return whole.dim(d);

proc CyclicMMAPDom.getLocDom(localeIdx) return locDoms(localeIdx);

override proc CyclicMMAPDom.dsiMyDist() return dist;



proc CyclicMMAPDom.dsiGetIndices() {
  return whole.getIndices();
}

proc CyclicMMAPDom.dsiSetIndices(x: domain) {
  whole = x;
  setup();
}

proc CyclicMMAPDom.dsiSetIndices(x) {
  whole.setIndices(x);
  setup();
}

proc CyclicMMAPDom.dsiAssignDomain(rhs: domain, lhsPrivate:bool) {
  chpl_assignDomainWithGetSetIndices(this, rhs);
}

proc CyclicMMAPDom.dsiSerialWrite(x) {
  if verboseCyclicMMAPDistWriters {
    x <~> this.type:string <~> "\n";
    x <~> "------\n";
    for loc in dist.targetLocDom {
      x <~> "[" <~> loc <~> "=" <~> dist.targetLocs(loc) <~> "] owns " <~>
        locDoms(loc).myBlock <~> "\n";
    }
  } else {
    x <~> whole;
  }
}

proc CyclicMMAPDom.doiToString() {
  return whole:string;
}

proc CyclicMMAPDom.dsiNumIndices return whole.sizeAs(uint);

iter CyclicMMAPDom.these() {
  for i in whole do
    yield i;
}

iter CyclicMMAPDom.these(param tag: iterKind) where tag == iterKind.leader {
  const maxTasks = dist.dataParTasksPerLocale;
  const ignoreRunning = dist.dataParIgnoreRunningTasks;
  const minSize = dist.dataParMinGranularity;
  const wholeLow = whole.low;
  const wholeStride = whole.stride;

  // If this is the only task running on this locale, we don't want to
  // count it when we try to determine how many tasks to use.  Here we
  // check if we are the only one running, and if so, use
  // ignoreRunning=true for this locale only.  Obviously there's a bit
  // of a race condition if some other task starts after we check, but
  // in that case there is no correct answer anyways.
  //
  // Note that this code assumes that any locale will only be in the
  // targetLocales array once.  If this is not the case, then the
  // tasks on this locale will *all* ignoreRunning, which may have
  // performance implications.
  const hereId = here.id;
  const hereIgnoreRunning = if here.runningTasks() == 1 then true
                            else ignoreRunning;
  coforall locDom in locDoms do on locDom {
    const myIgnoreRunning = if here.id == hereId then hereIgnoreRunning
      else ignoreRunning;

    // Forward to defaultRectangular to iterate over the indices we own locally
    for followThis in locDom.myBlock.these(iterKind.leader, maxTasks,
                                           myIgnoreRunning, minSize) do {
      // translate the 0-based indices yielded back to our indexing scheme
      const newFollowThis = chpl__followThisToOrig(idxType, followThis, locDom.myBlock);

      // translate the local indices back to 0-based global indices
      // note that we need to go back and forth in order to distinguish
      // between global strides and those that are due to the CyclicMMAP
      // distribution (at least, I couldn't figure out a way to not go
      // back and forth without breaking tests)
      const zeroShift = {(...newFollowThis)}.chpl__unTranslate(wholeLow);
      var result: rank*range(idxType=idxType, stridable=true);
      type strType = chpl__signedType(idxType);
      for param i in 0..rank-1 {
        const wholestride = chpl__tuplify(wholeStride)(i);
        const ref dim = zeroShift.dim(i);
        result(i) = (dim.first / wholestride:idxType)..(dim.last / wholestride:idxType) by (dim.stride:strType / wholestride);
      }
      yield result;
    }
  }
}

// Utility routine to convert 0-based indices back to the indexing scheme
// of 'whole'
private proc chpl__followThisToOrig(type idxType, followThis, whole) {
  param rank = followThis.size;
  var t: rank*range(idxType, stridable=true);
  if debugCyclicMMAPDist then
    writeln(here.id, ": follower whole is: ", whole,
                     " follower is: ", followThis);
  for param i in 0..rank-1 {
    // NOTE: unsigned idxType with negative stride will not work
    const wholestride = whole.dim(i).stride:chpl__signedType(idxType);
    t(i) = ((followThis(i).low*wholestride:idxType)..(followThis(i).high*wholestride:idxType) by (followThis(i).stride*wholestride)) + whole.dim(i).alignedLow;
  }
  return t;
}

iter CyclicMMAPDom.these(param tag: iterKind, followThis) where tag == iterKind.follower {
  const t = chpl__followThisToOrig(idxType, followThis, whole);
  if debugCyclicMMAPDist then
    writeln(here.id, ": follower maps to: ", t);
  for i in {(...t)} do
    yield i;
}

proc CyclicMMAPDom.chpl__serialize() {
  return pid;
}

// TODO: What happens when we try to deserialize on a locale that doesn't
// own a copy of the privatized class?  (can that happen?)  Could this
// be a way to lazily privatize by also making the originating locale part
// of the 'data'?
proc type CyclicMMAPDom.chpl__deserialize(data) {
  return chpl_getPrivatizedCopy(unmanaged CyclicMMAPDom(rank=this.rank,
                                                    idxType=this.idxType,
                                                    stridable=this.stridable),
                                data);
}

override proc CyclicMMAPDom.dsiSupportsPrivatization() param return true;

proc CyclicMMAPDom.dsiGetPrivatizeData() return 0;

proc CyclicMMAPDom.dsiPrivatize(privatizeData) {
  var privdist = chpl_getPrivatizedCopy(dist.type, dist.pid);
  return new unmanaged CyclicMMAPDom(rank, idxType, stridable,
                                 privdist, locDoms, whole);
}

proc CyclicMMAPDom.dsiGetReprivatizeData() return 0;

proc CyclicMMAPDom.dsiReprivatize(other, reprivatizeData) {
  locDoms = other.locDoms;
  whole = other.whole;
}

proc CyclicMMAPDom.dsiLocalSlice(param stridable: bool, ranges) {
  return whole((...ranges));
}


class LocCyclicMMAPDom {
  param rank: int;
  type idxType;

  // The local block type is always stridable
  // (because that's inherent to the CyclicMMAP distribution)
  var myBlock: domain(rank, idxType, stridable=true);
}

//
// Added as a performance stopgap to avoid returning a domain
//
proc LocCyclicMMAPDom.contains(i) return myBlock.contains(i);


class CyclicMMAPArr: BaseRectangularArr {
  var dom: unmanaged CyclicMMAPDom(rank, idxType, stridable);

  var locArr: [dom.dist.targetLocDom] unmanaged LocCyclicMMAPArr(eltType, rank, idxType);
  var myLocArr: unmanaged LocCyclicMMAPArr(eltType=eltType, rank=rank, idxType=idxType)?;
  const SENTINEL = max(rank*int);
}

pragma "no copy return"
proc CyclicMMAPArr.dsiLocalSlice(ranges) {
  var low: rank*idxType;
  for param i in 0..rank-1 {
    low(i) = ranges(i).alignedLow;
  }

  return locArr(dom.dist.targetLocsIdx(low)).myElems((...ranges));
}

override proc CyclicMMAPArr.dsiDisplayRepresentation() {
  for tli in dom.dist.targetLocDom {
    writeln("locArr[", tli, "].myElems = ", for e in locArr[tli].myElems do e);  }
  dom.dsiDisplayRepresentation();
}

override proc CyclicMMAPArr.dsiGetBaseDom() return dom;

override proc CyclicMMAPArr.dsiIteratorYieldsLocalElements() param {
  return true;
}


//
// NOTE: Each locale's myElems array be initialized prior to setting up
// the RAD cache.
//
proc CyclicMMAPArr.setupRADOpt() {
}

override proc CyclicMMAPArr.dsiElementInitializationComplete() {
  coforall localeIdx in dom.dist.targetLocDom {
    on dom.dist.targetLocs(localeIdx) {
      var arr = locArr(localeIdx);
    }
  }
}

override proc CyclicMMAPArr.dsiElementDeinitializationComplete() {
  coforall localeIdx in dom.dist.targetLocDom {
    on dom.dist.targetLocs(localeIdx) {
      var arr = locArr(localeIdx);
    }
  }
}

override proc CyclicMMAPArr.dsiDestroyArr(deinitElts:bool) {
  coforall localeIdx in dom.dist.targetLocDom {
    on dom.dist.targetLocs(localeIdx) {
      var arr = locArr(localeIdx);
      delete arr;
    }
  }
}

proc CyclicMMAPArr.chpl__serialize() {
  return pid;
}

proc type CyclicMMAPArr.chpl__deserialize(data) {
  return chpl_getPrivatizedCopy(unmanaged CyclicMMAPArr(rank=this.rank,
                                                    idxType=this.idxType,
                                                    stridable=this.stridable,
                                                    eltType=this.eltType),
                                data);
}

override proc CyclicMMAPArr.dsiSupportsPrivatization() param return true;

proc CyclicMMAPArr.dsiGetPrivatizeData() return 0;

proc CyclicMMAPArr.dsiPrivatize(privatizeData) {
  var privdom = chpl_getPrivatizedCopy(dom.type, dom.pid);
  var c = new unmanaged CyclicMMAPArr(eltType=eltType, rank=rank, idxType=idxType,
                              stridable=stridable, dom=privdom, locArr=locArr);
  for localeIdx in dom.dist.targetLocDom do
    if c.locArr(localeIdx).locale == here then
      c.myLocArr = c.locArr(localeIdx);
  return c;
}


inline proc _remoteAccessData.getDataIndex(
    param stridable,
    myStr: rank*chpl__signedType(idxType),
    ind: rank*idxType,
    startIdx,
    dimLen) {
  // modified from DefaultRectangularArr
  var sum = origin;
  if stridable {
    halt("RADOpt not supported for strided CyclicMMAP arrays.");
  } else {
    for param i in 0..rank-1 do {
      sum += (((ind(i) - off(i)):int * blk(i))-startIdx(i):int)/dimLen(i);
    }
  }
  return sum;
}

inline proc CyclicMMAPArr.dsiLocalAccess(i: rank*idxType) ref {
  return _to_nonnil(myLocArr).this(i);
}

proc CyclicMMAPArr.dsiAccess(i:rank*idxType) ref {
  local {
    if const myLocArrNN = myLocArr then
      if myLocArrNN.locDom.contains(i) then
        return myLocArrNN.this(i);
  }
  
  return locArr(dom.dist.targetLocsIdx(i))(i);
}

proc CyclicMMAPArr.dsiAccess(i: idxType...rank) ref
  return dsiAccess(i);

proc CyclicMMAPArr.dsiBoundsCheck(i: rank*idxType) {
  return dom.dsiMember(i);
}

iter CyclicMMAPArr.these() ref {
  foreach i in dom do
    yield dsiAccess(i);
}

iter CyclicMMAPArr.these(param tag: iterKind) where tag == iterKind.leader {
  for followThis in dom.these(tag) do
    yield followThis;
}

override proc CyclicMMAPArr.dsiStaticFastFollowCheck(type leadType) param {
  if isSubtype(leadType, CyclicMMAPArr) {
    var x : leadType?;
    return _to_borrowed(x!.dom.type) == _to_borrowed(this.dom.type);
  } else {
    return _to_borrowed(leadType) == _to_borrowed(this.dom.type);
  }
}

proc CyclicMMAPArr.dsiDynamicFastFollowCheck(lead: [])
  return this.dsiDynamicFastFollowCheck(lead.domain);

proc CyclicMMAPArr.dsiDynamicFastFollowCheck(lead: domain) {
  return lead.dist.dsiEqualDMaps(this.dom.dist) && lead._value.whole == this.dom.whole;
}

iter CyclicMMAPArr.these(param tag: iterKind, followThis, param fast: bool = false) ref where tag == iterKind.follower {
  if testFastFollowerOptimization then
    writeln((if fast then "fast" else "regular") + " follower invoked for CyclicMMAP array");

  var t: rank*range(idxType=idxType, stridable=true);
  for param i in 0..rank-1 {
    type strType = chpl__signedType(idxType);
    const wholestride = dom.whole.dim(i).stride:chpl__signedType(idxType);
    if wholestride < 0 && idxType != strType then
      halt("negative stride with unsigned idxType not supported");
    const iStride = wholestride:idxType;
    const      lo = (followThis(i).low * iStride):idxType,
               hi = (followThis(i).high * iStride):idxType,
           stride = (followThis(i).stride*wholestride):strType;
    t(i) = (lo..hi by stride) + dom.whole.dim(i).alignedLow;
  }
  const myFollowThisDom = {(...t)};
  if fast {
    const arrSection = locArr(dom.dist.targetLocsIdx(myFollowThisDom.low));

    //
    // Slicing arrSection.myElems will require reference counts to be updated.
    // If myElems is an array of arrays, the inner array's domain or dist may
    // live on a different locale and require communication for reference
    // counting. Simply put: don't slice inside a local block.
    //
    // TODO: Can myLocArr be used here to simplify things?
    //
    // MPF: Why doesn't this just slice the *domain* ?
    ref chunk = arrSection.myElems(myFollowThisDom);

    if arrSection.locale.id == here.id {
      local {
        foreach i in chunk do yield i;
      }
    } else {
      foreach i in chunk do yield i;
    }
  } else {
    proc accessHelper(i) ref {
      if const myLocArrNN = myLocArr then local {
        if myLocArrNN.locDom.contains(i) then
          return myLocArrNN.this(i);
      }
      return dsiAccess(i);
    }

    foreach i in myFollowThisDom {
      yield accessHelper(i);
    }
  }
}

proc CyclicMMAPArr.dsiSerialRead(f) {
  chpl_serialReadWriteRectangular(f, this);
}

proc CyclicMMAPArr.dsiSerialWrite(f) {
  chpl_serialReadWriteRectangular(f, this);
}

override proc CyclicMMAPArr.dsiReallocate(bounds:rank*range(idxType,BoundedRangeType.bounded,stridable)) {
  myLocArr!.dsiReallocate(bounds);
}

override proc CyclicMMAPArr.dsiPostReallocate() {
  // Call this *after* the domain has been reallocated
}

proc CyclicMMAPArr.setRADOpt(val=true) {

}

class LocCyclicMMAPArr {
  type eltType;
  param rank: int;
  type idxType;

  const locDom: unmanaged LocCyclicMMAPDom(rank, idxType);

  pragma "local field" pragma "unsafe"
  // may be initialized separately
  var myFName : string;
  var myFP : c_void_ptr;
  var myElems: c_ptr(eltType);
  var locRADLock: chpl_LocalSpinlock;

  proc init(type eltType,
            param rank: int,
            type idxType,
            const locDom: unmanaged LocCyclicMMAPDom(rank, idxType),
            param initElts: bool) {
    this.eltType = eltType;
    this.rank = rank;
    this.idxType = idxType;
    this.locDom = locDom;
    if this.locDom.myBlock.size == 0 {
      this.myElems = nil;
      this.complete();
    } else {
      var size = this.locDom.myBlock.size * c_sizeof(eltType):int;
      myFName = mmapFilePrefix + "tester-" + uniqueFileSuffix.fetchAdd(1):string + "-" +  here:string;
      myFP = fopen(myFName.c_str(), "wb+");
      assert(myFP != nil, "Failed to open: '", myFName, "'");
      var ret = ftruncate(fileno(myFP), size); // Make file exact requested size
      assert(ret == 0, "Unable to truncate file: '", myFName, "' to size of ", size, " bytes!");
      this.complete();
      var retptr = myElems:c_void_ptr;
      ret = sys_mmap(nil, size:uint, (PROT_READ | PROT_WRITE):int(32), MAP_SHARED:int(32), fileno(myFP):int(32), 0, retptr);
      myElems = retptr:c_ptr(eltType);
      assert(ret == 0, "Failed to mmap file: '", myFName, "' with error: ", ret);
    }
  }

  proc deinit() {
    // Elements in myElems are deinited in dsiDestroyArr if necessary.
    if myElems != nil then munmap(myElems, this.locDom.myBlock.size * c_sizeof(eltType):int);
    if myFP != nil {
      fclose(myFP);
      unlink(myFName.c_str());
    }
  }

  // guard against dynamic dispatch resolution trying to resolve
  // write()ing out an array of sync vars and hitting the sync var
  // type's compilerError()
  override proc writeThis(f) throws {
    halt("LocCyclicMMAPArr.writeThis() is not implemented / should not be needed");
  }

  proc dsiReallocate(bounds:rank*range) {
    if myElems != nil then munmap(myElems, this.locDom.myBlock.size * c_sizeof(eltType):int);
    if myFP != nil then fclose(myFP);
    var oldSize = this.locDom.myBlock.size * c_sizeof(eltType):int;
    this.locDom.myBlock = bounds;
    var size = this.locDom.myBlock.size * c_sizeof(eltType):int;
    myFP = fopen(myFName.c_str(), "rb+");
    assert(myFP != nil, "Failed to open: '", myFName, "'");

    use Time;

    writeln("Reallocating ", myFName, " from ", oldSize, " to ", size);
    if size > oldSize {
      var ret = fallocate(fileno(myFP), FALLOC_FL_ZERO_RANGE, oldSize, size - oldSize); // Make file exact requested size
      assert(ret == 0, "Unable to fallocate file: '", myFName, "' to size of ", size, " bytes! Error: ", ret);
      var retptr = myElems:c_void_ptr;
      ret = sys_mmap(nil, size:uint, (PROT_READ | PROT_WRITE):int(32), MAP_SHARED:int(32), fileno(myFP):int(32), 0, retptr);
      myElems = retptr:c_ptr(eltType);
      assert(ret == 0, "Failed to mmap file: '", myFName, "' with error: ", ret);
    } else {
      var ret = ftruncate(fileno(myFP), size); // Make file exact requested size
      assert(ret == 0, "Unable to truncate file: '", myFName, "' from size of ", oldSize, " bytes to ", size, " bytes! Error: ", ret);
      var retptr = myElems:c_void_ptr;
      ret = sys_mmap(nil, size:uint, (PROT_READ | PROT_WRITE):int(32), MAP_SHARED:int(32), fileno(myFP):int(32), 0, retptr);
      myElems = retptr:c_ptr(eltType);
      assert(ret == 0, "Failed to mmap file: '", myFName, "' with error: ", ret);
    }
  }
}

proc LocCyclicMMAPArr.this(i) ref where isTuple(i) {
  assert(i.size == 1, "Can only support 1-dimensional arrays!");
  return myElems[i[0]];
}

proc LocCyclicMMAPArr.this(i) ref {
  return myElems(i);
}

iter LocCyclicMMAPArr.these() ref {
  for elem in myElems {
    yield elem;
  }
}


// NOTE: I'd rather have this be a subclass of LocRADCache, but I
// couldn't find a way to use rank (param) and idxType (type) from the
// parent class in declarations in the subclass does not work.
class LocCyclicMMAPRADCache /* : LocRADCache */ {
  param rank: int;
  type idxType;
  var startIdx: rank*idxType;
  var targetLocDomDimLength: rank*int;

  proc init(param rank: int, type idxType, startIdx, targetLocDom) {
    this.rank = rank;
    this.idxType = idxType;

    this.complete();

    for param i in 0..rank-1 do
      // NOTE: Not bothering to check to see if length can fit into idxType
      targetLocDomDimLength(i) = targetLocDom.dim(i).sizeAs(int);
  }
}

private proc canDoAnyToCyclicMMAP(A, aView, B, bView) param : bool {
  return false;
}


proc CyclicMMAPArr.dsiTargetLocales() const ref {
  return dom.dist.targetLocs;
}
proc CyclicMMAPDom.dsiTargetLocales() const ref {
  return dist.targetLocs;
}
proc CyclicMMAP.dsiTargetLocales() const ref {
  return targetLocs;
}

// CyclicMMAP subdomains are represented as a single domain

proc CyclicMMAPArr.dsiHasSingleLocalSubdomain() param return true;
proc CyclicMMAPDom.dsiHasSingleLocalSubdomain() param return true;

proc CyclicMMAPArr.dsiLocalSubdomain(loc: locale) {
  if (loc == here) {
    // quick solution if we have a local array
    if const myLocArrNN = myLocArr then
      return myLocArrNN.locDom.myBlock;
    // if not, we must not own anything
    var d: domain(rank, idxType, stridable=true);
    return d;
  } else {
    return dom.dsiLocalSubdomain(loc);
  }
}
proc CyclicMMAPDom.dsiLocalSubdomain(loc: locale) {
  const (gotit, locid) = dist.chpl__locToLocIdx(loc);
  if (gotit) {
    return whole[(...(chpl__computeCyclicMMAP(this.idxType, locid, dist.targetLocDom.dims(), dist.startIdx)))];
  } else {
    var d: domain(rank, idxType, stridable=true);
    return d;
  }
}

proc newCyclicMMAPDom(dom: domain) {
  return dom dmapped CyclicMMAP(startIdx=dom.low);
}

proc newCyclicMMAPArr(dom: domain, type eltType) {
  var D = newCyclicMMAPDom(dom);
  var A: [D] eltType;
  return A;
}

proc newCyclicMMAPDom(rng: range...) {
  return newCyclicMMAPDom({(...rng)});
}

proc newCyclicMMAPArr(rng: range..., type eltType) {
  return newCyclicMMAPArr({(...rng)}, eltType);
}


proc main() : int {
  var S = {0..#1024};
  var D = S dmapped CyclicMMAP(startIdx=0);
  var A: [D] int;
  var A2 : [S] int;
  for i in D do A[i] = i;
  for i in S do A2[i] = i;
  assert((+ reduce A) == (+ reduce A2), (+ reduce A), "!=", (+ reduce A2));
  D = {0..#512};
  S = {0..#512};
  assert((+ reduce A) == (+ reduce A2), (+ reduce A), "!=", (+ reduce A2));
  D = {0..#2048};
  S = {0..#2048};
  for i in 512..2047 do A[i] = i;
  for i in 512..2047 do A2[i] = i;
  assert((+ reduce A) == (+ reduce A2), (+ reduce A), "!=", (+ reduce A2));
  return 0;
}