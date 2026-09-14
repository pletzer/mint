from ctypes import c_void_p, c_int, c_double, byref, POINTER
from . import MINTLIB
from . import error_handler


FILE = 'cell_locator.py'


class _VmtCellLocatorBase(object):
    """
    Shared base for LonLatCellLocator and XYZCellLocator -- both wrap a
    vmtCellLocator (the abstract C++ interface, see src/vmtCellLocator.h)
    instance and share the same setDataSet/setNumberOfCellsPerBucket/
    buildLocator operations. Construct a concrete subclass, not this class
    directly.

    Once built, hand the locator to VectorInterp.setLocator instead of
    calling VectorInterp.buildLocator -- e.g. to configure it beyond what
    buildLocator exposes, or to share one locator across several
    VectorInterp/PolylineIntegral objects. setLocator takes it as a
    borrowed reference (it does not take ownership), so keep this Python
    object alive for as long as anything is using it.
    """

    _NEW_FUNC_NAME = None  # set by subclasses

    def __init__(self):
        """
        Constructor.
        """
        self.ptr = c_void_p()
        self.obj = byref(self.ptr)

        newFunc = getattr(MINTLIB, self._NEW_FUNC_NAME)
        newFunc.argtypes = [POINTER(c_void_p)]
        ier = newFunc(self.obj)
        if ier:
            error_handler(FILE, '__init__', ier)

    def __del__(self):
        """
        Destructor.
        """
        MINTLIB.mnt_vmtcelllocator_del.argtypes = [POINTER(c_void_p)]
        ier = MINTLIB.mnt_vmtcelllocator_del(self.obj)
        if ier:
            error_handler(FILE, '__del__', ier)

    def setDataSet(self, grid):
        """
        Attach the grid this locator will search.

        :param grid: a mint.Grid instance
        """
        MINTLIB.mnt_vmtcelllocator_setDataSet.argtypes = [POINTER(c_void_p), c_void_p]
        ier = MINTLIB.mnt_vmtcelllocator_setDataSet(self.obj, grid.ptr)
        if ier:
            error_handler(FILE, 'setDataSet', ier)

    def setNumberOfCellsPerBucket(self, n):
        """
        Set the average number of cells per bucket.

        :param n: number
        """
        MINTLIB.mnt_vmtcelllocator_setNumberOfCellsPerBucket.argtypes = [POINTER(c_void_p), c_int]
        ier = MINTLIB.mnt_vmtcelllocator_setNumberOfCellsPerBucket(self.obj, n)
        if ier:
            error_handler(FILE, 'setNumberOfCellsPerBucket', ier)

    def buildLocator(self):
        """
        Build the locator. Call after setDataSet (and, for a
        LonLatCellLocator, any of setPeriodicityLengthX/enableFolding/
        setCubedSphere).
        """
        MINTLIB.mnt_vmtcelllocator_buildLocator.argtypes = [POINTER(c_void_p)]
        ier = MINTLIB.mnt_vmtcelllocator_buildLocator(self.obj)
        if ier:
            error_handler(FILE, 'buildLocator', ier)


class LonLatCellLocator(_VmtCellLocatorBase):
    """
    Cell locator for a (lon, lat[, elev=0]) grid, periodic in longitude and
    (optionally) folding across the poles -- see src/vmtLonLatCellLocator.h.
    """

    _NEW_FUNC_NAME = 'mnt_lonlatcelllocator_new'

    def setPeriodicityLengthX(self, periodX):
        """
        Set the periodicity length in x.

        :param periodX: length (0 means not periodic)
        """
        MINTLIB.mnt_vmtcelllocator_setPeriodicityLengthX.argtypes = [POINTER(c_void_p), c_double]
        ier = MINTLIB.mnt_vmtcelllocator_setPeriodicityLengthX(self.obj, periodX)
        if ier:
            error_handler(FILE, 'setPeriodicityLengthX', ier)

    def enableFolding(self):
        """
        Enable folding across the poles.
        """
        MINTLIB.mnt_vmtcelllocator_enableFolding.argtypes = [POINTER(c_void_p)]
        ier = MINTLIB.mnt_vmtcelllocator_enableFolding(self.obj)
        if ier:
            error_handler(FILE, 'enableFolding', ier)

    def setCubedSphere(self, isCubedSphere):
        """
        Declare whether the grid is a gnomonic (central-projection) cubed
        sphere, e.g. FV3's.

        :param isCubedSphere: True if the grid is a cubed sphere
        """
        MINTLIB.mnt_vmtcelllocator_setCubedSphere.argtypes = [POINTER(c_void_p), c_int]
        ier = MINTLIB.mnt_vmtcelllocator_setCubedSphere(self.obj, 1 if isCubedSphere else 0)
        if ier:
            error_handler(FILE, 'setCubedSphere', ier)


class XYZCellLocator(_VmtCellLocatorBase):
    """
    Cell locator for a genuinely 3D-embedded (x, y, z) surface mesh with no
    periodic seam -- a thin wrapper around VTK's vtkStaticCellLocator, see
    src/vmtXYZCellLocator.h. setPeriodicityLengthX/enableFolding/
    setCubedSphere make no sense for this locator and are simply not
    exposed here (unlike LonLatCellLocator) -- the underlying C++ class
    logs an error rather than silently doing the wrong thing if the shared
    (abstract-interface) versions of those calls ever reach it, e.g. via
    code that doesn't know which concrete locator it was handed.
    """

    _NEW_FUNC_NAME = 'mnt_xyzcelllocator_new'
