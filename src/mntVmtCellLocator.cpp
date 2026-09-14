#include <mntVmtCellLocator.h>
#include <vmtLonLatCellLocator.h>
#include <vmtXYZCellLocator.h>
#include <mntLogger.h>
#include <string>


LIBRARY_API
int mnt_lonlatcelllocator_new(vmtCellLocator** self) {
    *self = vmtLonLatCellLocator::New();
    return 0;
}


LIBRARY_API
int mnt_xyzcelllocator_new(vmtCellLocator** self) {
    *self = vmtXYZCellLocator::New();
    return 0;
}


LIBRARY_API
int mnt_vmtcelllocator_del(vmtCellLocator** self) {
    if (!*self) {
        std::string msg = "locator is NULL, cannot delete";
        mntlog::error(__FILE__, __func__, __LINE__, msg);
        return -1;
    }
    (*self)->Delete();
    *self = nullptr;
    return 0;
}


LIBRARY_API
int mnt_vmtcelllocator_setDataSet(vmtCellLocator** self, Grid_t* grid) {
    if (!*self) {
        std::string msg = "must construct the locator before calling this";
        mntlog::error(__FILE__, __func__, __LINE__, msg);
        return -1;
    }
    if (!grid) {
        std::string msg = "grid is NULL";
        mntlog::error(__FILE__, __func__, __LINE__, msg);
        return -2;
    }
    (*self)->SetDataSet(grid->grid);
    return 0;
}


LIBRARY_API
int mnt_vmtcelllocator_setNumberOfCellsPerBucket(vmtCellLocator** self, int n) {
    if (!*self) {
        std::string msg = "must construct the locator before calling this";
        mntlog::error(__FILE__, __func__, __LINE__, msg);
        return -1;
    }
    (*self)->SetNumberOfCellsPerBucket(n);
    return 0;
}


LIBRARY_API
int mnt_vmtcelllocator_setPeriodicityLengthX(vmtCellLocator** self, double periodX) {
    if (!*self) {
        std::string msg = "must construct the locator before calling this";
        mntlog::error(__FILE__, __func__, __LINE__, msg);
        return -1;
    }
    (*self)->setPeriodicityLengthX(periodX);
    return 0;
}


LIBRARY_API
int mnt_vmtcelllocator_enableFolding(vmtCellLocator** self) {
    if (!*self) {
        std::string msg = "must construct the locator before calling this";
        mntlog::error(__FILE__, __func__, __LINE__, msg);
        return -1;
    }
    (*self)->enableFolding();
    return 0;
}


LIBRARY_API
int mnt_vmtcelllocator_setCubedSphere(vmtCellLocator** self, int isCubedSphere) {
    if (!*self) {
        std::string msg = "must construct the locator before calling this";
        mntlog::error(__FILE__, __func__, __LINE__, msg);
        return -1;
    }
    (*self)->setCubedSphere(isCubedSphere != 0);
    return 0;
}


LIBRARY_API
int mnt_vmtcelllocator_buildLocator(vmtCellLocator** self) {
    if (!*self) {
        std::string msg = "must construct the locator before calling this";
        mntlog::error(__FILE__, __func__, __LINE__, msg);
        return -1;
    }
    (*self)->BuildLocator();
    return 0;
}
