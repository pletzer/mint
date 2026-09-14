#include <vmtXYZCellLocator.h>
#include <mntLogger.h>
#include <sstream>


vmtXYZCellLocator::vmtXYZCellLocator() {
    this->locator = vtkStaticCellLocator::New();
    this->scratchCell = vtkGenericCell::New();
}


vmtXYZCellLocator::~vmtXYZCellLocator() {
    this->locator->Delete();
    this->scratchCell->Delete();
}


void
vmtXYZCellLocator::SetDataSet(vtkUnstructuredGrid* grid) {
    this->locator->SetDataSet(grid);
}


void
vmtXYZCellLocator::SetNumberOfCellsPerBucket(int avgNumFacesPerBucket) {
    this->locator->SetNumberOfCellsPerNode(avgNumFacesPerBucket);
}


void
vmtXYZCellLocator::BuildLocator() {
    this->locator->BuildLocator();
}


vtkIdType
vmtXYZCellLocator::FindCell(const double point[3], double tol2, vtkGenericCell *cell,
                             double pcoords[3], double *weights) {
    // vtkStaticCellLocator::FindCell needs a real cell to write into, and
    // its own signature is not const-correct about `point` even though it
    // doesn't modify it -- copy into a local, mutable array
    double p[3] = {point[0], point[1], point[2]};
    int subId;
    vtkGenericCell* genCell = cell ? cell : this->scratchCell;
    return this->locator->FindCell(p, tol2, genCell, subId, pcoords, weights);
}


void
vmtXYZCellLocator::setPeriodicityLengthX(double periodX) {
    if (periodX > 0) {
        std::stringstream msg;
        msg << "vmtXYZCellLocator does not support periodicity (periodX=" << periodX
            << " > 0 requested) -- use vmtLonLatCellLocator for a periodic-longitude grid";
        mntlog::error(__FILE__, __func__, __LINE__, msg.str());
    }
}


void
vmtXYZCellLocator::setCubedSphere(bool isCubedSphere) {
    if (isCubedSphere) {
        mntlog::error(__FILE__, __func__, __LINE__,
            "vmtXYZCellLocator does not support the cubed-sphere spherical-patch model "
            "-- use vmtLonLatCellLocator");
    }
}
