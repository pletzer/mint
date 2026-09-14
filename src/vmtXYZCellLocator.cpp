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
    // doesn't modify it -- copy into a local, mutable array. Deliberately
    // using the 5-argument overload (no subId): some VTK versions (e.g.
    // 9.1, as shipped by several Linux distros) only expose that one on
    // vtkStaticCellLocator itself -- declaring any FindCell override hides
    // the base class's other overloads (including the 6-argument one with
    // subId) unless re-exposed with a `using` declaration, which VTK 9.1
    // does not do here. The 5-argument form is present in every VTK
    // version this library supports, so use that one instead of relying
    // on subId being available.
    double p[3] = {point[0], point[1], point[2]};
    vtkGenericCell* genCell = cell ? cell : this->scratchCell;
    return this->locator->FindCell(p, tol2, genCell, pcoords, weights);
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
