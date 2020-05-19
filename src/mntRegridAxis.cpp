#include <mntRegridAxis.h>

extern "C"
int mnt_regridaxis_new(RegridAxis_t** self) {

    *self = new RegridAxis_t();

    (*self)->spline = vtkCardinalSpline::New();

    // let the the spline estimate the second derivative at the boundaries
    (*self)->spline->SetLeftConstraint(3);
    (*self)->spline->SetRightConstraint(3);

    return 0;
}

extern "C"
int mnt_regridaxis_del(RegridAxis_t** self) {

    // destroy the spline object
    (*self)->spline->Delete();

    delete *self;

    return 0;
}


extern "C"
int mnt_regridaxis_build(RegridAxis_t** self, int numValues, const double srcValues[]) {

    (*self)->spline->RemoveAllPoints();

    if (numValues < 2) {
        // need at least two values (one interval)
        return 1;
    }

    (*self)->numValues = numValues;
    (*self)->tLo = 0;
    (*self)->tHi = numValues - 1;

    // axis is monotonically increasing
    (*self)->increasing = true;
    if (srcValues[numValues - 1] < srcValues[0]) {
        // axis is monotonically decreasing
        (*self)->increasing = false;
    }

    // construct the spline object. Note that the spline object maps the 
    // axis values to indices. This will allow us to quickly find the 
    // float index of a target axis value.
    int ier = 0;
    for (int i = 0; i < numValues - 1; ++i) {

        (*self)->spline->AddPoint(srcValues[i], double(i));

        if ((*self)->increasing && srcValues[i + 1] < srcValues[i]) {
            ier++;
        }
        else if (!(*self)->increasing && srcValues[i + 1] > srcValues[i]) {
            ier++;
        }
    }
    // add last value
    (*self)->spline->AddPoint(srcValues[numValues - 1], double(numValues - 1));

    // compute the spline coefficients
    (*self)->spline->Compute();

    return ier;
}

extern "C"
int mnt_regridaxis_getPointWeights(RegridAxis_t** self, double target, int indices[2], double weights[2]) {

    double t = (*self)->spline->Evaluate(target);

    // low side of the interval
    indices[0] = std::max(0, std::min((*self)->numValues - 2, int(floor(t))));

    // high side of the interval
    indices[1] = indices[0] + 1;

    // linear interpolation weights
    weights[0] = indices[1] - t;
    weights[1] = t - indices[0];

    return 0;
}


extern "C"
int mnt_regridaxis_getNumCellWeights(RegridAxis_t** self, const double targets[2], int* numCellWeights) {

    (*self)->tA = (*self)->spline->Evaluate(targets[0]);
    (*self)->tB = (*self)->spline->Evaluate(targets[1]);

    (*self)->indicesA[0] = std::max(0, std::min((*self)->numValues - 2, int(floor((*self)->tA))));
    (*self)->indicesA[1] = (*self)->indicesA[0] + 1;

    (*self)->indicesB[0] = std::max(0, std::min((*self)->numValues - 2, int(floor((*self)->tB))));
    (*self)->indicesB[1] = (*self)->indicesA[0] + 1;

    (*self)->numCellWeights = (*self)->indicesB[0] - (*self)->indicesA[1] + 2;

    *numCellWeights = (*self)->numCellWeights;

    return 0;
}

extern "C"
int mnt_regridaxis_getCellWeights(RegridAxis_t** self, int indices[], double weights[]) {

    // right weight on the low side
    indices[0] = (*self)->indicesA[0];
    weights[0] = (*self)->tA - (*self)->indicesA[0];

    // wwight of one for all the fully covered intervals
    for (int i = 1; i < (*self)->numCellWeights - 2; ++i) {
        indices[i] = (*self)->indicesA[0] + i;
        weights[i] = 1.0;
    }

    // left weight for the right side
    indices[(*self)->numCellWeights - 1] = (*self)->indicesB[0];
    weights[(*self)->numCellWeights - 1] = (*self)->indicesB[1] - (*self)->tB;

    return 0;
}

