#include <mntPolylineParser.h>
#include <mntGrid.h>
#include <mntPolysegmentIter3d.h>
#include <CmdLineArgParser.h>
#include <vtkUnstructuredGrid.h>
#include <vtkCellLocator.h>
#include <vtkCellData.h>
#include <string>
#include <iostream>
#include <cstdio>


/**
 * Compute line integral of a 3d edge field
 *
 */

int main(int argc, char** argv) {

    int ier;
    CmdLineArgParser args;
    args.setPurpose("Project a line segment.");
    args.set("-i", std::string(""), "Source grid file in VTK format");
    args.set("-p", std::string("(0., 0.),(6.283185307179586, 0.)"), "Points defining the path.");
    args.set("-v", std::string("edgeData"), "Edge variable name.");
    args.set("-N", 128, "Average number of cells per bucket.");
    args.set("-verbose", false, "Verbose mode.");

    bool success = args.parse(argc, argv);
    bool help = args.get<bool>("-h");

    if (success && !help) {
        std::string srcFile = args.get<std::string>("-i");
        PolylineParser pp(3); // 3d
        pp.parse(args.get<std::string>("-p"));
        const std::vector< Vector<double> >& points = pp.getPoints();
        size_t npts = points.size();
        std::cout << "Path:\n";
        for (size_t i = 0; i< npts; ++i) {
            std::cout << i << ": ";
            for (size_t j = 0; j < points[i].size(); ++j) std::cout << points[i][j] << ", ";
            std::cout << '\n';
        }

        // checks
        if (srcFile.size() == 0) {
            std::cerr << "ERROR: must specify a source grid file (-i)\n";
            return 1;
        }
        if (npts < 2) {
            std::cerr << "ERROR: must have at least two points\n";
            return 2;
        }
        // compute length of path and check it is != 0
        double pathLengthSq = 0;
        for (size_t i = 0; i < npts - 1; ++i) {
            double dx = points[i + 1][0] - points[i + 0][0];
            double dy = points[i + 1][1] - points[i + 0][1];
            double dz = points[i + 1][2] - points[i + 0][2];
            pathLengthSq += dx*dx + dy*dy + dz*dz;
        }
        if (pathLengthSq == 0) {
            std::cerr << "ERROR: path length must be > 0\n";
            return 3;
        }

        // read/build the src grid
        Grid_t* srcGrid = NULL;
        mnt_grid_new(&srcGrid);
        mnt_grid_load(&srcGrid, srcFile.c_str());

        // get the grid pointers
        vtkUnstructuredGrid* grid = NULL;
        mnt_grid_get(&srcGrid, &grid);

        std::cout << "no of cells " << grid->GetNumberOfCells() << " no of points " << grid->GetNumberOfPoints() << '\n';

        // build locator
        vtkCellLocator* loc = vtkCellLocator::New();
        loc->SetDataSet(grid);
        loc->SetNumberOfCellsPerBucket(args.get<int>("-N"));
        loc->BuildLocator();

        vtkDataArray* arr = grid->GetCellData()->GetArray(args.get<std::string>("-v").c_str());

        // iterate over segments
        double lineIntegral = 0.0;
        size_t nsegs = npts - 1;
        for (size_t iseg0 = 0; iseg0 < nsegs; ++iseg0) {

            size_t iseg1 = iseg0 + 1;
            std::cout << "Segment " << iseg0 << " (";
            for (size_t j = 0; j < points[iseg0].size(); ++j) {
                std::cout << points[iseg0][j] << ", ";
            }
            std::cout << ") -> (";
            for (size_t j = 0; j < points[iseg1].size(); ++j) {
                std::cout << points[iseg1][j] << ", ";
            }
            std::cout << ")\n";

            PolysegmentIter3d polyseg(grid, loc, &points[iseg0][0], &points[iseg1][0]);

            double integralFromSegment = 0;
            size_t numSegs = polyseg.getNumberOfSegments();
            polyseg.reset();
            for (size_t i = 0; i < numSegs; ++i) {

                vtkIdType cellId = polyseg.getCellId();
                double ta = polyseg.getBegLineParamCoord();
                double tb = polyseg.getEndLineParamCoord();
                const Vector<double>& xia = polyseg.getBegCellParamCoord();
                const Vector<double>& xib = polyseg.getEndCellParamCoord();
                double coeff = polyseg.getCoefficient();

                if (args.get<bool>("-verbose")) {
                    std::cout << "\tsub segment: " << i << " cell=" << cellId << " coeff=" << coeff << " ta=" << ta << " tb=" << tb 
                              << " xia = " << xia << " xib = " << xib << '\n';
                }

                // difference from start to end
                Vector<double> dxi = xib - xia;
                // mid point of the line
                Vector<double> xiMid = 0.5*(xia + xib);

                // this cell
                vtkCell* cell = grid->GetCell(cellId);

                double pointEdge0[3];
                double pointEdge1[3];
                double paramEdge0[3];
                double paramEdge1[3];
                int subId;
                double dist2;
                double weights[8];
                vtkIdType pointId0, pointId1;
                vtkPoints* gridPoints = grid->GetPoints();

                double* fieldVals = arr->GetTuple(cellId);

                // iterate over the edges
                for (int iEdge = 0; iEdge < cell->GetNumberOfEdges(); ++iEdge) {

                    vtkCell* edge = cell->GetEdge(iEdge);
                    pointId0 = edge->GetPointId(0);
                    pointId1 = edge->GetPointId(1);
                    gridPoints->GetPoint(pointId0, pointEdge0);
                    gridPoints->GetPoint(pointId1, pointEdge1);
                    // get the parametric edge positions
                    edge->EvaluatePosition(pointEdge0, NULL, subId, paramEdge0, dist2, weights);
                    edge->EvaluatePosition(pointEdge1, NULL, subId, paramEdge1, dist2, weights);

                    // compute the projection of the line onto the edge basis function
                    double basisIntegral = 1.0;
                    for ( size_t j = 0; j < 3; ++j) {
                        double paramEdgeMid = 0.5*(paramEdge0[j] + paramEdge1[j]);
                        double paramEdgeMidBar = 1.0 - paramEdgeMid;
                        double xiM = xiMid[j];
                        double xiMBar = 1.0 - xiM;
                        basisIntegral *= paramEdgeMidBar*xiMBar + paramEdgeMid*xiM + 4.0*paramEdgeMidBar*paramEdgeMid*dxi[j];
                    }

                    lineIntegral += fieldVals[iEdge] * basisIntegral;

                    if (args.get<bool>("-verbose")) {
                        std::cout << "\t\tedge: " << iEdge << " weight: " << basisIntegral << " field: " << fieldVals[iEdge] << '\n';
                    }

                }

                polyseg.next();
            }

            // check if the line segment is fully accounted for
            double tTotal = polyseg.getIntegratedParamCoord();
            const double tol = 1.e-10;
            if (std::abs(tTotal - 1.0) > tol) {
                std::cout << "Warning: total integrated length for segment " << iseg0 << " is " << tTotal << " != 1 (diff=" << tTotal - 1. << ")\n";
            }
        }

        std::cout << "total line integral is: " << lineIntegral << '\n';

        // cleanup
        mnt_grid_del(&srcGrid);
        loc->Delete();


    }
    else if (help) {
        args.help();
    }
    else {
        std::cerr << "ERROR when parsing command line arguments\n";
        return 3;
    }

    return 0;
}
