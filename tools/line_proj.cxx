#include <mntVecN.h>
#include <mntPolylineParser.h>
#include <mntGrid.h>
#include <mntPolysegmentIter.h>
#include <CmdLineArgParser.h>
#include <vtkUnstructuredGrid.h>
#include <vmtCellLocator.h>
#include <vtkCellData.h>
#include <string>
#include <iostream>
#include <cstdio>

/**
 * Compute line integral of a 2d edge field
 *
 */

int main(int argc, char** argv) {

    int ier;
    CmdLineArgParser args;
    args.setPurpose("Project a line segment.");
    args.set("-i", std::string(""), "Source grid file in VTK format");
    args.set("-p", std::string("(0., 0.),(360., 0.)"), "Points defining the path.");
    args.set("-v", std::string("edgeData"), "Edge variable name.");
    args.set("-N", 1, "Average number of cells per bucket.");
    args.set("-P", 0.0, "Periodicity length of x (0 if not periodic)");
    args.set("-verbose", false, "Verbose mode.");

    bool success = args.parse(argc, argv);
    bool help = args.get<bool>("-h");

    if (success && !help) {
        std::string srcFile = args.get<std::string>("-i");
        PolylineParser pp(2);
        pp.parse(args.get<std::string>("-p"));
        const std::vector< Vec3 >& points = pp.getPoints();
        size_t npts = points.size();
        std::cout << "Path:\n";
        for (size_t i = 0; i< npts; ++i) {
            std::cout << i << ": " << points[i] << '\n';
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
        vmtCellLocator* loc = vmtCellLocator::New();
        loc->SetDataSet(grid);
        loc->SetNumberOfCellsPerBucket(args.get<int>("-N"));
        loc->BuildLocator();

        double totFlux = 0.0;
        vtkDataArray* arr = grid->GetCellData()->GetArray(args.get<std::string>("-v").c_str());
        if (! arr) {
            std::cerr << "ERROR: could not find edge field \"" << args.get<std::string>("-v") << "\"\n";
            return 4;
        }

        // iterate over segments
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

            PolysegmentIter polyseg(grid, loc, &points[iseg0][0], &points[iseg1][0], args.get<double>("-P"));

            double fluxFromSegment = 0;
            size_t numSegs = polyseg.getNumberOfSegments();
            polyseg.reset();
            for (size_t i = 0; i < numSegs; ++i) {

                vtkIdType cellId = polyseg.getCellId();
                double ta = polyseg.getBegLineParamCoord();
                double tb = polyseg.getEndLineParamCoord();
                const Vec3& xia = polyseg.getBegCellParamCoord();
                const Vec3& xib = polyseg.getEndCellParamCoord();
                double coeff = polyseg.getCoefficient();

                if (args.get<bool>("-verbose")) {
                    std::cout << "\tsub segment: " << i << " cell=" << cellId << " coeff=" << coeff << " ta=" << ta << " tb=" << tb 
                              << " xia=" << xia[0] << "," << xia[1] << " xib=" << xib[0] << "," << xib[1] << '\n';
                }

                std::vector<double> dxi({xib[0] - xia[0], xib[1] - xia[1]});
                std::vector<double> xiMid({0.5*(xia[0] + xib[0]), 0.5*(xia[1] + xib[1])});

                // interpolation weights
                double ws[] = {+ dxi[0] * (1.0 - xiMid[1]) * coeff,
                               + dxi[1] * (0.0 + xiMid[0]) * coeff,
                               + dxi[0] * (0.0 + xiMid[1]) * coeff,
                               + dxi[1] * (1.0 - xiMid[0]) * coeff};

                double* fieldVals = arr->GetTuple(cellId);

                if (args.get<bool>("-verbose")) {
                    std::cout << "\t             weights        : " << ws[0] << "," << ws[1] << "," << ws[2] << "," << ws[3] << '\n';
                    std::cout << "\t             edge field vals: " << fieldVals[0] << "," << fieldVals[1] << "," << fieldVals[2] << "," << fieldVals[3] << '\n';
                }

                // iterate over edges
                for (size_t ie = 0; ie < 4; ++ie) {
                    fluxFromSegment += ws[ie] * fieldVals[ie];
                }

                polyseg.next();
            }

            // check if the line segment is fully accounted for
            double tTotal = polyseg.getIntegratedParamCoord();
            const double tol = 1.e-10;
            if (std::abs(tTotal - 1.0) > tol) {
                std::cout << "Warning: total integrated length for segment " << iseg0 << " is " << tTotal << " != 1 (diff=" << tTotal - 1. << ")\n";
            }

            // add this cell contribution
            totFlux += fluxFromSegment;
        }

        std::cout << "total flux is: " << totFlux << '\n';

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
