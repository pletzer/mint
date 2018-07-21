#include <mntGrid.h>
#include <mntPolysegmentIter.h>
#include <CmdLineArgParser.h>
#include <vtkUnstructuredGrid.h>
#include <vtkCellLocator.h>
#include <vtkCellData.h>
#include <string>
#include <iostream>
#include <cstdio>

std::vector<double> parsePosition(const std::string& posStr) {
    std::vector<double> res(3, 0);
    size_t commaPos = posStr.find(',');
    res[0] = atof(posStr.substr(0, commaPos).c_str());
    res[1] = atof(posStr.substr(commaPos + 1).c_str());
    return res;
}

std::vector< std::vector<double> > parsePoints(const std::string& pointsStr) {
    std::vector< std::vector<double> > res;
    const char leftDelim = '(';
    const char rghtDelim = ')';
    size_t leftPos = 0;
    size_t rghtPos = pointsStr.size();
    while(true) {
        leftPos = pointsStr.find(leftDelim, leftPos);
        rghtPos = pointsStr.find(rghtDelim, leftPos);
        if (leftPos == std::string::npos || rghtPos == std::string::npos) {
            // could not find a point or parentheses don't match
            break;
        }
        leftPos++;
        size_t n = rghtPos - leftPos;
        std::string posStr = pointsStr.substr(leftPos, n);
        res.push_back(parsePosition(posStr));
    }
    return res;
}

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
        std::vector< std::vector<double> > points = parsePoints(args.get<std::string>("-p"));
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
        double pathLength = 0;
        for (size_t i = 0; i < npts - 1; ++i) {
            double dx = points[i + 1][0] - points[i + 0][0];
            double dy = points[i + 1][1] - points[i + 0][1];
            double dz = points[i + 1][2] - points[i + 0][2];
            pathLength += sqrt(dx*dx + dy*dy + dz*dz);
        }
        if (pathLength == 0) {
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

        double totFlux = 0.0;
        vtkDataArray* arr = grid->GetCellData()->GetArray(args.get<std::string>("-v").c_str());

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

            PolysegmentIter polyseg(grid, loc, &points[iseg0][0], &points[iseg1][0]);

            double fluxFromSegment = 0;
            size_t numSegs = polyseg.getNumberOfSegments();
            polyseg.reset();
            for (size_t i = 0; i < numSegs; ++i) {

                vtkIdType cellId = polyseg.getCellId();
                double ta = polyseg.getBegLineParamCoord();
                double tb = polyseg.getEndLineParamCoord();
                const std::vector<double>& xia = polyseg.getBegCellParamCoord();
                const std::vector<double>& xib = polyseg.getEndCellParamCoord();
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
