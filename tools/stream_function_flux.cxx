#include <mntGrid.h>
#include <mntPolysegmentIter.h>
#include <CmdLineArgParser.h>
#include <vtkUnstructuredGrid.h>
#include <vtkCellLocator.h>
#include <string>
#include <iostream>
#include <cstdio>

#define NDIMS 3

std::vector<double> parsePosition(const std::string& posStr) {
    std::vector<double> res(NDIMS, 0);
    size_t commaPos = posStr.find(',');
    res[0] = atof(posStr.substr(0, commaPos).c_str());
    res[1] = atof(posStr.substr(commaPos + 1).c_str());
    return res;
}

std::vector<double> parsePointList(const std::string& pointListStr) {

    size_t len = pointListStr.size();

    std::cerr << pointListStr << '\n';
    std::vector<double> res;
    size_t startPos = 0; 
    while (startPos < pointListStr.size()) {

        size_t lftParenPos = pointListStr.find('(', startPos);
        if (lftParenPos == std::string::npos) {
            // ( not 
            break;
        }

        size_t rgtParenPos = pointListStr.find(')', lftParenPos + 1);
        if (rgtParenPos == std::string::npos) {
            // ) not found
            std::cerr << "Warning: mismatch in parentheses?\n" << pointListStr << '\n';
            break;
        }

        // parse the position
        size_t nchars = rgtParenPos - lftParenPos - 1;
        const std::string p = pointListStr.substr(lftParenPos + 1, nchars);
        std::vector<double> point = parsePosition(p);
        for (size_t i = 0; i < point.size(); ++i) {
            res.push_back(point[i]);
        }

        // start new search on the right of )
        startPos = rgtParenPos + 1;

    }

    return res;
}

int main(int argc, char** argv) {

    int ier;
    CmdLineArgParser args;
    args.setPurpose("Project a line segment.");
    args.set("-i", std::string(""), "Source grid file in VTK format");
    args.set("-points", std::string("(0., 0.), (6.283185307179586, 0.)"), "Points");

    bool success = args.parse(argc, argv);
    bool help = args.get<bool>("-h");

    if (success && !help) {

        std::string srcFile = args.get<std::string>("-i");
        if (srcFile.size() == 0) {
            std::cerr << "ERROR: must specify a source grid file (-i)\n";
            return 1;
        }

        std::string pointListStr = args.get<std::string>("-points");
        std::vector<double> points = parsePointList(pointListStr);
        size_t npts = points.size() / NDIMS;
        std::cout << "Polyline points (lon/lat):\n";
        for (size_t i = 0; i < npts; ++i) {
            std::cout << i << " " << points[NDIMS*i + 0] << ", " << points[NDIMS*i + 1] << '\n';
        }
        if (npts < 2) {
            std::cerr << "ERROR: need at least two points to create a path\n";
            return 2;
        }

        // read/build the src grid
        Grid_t* srcGrid = NULL;
        mnt_grid_new(&srcGrid);
        mnt_grid_load(&srcGrid, srcFile.c_str());

        // get the grid pointers
        vtkUnstructuredGrid* grid = NULL;
        mnt_grid_get(&srcGrid, &grid);

        std::cout << "Grid: no of cells " << grid->GetNumberOfCells() << " no of points " << grid->GetNumberOfPoints() << '\n';

        // build locator
        vtkCellLocator* loc = vtkCellLocator::New();
        loc->SetDataSet(grid);
        loc->BuildLocator();

        // iterate over segments
        size_t nsegs = npts - 1;
        for (size_t iseg = 0; iseg < nsegs; ++iseg) {

            PolysegmentIter polyseg(grid, loc, &points[NDIMS*iseg], &points[NDIMS*(iseg + 1)]);

            size_t numSubSegs = polyseg.getNumberOfSegments();
            polyseg.reset();

            // iterate over the sub-segments 
            for (size_t isub = 0; isub < numSubSegs; ++isub) {
                vtkIdType cellId = polyseg.getCellId();
                double ta = polyseg.getBegLineParamCoord();
                double tb = polyseg.getEndLineParamCoord();
                const std::vector<double>& xia = polyseg.getBegCellParamCoord();
                const std::vector<double>& xib = polyseg.getEndCellParamCoord();
                double coeff = polyseg.getCoefficient();

                std::cout << "cell " << cellId << " t = " << ta  << " -> " << tb << 
                            " xi = " << xia[0] << ',' << xia[1] << " -> " << xib[0] << ',' << xib[1] << 
                            " coeff = " << coeff << '\n';
                polyseg.next();
            }
            double tTotal = polyseg.getIntegratedParamCoord();

        }

        // cleanup
        mnt_grid_del(&srcGrid);
        loc->Delete();


    }
    else if (help) {
        args.help();
    }
    else {
        std::cerr << "ERROR when parsing command line arguments\n";
    }

    return 0;
}