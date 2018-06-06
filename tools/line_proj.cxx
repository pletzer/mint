#include <mntGrid.h>
#include <mntPolysegmentIter.h>
#include <CmdLineArgParser.h>
#include <vtkUnstructuredGrid.h>
#include <vtkCellLocator.h>
#include <string>
#include <iostream>

std::vector<double> parsePosition(const std::string& posStr) {
	std::vector<double> res(2);
	size_t commaPos;
	res[0] = std::stod(posStr, &commaPos);
	res[1] = std::stod(posStr.substr(commaPos + 1));
	return res;
}

int main(int argc, char** argv) {

    int ier;
    CmdLineArgParser args;
    args.setPurpose("Project a line segment.");
    args.set("-i", std::string(""), "Source grid file in VTK format");
    args.set("-p0", std::string("0., 0."), "Start point");
    args.set("-p1", std::string("2*pi, 0."), "End point");

    bool success = args.parse(argc, argv);
    bool help = args.get<bool>("-h");

    if (success && !help) {
        std::string srcFile = args.get<std::string>("-i");
        std::vector<double> p0 = parsePosition(args.get<std::string>("-p0"));
        std::vector<double> p1 = parsePosition(args.get<std::string>("-p1"));
        std::cout << "p0 = " << p0[0] << ", " << p0[1] << '\n';
        std::cout << "p1 = " << p1[0] << ", " << p1[1] << '\n';

        if (srcFile.size() == 0) {
            std::cerr << "ERROR: must specify a source grid file (-i)\n";
            return 1;
        }


        // read/build the src grid
        Grid_t* srcGrid = NULL;
        mnt_grid_new(&srcGrid);
        mnt_grid_load(&srcGrid, srcFile.c_str());

        // get the grid pointers
        vtkUnstructuredGrid* grid = NULL;
        mnt_grid_get(&srcGrid, &grid);

        // build locator
        vtkCellLocator* loc = vtkCellLocator::New();
        loc->SetDataSet(grid);

        PolysegmentIter polyseg(grid, loc, &p0[0], &p1[0]);

        size_t numSegs = polyseg.getNumberOfSegments();
        polyseg.reset();
        for (size_t i = 0; i < numSegs; ++i) {
        	vtkIdType cellId = polyseg.getCellId();
        	double ta = polyseg.getBegLineParamCoord();
        	double tb = polyseg.getEndLineParamCoord();
        	const std::vector<double>& xia = polyseg.getBegCellParamCoord();
        	const std::vector<double>& xib = polyseg.getEndCellParamCoord();
        	double coeff = polyseg.getCoefficient();
        	polyseg.next();
        }
        double tTotal = polyseg.getIntegratedParamCoord();

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