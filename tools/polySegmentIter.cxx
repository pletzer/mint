#include <mntGrid.h>
#include <CmdLineArgParser.h>
#include <mntVecN.h>
#include <vmtCellLocator.h>
#include <mntPolysegmentIter.h>
#include <vtkUnstructuredGrid.h>


Vec3 getVectorFromString(const std::string& expr) {
    Vec3 res(0.);
    size_t commaPos = expr.find(',');
    if (commaPos != std::string::npos) {
        res[0] = std::stod(expr.substr(0, commaPos));
        res[1] = std::stod(expr.substr(commaPos + 1, std::string::npos));
    }
    return res;
}

int main(int argc, char** argv) {

    int ier;
    CmdLineArgParser args;
    args.setPurpose("Break a segment into subsegments.");
    args.set("-s", std::string(""), "UGRID source grid file and mesh name, specified as \"filename:meshname\"");
    args.set("-p0", std::string("0., 0."), "lon,lat start point");
    args.set("-p1", std::string("360., 0."), "lon,lat end point");
    args.set("-S", 1, "Set to zero to disable source grid regularization, -S 0 is required for uniform lon-lat grid");

    bool success = args.parse(argc, argv);
    bool help = args.get<bool>("-h");

    if (!success) {
        std::cout << help << '\n';
        // error
        return 1;
    }

    std::string srcFileMesh = args.get<std::string>("-s");

    Grid_t* grd;
    ier = mnt_grid_new(&grd);
    ier = mnt_grid_loadFrom2DUgrid(&grd, srcFileMesh.c_str());
    if (ier != 0) {
        std::cerr << "ERROR: failed to read " << srcFileMesh << '\n';
        return 2;        
    }

    // defaults are suitable for cubed-sphere 
    int fixLonAcrossDateline = 1;
    int averageLonAtPole = 1;
    if (args.get<int>("-S") == 0) {
        fixLonAcrossDateline = 0;
        averageLonAtPole = 0;
        std::cout << "info: no regularization applied to grid\n";
    }
    ier = mnt_grid_setFlags(&grd, fixLonAcrossDateline, averageLonAtPole);


    vtkUnstructuredGrid* ugrid;
    ier = mnt_grid_get(&grd, &ugrid);

    if (args.get<std::string>("-p0").find(',') == std::string::npos) {
        std::cerr << "ERROR: must provide -p0 lon,lat\n";
        return 3;
    }

    if (args.get<std::string>("-p1").find(',') == std::string::npos) {
        std::cerr << "ERROR: must provide -p1 lon,lat\n";
        return 4;
    }

    Vec3 p0 = getVectorFromString(args.get<std::string>("-p0"));
    std::cout << "start point: " << p0 << '\n';
    Vec3 p1 = getVectorFromString(args.get<std::string>("-p1"));
    std::cout << "end   point: " << p1 << '\n';

    vmtCellLocator* loc = vmtCellLocator::New();
    loc->SetDataSet(ugrid);
    loc->BuildLocator();

    PolysegmentIter psi(ugrid, loc, &p0[0], &p1[0]);
    size_t numSegs = psi.getNumberOfSegments();
    psi.reset();
    for (size_t i = 0; i < numSegs; ++i) {
        vtkIdType cellId = psi.getCellId();
        const Vec3& xia = psi.getBegCellParamCoord();
        const Vec3& xib = psi.getEndCellParamCoord();
        double ta = psi.getBegLineParamCoord();
        double tb = psi.getEndLineParamCoord();
        double coeff = psi.getCoefficient();
        std::cout << "seg " << i << " cell=" << cellId \
                            << " ta=" << ta << " xia=" << xia[0] << ',' << xia[1] 
                            << " tb=" << tb << " xib=" << xib[0] << ',' << xib[1] 
                            << '\n';
        psi.next();
    }
    std::cout << "num segments = " << numSegs << " integrated param coord = " << psi.getIntegratedParamCoord() << '\n';

    // cleanup
    loc->Delete();
    ier = mnt_grid_del(&grd);

    return 0;
}