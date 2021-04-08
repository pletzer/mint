#include <mntRegridEdges.h>
#include <mntNcAttributes.h>
#include <mntNcDimensions.h>
#include <mntMultiArrayIter.h>
#include <mntNcFieldRead.h>
#include <mntNcFieldWrite.h>
#include <mntGrid.h>
#include <CmdLineArgParser.h>
#include <vtkUnstructuredGrid.h>
#include <vtkAbstractArray.h>
#include <vtkCellData.h>
#include <iostream>
#include <limits>
#include <string>
#include <sstream>
#include <cmath>
#include <regex>
#include <netcdf.h>

#define NUM_EDGES_PER_CELL 4

/**
 * Convert a number to a string
 * @param arg number
 * @return string
 */
template <typename T>
std::string toString (T arg) {
    std::stringstream ss;
    ss << arg;
    return ss.str ();
}

/**
 * Split string by separator
 * @param fmname string, eg "filename:meshname"
 * @param separator separator, eg ':'
 * @return filename, meshname pair
 */
std::pair<std::string, std::string> split(const std::string& fmname, char separator) {
    std::pair<std::string, std::string> res;
    size_t pos = fmname.find(separator);
    res.first = fmname.substr(0, pos);
    if (pos < std::string::npos) {
        // file name substring
        // mesh name substring
        res.second = fmname.substr(pos + 1, std::string::npos);
    }
    return res;
}

/**
 * Compute the loop integrals for each cell
 * @param gridObj grid object
 * @param edge data
 * @param avgAbsLoop abs of average loop integral (output)
 * @param minAbsLoop abs of min loop integral (output)
 * @param maxAbsLoop abs of max loop integral (output)
 * @param loop_integrals loop integrals, array must have number of cells elements
 */
void computeLoopIntegrals(Grid_t* grd, const std::vector<double>& edgeData,
                          double* avgAbsLoop, double* minAbsLoop, double* maxAbsLoop,
                          std::vector<double>& loop_integrals) {
    size_t numCells, edgeId;
    int edgeSign, ier;
    mnt_grid_getNumberOfCells(&grd, &numCells);
    *minAbsLoop = + std::numeric_limits<double>::max();
    *maxAbsLoop = - std::numeric_limits<double>::max();
    *avgAbsLoop = 0.0;
    for (size_t cellId = 0; cellId < numCells; ++cellId) {
        double loop = 0.0;
        for (int ie = 0; ie < NUM_EDGES_PER_CELL; ++ie) {

            ier = mnt_grid_getEdgeId(&grd, cellId, ie, &edgeId, &edgeSign);
            assert(ier == 0);

            // +1 for ie = 0, 1; -1 for ie = 2, 3
            int sgn = 1 - 2*(ie/2);
            
            loop += sgn * edgeSign * edgeData[edgeId];
        }

        loop_integrals[cellId] = loop;
        loop = std::abs(loop);
        *minAbsLoop = std::min(loop, *minAbsLoop);
        *maxAbsLoop = std::max(loop, *maxAbsLoop);
        *avgAbsLoop += loop;
    }
    *avgAbsLoop /= double(numCells);
}

int main(int argc, char** argv) {

    int ier;
    NcFieldRead_t* reader = NULL;
    NcFieldWrite_t* writer = NULL;

    CmdLineArgParser args;
    args.setPurpose("Regrid an edge centred field.");
    args.set("-s", std::string(""), "UGRID source grid file and mesh name, specified as \"filename:meshname\"");
    args.set("-v", std::string(""), "Specify edge staggered field variable name in source UGRID file, varname[@filename:meshname]");
    args.set("-P", 0.0, "Specify the periodicity length in longitudes (default is non-periodic)");
    args.set("-d", std::string(""), "UGRID destination grid file name");
    args.set("-w", std::string(""), "Write interpolation weights to file");
    args.set("-W", std::string(""), "Load interpolation weights from file");
    args.set("-o", std::string(""), "Specify output VTK file where regridded edge data are saved");
    args.set("-O", std::string(""), "Specify output 2D UGRID file where regridded edge data are saved");
    args.set("-S", 1, "Set to zero to disable source grid regularization, -S 0 is required for uniform lon-lat grid");
    args.set("-D", 1, "Set to zero to disable destination grid regularization, -S 0 is required for uniform lon-lat grid");
    args.set("-N", 128, "Average number of cells per bucket");
    args.set("-debug", 1, "0=no checks, 1=print outside segments, 2=save outside segments");

    bool success = args.parse(argc, argv);
    bool help = args.get<bool>("-h");

    if (help) {
        args.help();
        return 0;
    }

    if (!success) {
        std::cerr << "ERROR when parsing command line arguments\n";
        return 1;
    }

    // extract the command line arguments
    std::string srcFile = args.get<std::string>("-s");
    std::string dstFile = args.get<std::string>("-d");
    std::string weightsFile = args.get<std::string>("-w");
    std::string loadWeightsFile = args.get<std::string>("-W");
    std::string vtkOutputFile = args.get<std::string>("-o");
    std::string dstEdgeDataFile = args.get<std::string>("-O");

    // run some checks
    if (srcFile.size() == 0) {
        std::cerr << "ERROR: must specify a source grid file (-s)\n";
        return 2;
    }
    if (dstFile.size() == 0) {
        std::cerr << "ERROR: must specify a destination grid file (-d)\n";
        return 3;
    }

    // create regridder 
    RegridEdges_t* rg;
    mnt_regridedges_new(&rg);

    // defaults are suitable for cubed-sphere 
    int fixLonAcrossDateline = 1;
    int averageLonAtPole = 1;
    if (args.get<int>("-S") == 0) {
        fixLonAcrossDateline = 0;
        averageLonAtPole = 0;
        std::cout << "info: no regularization applied to source grid\n";
    }
    ier = mnt_regridedges_setSrcGridFlags(&rg, fixLonAcrossDateline, averageLonAtPole);

    // ...destination grid
    fixLonAcrossDateline = 1;
    averageLonAtPole = 1;
    if (args.get<int>("-D") == 0) {
        fixLonAcrossDateline = 0;
        averageLonAtPole = 0;
        std::cout << "info: no regularization applied to destination grid\n";
    }
    ier = mnt_regridedges_setDstGridFlags(&rg, fixLonAcrossDateline, averageLonAtPole);

    // read the source grid
    ier = mnt_regridedges_loadSrcGrid(&rg, srcFile.c_str(), srcFile.size());
    if (ier != 0) {
        std::cerr << "ERROR: could not read file \"" << srcFile << "\"\n";
        return 4;
    }

    // read the destination grid
    ier = mnt_regridedges_loadDstGrid(&rg, dstFile.c_str(), dstFile.size());
    if (ier != 0) {
        std::cerr << "ERROR: could not read file \"" << dstFile << "\"\n";
        return 5;
    }

    if (loadWeightsFile.size() == 0) {

        // compute the weights
        std::cout << "info: computing weights\n";
        ier = mnt_regridedges_build(&rg, args.get<int>("-N"), 
                                         args.get<double>("-P"), args.get<int>("-debug"));
        if (ier != 0) {
            return 6;
        }
    
        // save the weights to file
        if (weightsFile.size() != 0) {
            std::cout << "info: saving weights in file " << weightsFile << '\n';
            ier = mnt_regridedges_dumpWeights(&rg, weightsFile.c_str(), (int) weightsFile.size());
        }

    }
    else {

        // weights have been pre-computed, just load them
        std::cout << "info: loading weights from file " << loadWeightsFile << '\n';
        ier = mnt_regridedges_loadWeights(&rg, loadWeightsFile.c_str(), (int) loadWeightsFile.size());
        if (ier != 0) {
            return 7;
        }

    }

    std::string varAtFileMesh = args.get<std::string>("-v");
    if (varAtFileMesh.size() > 0) {

        // get the variable name and the source file/mesh names
        std::pair<std::string, std::string> vfm = split(varAtFileMesh, '@');
        std::string vname = vfm.first;
        // by default the variable is stored in srcFile
        std::string srcFileMeshName = srcFile;
        if (vfm.second.size() > 0) {
            srcFileMeshName = vfm.second;
        }
        std::pair<std::string, std::string> fm = split(srcFileMeshName, ':');
        std::string srcFileName = fm.first;

        // get the ncid and varid's so we can read the attributes and dimensions
        int srcNcid;
        ier = nc_open(srcFileName.c_str(), NC_NOWRITE, &srcNcid);
        if (ier != 0) {
            std::cerr << "ERROR: could not open file \"" << srcFileName << "\"\n";
            return 8;
        }

        int srcVarid;
        ier = nc_inq_varid(srcNcid, vname.c_str(), &srcVarid);
        if (ier != 0) {
            std::cerr << "ERROR: could not find variable \"" << vname << "\" in file \"" << srcFileName << "\"\n";
            return 9;
        }

        // get the attributes of the variable from the netcdf file
        NcAttributes_t* attrs = NULL;
        ier = mnt_ncattributes_new(&attrs);
        // read the attributes
        ier = mnt_ncattributes_read(&attrs, srcNcid, srcVarid);
        if (ier != 0) {
            std::cerr << "ERROR: could not extract attributes for variable \"" 
                      << vname << "\" in file \"" << srcFileName << "\"\n";
            return 11;
        }

        // read the dimensions
        NcDimensions_t* srcVarDimsObj = NULL;
        ier = mnt_ncdimensions_new(&srcVarDimsObj);
        ier = mnt_ncdimensions_read(&srcVarDimsObj, srcNcid, srcVarid);
        int srcNdims;
        ier = mnt_ncdimensions_getNumDims(&srcVarDimsObj, &srcNdims);
        size_t srcDims[srcNdims];
        for (int i = 0; i < srcNdims; ++i) {
            ier = mnt_ncdimensions_get(&srcVarDimsObj, i, &srcDims[i]); 
        }      
        ier = mnt_ncdimensions_del(&srcVarDimsObj);

        // prepare to read the field 
        ier = mnt_ncfieldread_new(&reader, srcNcid, srcVarid);

        // get the number of edges and allocate src/dst data
        size_t numSrcEdges, numDstEdges;
        mnt_regridedges_getNumSrcEdges(&rg, &numSrcEdges);
        mnt_regridedges_getNumDstEdges(&rg, &numDstEdges);
        std::cout << "info: number of src edges: " << numSrcEdges << '\n';
        std::cout << "info: number of dst edges: " << numDstEdges << '\n';
        std::vector<double> srcEdgeData(numSrcEdges);
        std::vector<double> dstEdgeData(numDstEdges);

        size_t numSrcCells, numDstCells;
        mnt_regridedges_getNumSrcCells(&rg, &numSrcCells);
        mnt_regridedges_getNumDstCells(&rg, &numDstCells);

        if (dstEdgeDataFile.size() > 0) {
            // user provided a file name to store the regridded data

            std::pair<std::string, std::string> fm = split(dstEdgeDataFile, ':');
            // get the dst file name
            std::string dstFileName = fm.first;

            int n1 = dstFileName.size();
            int n2 = vname.size();
            const int append = 0; // new file
            ier = mnt_ncfieldwrite_new(&writer, dstFileName.c_str(), n1, vname.c_str(), n2, append);
            if (ier != 0) {
                std::cerr << "ERROR: create file " << dstFileName << " with field " 
                          << vname << " in append mode " << append << '\n';
                return 14;
            }

            ier = mnt_ncfieldwrite_setNumDims(&writer, srcNdims); // matches the number of source field dimensions
            if (ier != 0) {
                std::cerr << "ERROR: cannot set the number of dimensions for field " << vname << " in file " << dstFileName << '\n';
                ier = mnt_ncfieldwrite_del(&writer);
                return 15;
            }

            // add the field's axes. Assume the dst field dimensions are the same as the src field except for the last
            // num edges dimension


            // add num_edges axis. WE SHOULD GET THIS FROM THE DEST FILE?
            std::string axname = "num_edges";
            int n3 = axname.size();
            ier = mnt_ncfieldwrite_setDim(&writer, srcNdims - 1, axname.c_str(), n3, numDstEdges);
            if (ier != 0) {
                std::cerr << "ERROR: setting dimension 0 (" << axname << ") to " << numDstEdges
                          << " for field " << vname << " in file " << dstFileName << '\n';
                ier = mnt_ncfieldwrite_del(&writer);
                return 16;
            }

            // add the remaining axes, ASSUME THE ADDITIONAL DST AXES TO MATCH THE SRC AXES
            for (int i = 0; i < srcNdims - 1; ++i) {
                axname = "n_" + toString(srcDims[i]);
                ier = mnt_ncfieldwrite_setDim(&writer, i, axname.c_str(), axname.size(), srcDims[i]);
            }

            // add the attributes
            ier = mnt_ncattributes_write(&attrs, writer->ncid, writer->varid);
            if (ier != 0) {
                std::cerr << "ERROR: writing attributes for field " << vname << " in file " << dstFileName << '\n';
                ier = mnt_ncfieldwrite_del(&writer);
                return 17;
            }
        }

        std::vector<double> loop_integrals;
        std::vector<double> dstCellByCellData;

        // allocate
        loop_integrals.resize(numDstCells);
        dstCellByCellData.resize(numDstCells * NUM_EDGES_PER_CELL);

        // attach field to grid so we can save the data to file
        mnt_grid_attach(&rg->dstGridObj, vname.c_str(), NUM_EDGES_PER_CELL, &dstCellByCellData[0]);

        std::string loop_integral_varname = std::string("loop_integrals_of_") + vname;
        mnt_grid_attach(&rg->dstGridObj, loop_integral_varname.c_str(), 1, &loop_integrals[0]);

        //
        // iterate over the axes other than edge. ASSUMES num edges is the last dimension!!!
        //

        // leading indices into the src/dst array for each slice
        std::vector<size_t> srcIndices(srcNdims, 0);
        std::vector<size_t> dstIndices(srcNdims, 0);

        // number of data values to read, regrid and write
        std::vector<size_t> srcCounts(srcNdims, 1);
        srcCounts[srcNdims - 1] = numSrcEdges;
        std::vector<size_t> dstCounts(srcNdims, 1);
        dstCounts[srcNdims - 1] = numDstEdges;

        MultiArrayIter_t* mai = NULL;
        // iterate over the non-edge indices only, hence srcNdims - 1
        ier = mnt_multiarrayiter_new(&mai, srcNdims - 1, srcDims);

        // total number of elevation * time values
        size_t numIters;
        ier = mnt_multiarrayiter_getNumIters(&mai, &numIters);
        for (size_t iter = 0; iter < numIters; ++iter) {

            // same indices in this version, in later versions srcIndices and dstIndices 
            // might be different
            ier = mnt_multiarrayiter_getIndices(&mai, &srcIndices[0]);
            ier = mnt_multiarrayiter_getIndices(&mai, &dstIndices[0]);

            // read a slice of the data from file
            std::cout << "info: reading slice " << iter << " of field " << vname << " from file \"" << srcFileName << "\"\n";
            ier = mnt_ncfieldread_dataSlice(&reader, &srcIndices[0], &srcCounts[0], &srcEdgeData[0]);
            if (ier != 0) {
                std::cerr << "ERROR: could not read variable \"" << vname << "\" from file \"" << srcFileName << "\"\n";
                return 12;
            }

            // apply the weights to the src field
            ier = mnt_regridedges_apply(&rg, &srcEdgeData[0], &dstEdgeData[0]);
            if (ier != 0) {
                std::cerr << "ERROR: failed to apply weights to dst field \"" << vname << "\"\n";
                return 13;
            }

            // compute loop integrals for each cell
            double avgAbsLoop, minAbsLoop, maxAbsLoop;
            computeLoopIntegrals(rg->dstGridObj, dstEdgeData, &avgAbsLoop, &minAbsLoop, &maxAbsLoop, loop_integrals);
            std::cout << "Min/avg/max cell loop integrals: " << minAbsLoop << "/" << avgAbsLoop << "/" << maxAbsLoop << '\n';

            if (vtkOutputFile.size() > 0) {

                // new file name for each elevation, time, etc
                std::string vtkFilename = std::regex_replace(vtkOutputFile, 
                                          std::regex(".vtk"), "_" + toString(iter) + ".vtk");

                // compute the cell by cell data
                size_t dstEdgeId;
                int dstEdgeSign;
                for (size_t dstCellId = 0; dstCellId < numDstCells; ++dstCellId) {
                    for (int ie = 0; ie < 4; ++ie) {
                        ier = mnt_grid_getEdgeId(&rg->dstGridObj, dstCellId, ie, &dstEdgeId, &dstEdgeSign);
                        size_t k = dstCellId*NUM_EDGES_PER_CELL + ie;
                        dstCellByCellData[k] = dstEdgeData[dstEdgeId] * dstEdgeSign;
                    }
                }

                std::cout << "info: writing \"" << vname << "\" to " << vtkFilename << '\n';
                mnt_grid_dump(&rg->dstGridObj, vtkFilename.c_str());
            }

            if (dstEdgeDataFile.size() > 0) {

                // write the slice of data to a netcdf file
                ier = mnt_ncfieldwrite_dataSlice(&writer, &dstIndices[0], &dstCounts[0], &dstEdgeData[0]);
                if (ier != 0) {
                    std::cerr << "ERROR: writing slice " << iter << " of data for field " 
                              << vname << " in file " << dstEdgeDataFile << '\n';
                    ier = mnt_ncfieldwrite_del(&writer);
                    return 18;
                }

            }

            // increment the iterator (next time, elevation....)
            ier = mnt_multiarrayiter_next(&mai);

        }

        ier = mnt_multiarrayiter_del(&mai);

        if (dstEdgeDataFile.size() > 0) {
            ier = mnt_ncfieldwrite_del(&writer);
        }

        // must destroy before closing the file
        ier = mnt_ncfieldread_del(&reader);
        ier = mnt_ncattributes_del(&attrs);

        // done with reading the attributes
        ier = nc_close(srcNcid);

    } // has variable 
    else {
        std::cout << "info: no variable name was provided, thus only computing weights\n";
    }

    // clean up
    mnt_regridedges_del(&rg);

    return 0;
}
