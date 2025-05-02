//
//  blob.cpp
//  AMPLIFI
//
//  Created by Ningyu Liu on 5/1/25.
//

#include "ParmParse.H"

#include "blob.hpp"
#include "globalVariables.h"

#include "NamespaceHeader.H"  // This is often required for Chombo code

void parseBlobsFromParmParse(MultiBlob& multiBlob)
{
    ParmParse pp("iniPlaCloud");
    int blobNum;
    pp.get("number", blobNum);

    // Initialize vectors to store blob parameters
    Vector<Real> centers(blobNum * SpaceDim, 0.0);
    Vector<Real> mags(blobNum, 1.0);
    Vector<Real> sharpness(blobNum, 20.0); // Optional sharpness parameter
    Vector<std::string> types(blobNum, "spheroid"); // Default to spheroid
    Vector<Real> axisDirs(blobNum * SpaceDim, 0.0); // Axis direction for each blob
    Vector<Real> axisLengths(blobNum, 0.0); // Axis length for spheroids and cylinders
    Vector<Real> radius(blobNum, 0.0); // Single radius for spheroid, radial radius for cylinder
    
    pp.getarr("center", centers, 0, blobNum * SpaceDim);
    pp.getarr("radius", radius, 0, blobNum);  // One radius per blob (either radial or spherical)
    pp.queryarr("mag", mags, 0, blobNum);
    pp.queryarr("sharpness", sharpness, 0, blobNum);
    pp.queryarr("type", types, 0, blobNum); // spheroid or cylinder
    pp.queryarr("axisDirection", axisDirs, 0, blobNum * SpaceDim); // Axis direction for each blob
    pp.queryarr("axisLength", axisLengths, 0, blobNum); // Axis length for both spheroids and cylinders
    
  // Apply normalization
    for (int i = 0; i < blobNum; i++) {
        for (int d = 0; d < SpaceDim; d++) {
            centers[d + i * SpaceDim] /= normalization::scalingFactor / normalization::lBar;
        }
        radius[i] /= normalization::scalingFactor / normalization::lBar;
        mags[i] *= normalization::scalingFactor * normalization::scalingFactor / normalization::nBar;
        axisLengths[i] /= normalization::scalingFactor / normalization::lBar;
        sharpness[i] *= normalization::scalingFactor / normalization::lBar;  // If sharpness is scaled by length
    }
  
    for (int i = 0; i < blobNum; ++i)
    {
        std::vector<Real> center(SpaceDim), axisDir(SpaceDim);
        for (int d = 0; d < SpaceDim; ++d)
        {
            center[d] = centers[i * SpaceDim + d];
            axisDir[d] = axisDirs[i * SpaceDim + d];
        }

        // Normalize the axis direction
        Real norm = 0.0;
        for (int d = 0; d < SpaceDim; ++d)
        {
            norm += axisDir[d] * axisDir[d];
        }
        norm = std::sqrt(norm);
        assert(norm > 0.0);  // Ensure the axis direction is not zero
        for (int d = 0; d < SpaceDim; ++d)
        {
            axisDir[d] /= norm;  // Normalize the axis direction
        }

        // Construct the correct blob type based on the input
        RefCountedPtr<Blob> blob;

        if (types[i] == "spheroid")
        {
            if (SpaceDim == 2)
            {
                // 2D spheroid (circle)
                Real radiusValue = radius[i]; // Only one radius value for 2D
                blob = RefCountedPtr<Blob>(
                    new SpheroidalBlob(center, axisDir, axisLengths[i], radiusValue,
                                       sharpness[i], mags[i]));
            }
            else
            {
                // 3D spheroid (ellipsoid)
                Real radiusValue = radius[i]; // Use the first radius value for 3D
                blob = RefCountedPtr<Blob>(
                    new SpheroidalBlob(center, axisDir, axisLengths[i], radiusValue,
                                       sharpness[i], mags[i]));
            }
        }
        else if (types[i] == "cylinder")
        {
            if (SpaceDim == 2)
            {
                MayDay::Error("Cylinder blobs are not supported in 2D.");
            }
            else
            {
                Real radialRadius = radius[i];  // Radial radius for the cylinder

                blob = RefCountedPtr<Blob>(
                    new CylinderBlob(center, axisDir, axisLengths[i], radialRadius,
                                     sharpness[i], mags[i]));
            }
        }
        else
        {
            MayDay::Error("Unknown blob type in input.");
        }

        multiBlob.addBlob(blob);
    }
}

#include "NamespaceFooter.H"
