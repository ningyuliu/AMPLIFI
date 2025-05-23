//
//  blob.cpp
//  AMPLIFI
//
//  Created by Ningyu Liu on 5/1/25.
//

#include <string>

#include "SPMD.H"
#include "ParmParse.H"

#include "blob.hpp"
#include "globalVariables.h"

#include "NamespaceHeader.H"  // This is often required for Chombo code

void generateRandomBlobs(MultiBlob& multiBlob) {
  ParmParse pp("iniPlaCloud");
  int blobNum;
  pp.get("randBlobNum", blobNum);

  Vector<Real> center(SpaceDim * blobNum, 0.0);
  Vector<Real> axis(SpaceDim * blobNum, 0.0);
  Vector<Real> radius(blobNum, 0.0);
  Vector<Real> axisLength(blobNum, 0.0);
  Vector<Real> mag(blobNum, 0.0);
  Vector<Real> sharpness(blobNum, 0.0);

  Vector<Real> distCenter(SpaceDim, 0.0);
  Vector<Real> distBoxLength(SpaceDim, 0.0);
  Vector<Real> radLim(2), axisLengthLim(2), magLim(2), sharpLim(2);
  pp.getarr("distCenter", distCenter, 0, SpaceDim);
  pp.getarr("distBoxLength", distBoxLength, 0, SpaceDim);
  pp.getarr("radLim", radLim, 0, 2);
  pp.getarr("axisLengthLim", axisLengthLim, 0, 2);
  pp.getarr("magLim", magLim, 0, 2);
  pp.getarr("sharpLim", sharpLim, 0, 2);

  for (int d = 0; d < SpaceDim; d++) {
    distCenter[d] /= normalization::scalingFactor * normalization::lBar;
    distBoxLength[d] /= normalization::scalingFactor * normalization::lBar;
  }
  for (int i = 0; i < 2; i++) {
    radLim[i] /= normalization::scalingFactor * normalization::lBar;
    axisLengthLim[i] /= normalization::scalingFactor * normalization::lBar;
    magLim[i] *= normalization::scalingFactor * normalization::scalingFactor / normalization::nBar;
  }

  const int masterProc = uniqueProc(SerialTask::compute);

  // String buffers for broadcasting (only filled by master)
  std::string centerStr, axisStr, radiusStr, axisLengthStr, magStr, sharpStr;

  if (procID() == masterProc) {
    for (int i = 0; i < blobNum; i++) {
      for (int d = 0; d < SpaceDim; d++) {
        center[i * SpaceDim + d] = ((1.0 * rand() / RAND_MAX) - 0.5) * distBoxLength[d] + distCenter[d];
        centerStr += std::to_string(center[i * SpaceDim + d]) + " ";
      }

      radius[i] = (rand() / static_cast<double>(RAND_MAX)) * (radLim[1] - radLim[0]) + radLim[0];
      radiusStr += std::to_string(radius[i]) + " ";

      axisLength[i] = (rand() / static_cast<double>(RAND_MAX)) * (axisLengthLim[1] - axisLengthLim[0]) + axisLengthLim[0];
      axisLengthStr += std::to_string(axisLength[i]) + " ";

      mag[i] = (rand() / static_cast<double>(RAND_MAX)) * (magLim[1] - magLim[0]) + magLim[0];
      magStr += std::to_string(mag[i]) + " ";

      sharpness[i] = (rand() / static_cast<double>(RAND_MAX)) * (sharpLim[1] - sharpLim[0]) + sharpLim[0];
      sharpStr += std::to_string(sharpness[i]) + " ";

      // Axis
      double norm = 0.0;
      for (int d = 0; d < SpaceDim; d++) {
        axis[i * SpaceDim + d] = ((1.0 * rand() / RAND_MAX) - 0.5);
        norm += axis[i * SpaceDim + d] * axis[i * SpaceDim + d];
      }
      norm = std::sqrt(norm);
      for (int d = 0; d < SpaceDim; d++) {
        axis[i * SpaceDim + d] /= norm;
        axisStr += std::to_string(axis[i * SpaceDim + d]) + " ";
      }
    }
  }

  // Broadcast the data
  broadcast(center, masterProc);
  broadcast(axis, masterProc);
  broadcast(radius, masterProc);
  broadcast(axisLength, masterProc);
  broadcast(mag, masterProc);
  broadcast(sharpness, masterProc);

  // Construct blobs
  for (int i = 0; i < blobNum; i++) {
    std::vector<double> blobCenter(SpaceDim), blobAxis(SpaceDim);
    for (int d = 0; d < SpaceDim; d++) {
      blobCenter[d] = center[i * SpaceDim + d];
      blobAxis[d] = axis[i * SpaceDim + d];
    }
    RefCountedPtr<Blob> blob(new SpheroidBlob(blobCenter, blobAxis, axisLength[i], radius[i], sharpness[i], mag[i]));
    multiBlob.addBlob(blob);
  }
}

void parseBlobsFromParmParse(MultiBlob& multiBlob) {
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
      centers[d + i * SpaceDim] = centers[d + i * SpaceDim] / normalization::scalingFactor / normalization::lBar;
    }
    radius[i] = radius[i] / normalization::scalingFactor / normalization::lBar;
    mags[i] = mags[i] * normalization::scalingFactor * normalization::scalingFactor / normalization::nBar;
    axisLengths[i] = axisLengths[i] / normalization::scalingFactor / normalization::lBar;
  }
  
  for (int i = 0; i < blobNum; ++i) {
    std::vector<Real> center(SpaceDim), axisDir(SpaceDim);
    for (int d = 0; d < SpaceDim; ++d) {
      center[d] = centers[i * SpaceDim + d];
      axisDir[d] = axisDirs[i * SpaceDim + d];
    }
    
    // Normalize the axis direction
    Real norm = 0.0;
    for (int d = 0; d < SpaceDim; ++d)
      norm += axisDir[d] * axisDir[d];
  
    norm = std::sqrt(norm);
    assert(norm > 0.0);  // Ensure the axis direction is not zero
    for (int d = 0; d < SpaceDim; ++d)
      axisDir[d] /= norm;  // Normalize the axis direction
    
    // Construct the correct blob type based on the input
    RefCountedPtr<Blob> blob;
    
    if (types[i] == "spheroid")
      if (SpaceDim == 2) {
        // 2D spheroid (circle)
        Real radiusValue = radius[i]; // Only one radius value for 2D
        blob = RefCountedPtr<Blob>(
                                   new SpheroidBlob(center, axisDir, axisLengths[i], radiusValue,
                                                      sharpness[i], mags[i]));
      }
      else {
        // 3D spheroid (ellipsoid)
        Real radiusValue = radius[i]; // Use the first radius value for 3D
        blob = RefCountedPtr<Blob>(
                                   new SpheroidBlob(center, axisDir, axisLengths[i], radiusValue,
                                                      sharpness[i], mags[i]));
      }
    else if (types[i] == "cylinder")
      if (SpaceDim == 2)
        MayDay::Error("Cylinder blobs are not supported in 2D.");
      else {
        Real radialRadius = radius[i];  // Radial radius for the cylinder
        
        blob = RefCountedPtr<Blob>(
                                   new CylinderBlob(center, axisDir, axisLengths[i], radialRadius,
                                                    sharpness[i], mags[i]));
      }
    else
      MayDay::Error("Unknown blob type in input.");
    
    multiBlob.addBlob(blob);
  }
}

void outputBlobs(const MultiBlob& multiBlob) {
  for (int i = 0; i < multiBlob.numBlobs(); ++i) {
    const RefCountedPtr<Blob>& blobPtr = multiBlob.m_blobs[i];

    if (const SpheroidBlob* sphBlob = dynamic_cast<const SpheroidBlob*>(blobPtr.operator->())) {
      std::cout << "Blob " << i + 1 << ": Spheroid\n";
      std::cout << "  Center: ";
      for (double c : sphBlob->m_center) std::cout << c << " ";
      std::cout << "\n  Axis: ";
      for (double a : sphBlob->m_axis) std::cout << a << " ";
      std::cout << "\n  Axis Length: " << sphBlob->m_axisLength
                << "\n  Radius: " << sphBlob->m_radius
                << "\n  Sharpness: " << sphBlob->m_sharpness
                << "\n  Amplitude: " << sphBlob->m_amplitude << "\n";

    } else if (const CylinderBlob* cylBlob = dynamic_cast<const CylinderBlob*>(blobPtr.operator->())) {
      std::cout << "Blob " << i + 1 << ": Cylinder\n";
      std::cout << "  Center: ";
      for (double c : cylBlob->m_center) std::cout << c << " ";
      std::cout << "\n  Axis: ";
      for (double a : cylBlob->m_axis) std::cout << a << " ";
      std::cout << "\n  Axis Length: " << cylBlob->m_axisLength
                << "\n  Radius: " << cylBlob->m_radius
                << "\n  Amplitude: " << cylBlob->m_amplitude << "\n";

    } else {
      std::cout << "Blob " << i + 1 << ": Unknown blob type\n";
    }
  }
}



// Simple function to print blob value at a point
void printValueAtPoint(const Blob& blob, const std::vector<double>& x, const std::string& label) {
  std::cout << label << " at [";
  for (int i = 0; i < SpaceDim; ++i) {
    std::cout << x[i] << (i < SpaceDim - 1 ? ", " : "");
  }
  std::cout << "] = " << blob.value(x) << "\n";
}

int testBlob() {
  // Center at origin
  std::vector<double> center = {0.0, 0.0, 0.0};
  
  // Axis aligned along z-direction
  std::vector<double> axis = {0.0, 0.0, 1.0};
  
  // Parameters
  double axisLength = 2.0;
  double radius = 1.0;
  double sharpness = 2.0;
  double amplitude = 1.0;
  
  // Create a spheroidal blob
  RefCountedPtr<Blob> spheroid(new SpheroidBlob(center, axis, axisLength, radius, sharpness, amplitude));
  
  // Create a cylindrical blob
  RefCountedPtr<Blob> cylinder(new CylinderBlob(center, axis, axisLength, radius, sharpness, amplitude));
  
  // Print test values at some points
  std::vector<std::vector<double>> testPoints = {
    {0.0, 0.0, 0.0},   // center
    {0.5, 0.0, 0.0},   // radial
    {0.0, 0.0, 1.0},   // axial
    {1.0, 0.0, 1.0},   // edge
    {2.0, 0.0, 0.0}    // outside
  };
  
  //    std::cout << "SpheroidBlob values:\n";
  //    for (const auto& pt : testPoints) {
  //        printValueAtPoint(*spheroid, pt, "  Value");
  //    }
  //
  //    std::cout << "\nCylinderBlob values:\n";
  //    for (const auto& pt : testPoints) {
  //        printValueAtPoint(*cylinder, pt, "  Value");
  //    }
  //
  //    // MultiBlob example
  //    MultiBlob multiBlob;
  //    multiBlob.addBlob(spheroid);
  //    multiBlob.addBlob(cylinder);
  //
  //    std::cout << "\nMultiBlob (Spheroid + Cylinder) values:\n";
  //    for (const auto& pt : testPoints) {
  //        printValueAtPoint(multiBlob, pt, "  Value");
  //    }
  
  center = {33070, 33070, 66140};
  axisLength = 6614;
  radius = 6614;
  sharpness = 0.001;
  SpheroidBlob sph(center, axis, axisLength, radius, sharpness, amplitude);
  testPoints = {
    //        {33070, 33070, 66140},   // center
    //        {33070+0.5*radius, 33070, 66140},   // radial
    //        {33070, 33070, 66140+0.5*axisLength+0.5*radius},   // axial
    //        {33070+0.5*radius, 33070, 66140+0.5*axisLength+0.5*radius},   // edge
    //        {33070+2*33070, 0.0, 0.0},    // outside
    {1111, 1111, 1111}
  };
  std::cout << "\nSpheroidBlob values:\n";
  for (const auto& pt : testPoints) {
    printValueAtPoint(sph, pt, "  Value");
  }
  
  return 0;
}

#include "NamespaceFooter.H"
