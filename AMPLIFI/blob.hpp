//
//  blob.hpp
//  AMPLIFI
//
//  Created by Ningyu Liu on 5/1/25.
//

#ifndef blob_hpp
#define blob_hpp

#include <array>
#include <cmath>
#include <cassert>
#include <vector>
#include <memory>

#include "Vector.H"
#include "RefCountedPtr.H"
#include "SPACE.H"

#include "NamespaceHeader.H"

class Blob {
public:
    virtual ~Blob() = default;
    virtual double value(const std::vector<double>& x) const = 0;
};

/* tanh profile: value(x) = amplitude * 0.5 * {1 - tanh[sharpness * (distance - radius)]}, where
the distance is to the respective focus. Rewrite it as
value(x) = amplitude * 0.5 * {1 - tanh[sharpness/radius * (distance/radius - 1.0)]} or
value(x) = amplitude * 0.5 * {1 - tanh[sharpnessNorm * (distance/radius - 1.0)]}, i.e., sharpnessNorm is
the sharpness normalized to radius. sharpnessNorm is what should be given as an input. */

class SpheroidalBlob : public Blob {
public:
    SpheroidalBlob(
        const std::vector<double>& center,
        const std::vector<double>& axisDirection,
        double axisLength,
        double radius,
        double sharpness,
        double amplitude = 1.0
    ) : m_center(center), m_axisLength(axisLength), m_radius(radius),
        m_sharpness(sharpness), m_amplitude(amplitude)
    {
        assert(center.size() == SpaceDim);
        assert(axisDirection.size() == SpaceDim);

        double norm = 0.0;
        m_axis.resize(SpaceDim);
        for (int d = 0; d < SpaceDim; ++d) norm += axisDirection[d] * axisDirection[d];
        norm = std::sqrt(norm);
        assert(norm > 0.0);

        for (int d = 0; d < SpaceDim; ++d)
            m_axis[d] = axisDirection[d] / norm;
    }

    double value(const std::vector<double>& x) const override {
        assert(x.size() == SpaceDim);
        std::vector<double> delta(SpaceDim);
        for (int d = 0; d < SpaceDim; ++d)
            delta[d] = x[d] - m_center[d];

        double l = 0.0;
        for (int d = 0; d < SpaceDim; ++d)
            l += delta[d] * m_axis[d];

      double s2 = 0.0;
      for (int d = 0; d < SpaceDim; ++d) {
          double diff = delta[d] - l * m_axis[d];
          s2 += diff * diff;
      }

        double normalizedDist = sqrt((l * l) / (m_axisLength * m_axisLength) + s2 / (m_radius * m_radius));
      
        return m_amplitude * 0.5 * (1.0 - std::tanh(m_sharpness * (normalizedDist - 1)));
    }

private:
    std::vector<double> m_center;
    std::vector<double> m_axis;
    double m_axisLength;
    double m_radius;
    double m_sharpness;
    double m_amplitude;
};

class CylinderBlob : public Blob {
public:
    CylinderBlob(
        const std::vector<double>& center,
        const std::vector<double>& axisDirection,
        double axisLength,
        double radius,
        double sharpness,
        double amplitude = 1.0
    ) : m_center(center), m_axisLength(axisLength), m_radius(radius),
        m_sharpness(sharpness), m_amplitude(amplitude)
    {
        assert(center.size() == SpaceDim);
        assert(axisDirection.size() == SpaceDim);

        double norm = 0.0;
        m_axis.resize(SpaceDim);
        for (int d = 0; d < SpaceDim; ++d) norm += axisDirection[d] * axisDirection[d];
        norm = std::sqrt(norm);
        assert(norm > 0.0);

        for (int d = 0; d < SpaceDim; ++d)
            m_axis[d] = axisDirection[d] / norm;
    }

    double value(const std::vector<double>& x) const override {
        assert(x.size() == SpaceDim);
        std::vector<double> delta(SpaceDim);
        for (int d = 0; d < SpaceDim; ++d)
            delta[d] = x[d] - m_center[d];

        double l = 0.0;
        for (int d = 0; d < SpaceDim; ++d)
            l += delta[d] * m_axis[d];

      double s2 = 0.0;
      for (int d = 0; d < SpaceDim; ++d) {
          double diff = delta[d] - l * m_axis[d];
          s2 += diff * diff;
      }

        double dist;
        double halfLength = 0.5 * m_axisLength;
        if (std::fabs(l) <= halfLength) {
            dist = std::sqrt(s2);
        } else {
            double dl = std::fabs(l) - halfLength;
            dist = std::sqrt(s2 + dl * dl);
        }
      double normalizedDist = dist/m_radius;
    
      return m_amplitude * 0.5 * (1.0 - std::tanh(m_sharpness * (normalizedDist - 1)));
    }

private:
    std::vector<double> m_center;
    std::vector<double> m_axis;
    double m_axisLength;
    double m_radius;
    double m_sharpness;
    double m_amplitude;
};

class MultiBlob : public Blob {
public:
    MultiBlob() = default;

    void addBlob(const RefCountedPtr<Blob>& blob)
    {
        m_blobs.push_back(blob);
    }

    int numBlobs() const
    {
        return m_blobs.size();
    }

    double value(const std::vector<double>& x) const override
    {
        CH_assert(x.size() == SpaceDim);

        double sum = 0.0;
        for (int i = 0; i < m_blobs.size(); ++i) {
            sum += m_blobs[i]->value(x);
        }
        return sum;
    }

private:
    Vector<RefCountedPtr<Blob> > m_blobs;
};

extern void parseBlobsFromParmParse(MultiBlob& multiBlob);
extern int testBlob();

#include "NamespaceFooter.H"
#endif /* blob_hpp */
