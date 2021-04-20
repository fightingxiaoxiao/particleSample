#include "particleSampleContainer.H"

namespace Foam
{

    void particleSampleContainer::classifyDiameterAlongHeight(scalar startHeight,
                                                              scalar deltaH)
    {
        for (auto &p : particleStorage)
        {
            label hIndex = static_cast<label>((p.second.position[2] - startHeight) / deltaH);
            if (hIndex < 0)
                continue;
            diameterList[hIndex].push_back(p.second.d);
        }
    }

    void particleSampleContainer::classifyVelocityAlongHeight(scalar startHeight,
                                                              scalar deltaH)
    {
        for (auto &p : particleStorage)
        {
            label hIndex = static_cast<label>((p.second.position[2] - startHeight) / deltaH);
            if (hIndex < 0)
                continue;
            velocityList[hIndex].push_back(mag(p.second.U));
        }
    }

    void particleSampleContainer::classifyFlowRateAlongHeight(scalar startHeight,
                                                              scalar deltaH,
                                                              scalar directionIndex,
                                                              scalar samplePosition,
                                                              scalar limitMoveDistanceInOneSample)
    {
        for (auto &p : particleStorage)
        {
            label hIndex = static_cast<label>((p.second.position[2] - startHeight) / deltaH);
            if (hIndex < 0)
                continue;
            if (massRateList.find(hIndex) == massRateList.end())
            {
                massRateList[hIndex] = 0.;
            }

            auto pOld = _particleStorage[p.first];
            auto _position = pOld.position;
            auto position = p.second.position;
            if (mag(position - _position) > limitMoveDistanceInOneSample)
            {
                continue;
            }

            if (_position[directionIndex] < samplePosition && position[directionIndex] >= samplePosition)
            {
                massRateList[hIndex] += p.second.nParticle * p.second.mass();
            }
            else if (_position[directionIndex] > samplePosition && position[directionIndex] <= samplePosition)
            {
                massRateList[hIndex] -= p.second.nParticle * p.second.mass();
            }
        }
    }

    word particleSampleContainer::writeDiameterInfo(scalar startHeight,
                                                    scalar deltaH)
    {
        std::stringstream diameterInfo("");
        diameterInfo << "height, diameter0, diameter1..." << std::endl;
        for (auto &d : diameterList)
        {
            diameterInfo << startHeight + d.first * deltaH;

            for (auto &data : d.second)
            {
                diameterInfo << "," << data;
            }

            diameterInfo << std::endl;
        }
        return diameterInfo.str();
    }

    word particleSampleContainer::writeVelocityInfo(scalar startHeight,
                                                    scalar deltaH)
    {
        std::stringstream velocityInfo("");
        velocityInfo << "height, velmag0, velmag1..." << std::endl;
        for (auto &v : velocityList)
        {
            velocityInfo << startHeight + v.first * deltaH;

            for (auto &data : v.second)
            {
                velocityInfo << "," << data;
            }

            velocityInfo << std::endl;
        }
        return velocityInfo.str();
    }

}
