using System.Diagnostics;
using System.Numerics;
using System.Text.Json;
using System.Threading;

internal class Program
{
    static void Main()
    {
        Orbit[] orbits;
        try
        {
            orbits = GetOrbitsFromJSON();
        }
        catch (Exception e)
        {
            Console.WriteLine(e);
            return;
        }
        foreach(Orbit o in orbits)
        {
            Console.WriteLine(o.ToStringDetail());
        }

        double nearestApproach = GetNearestApproach(orbits[0], orbits[1]);
        Console.WriteLine($"Nearest approach of {nearestApproach} calculated");
    }

    static Orbit[] GetOrbitsFromJSON()
    {
        string filePath = Path.Combine(AppDomain.CurrentDomain.BaseDirectory, "Orbits.json");
        string jsonContent = File.ReadAllText(filePath);

        return JsonSerializer.Deserialize<OrbitArray>(jsonContent).orbits;
    }
    static double GetNearestApproach(Orbit a, Orbit b, double targetRange = 1e-3, (double, double, double)[]? aSamples = null, (double, double, double)[]? bSamples = null)
    {
        if (aSamples == null)
            aSamples = a.SampleAngleBasedPoints(360);
        if (bSamples == null)
            bSamples = b.SampleAngleBasedPoints(360);

        double previousSquaredDistance = GetShortestDistanceBetweenPointAndOrbit(aSamples[aSamples.Length - 1], b, 1e-5, bSamples);
        double secondPreviousSquaredDistance = GetShortestDistanceBetweenPointAndOrbit(aSamples[aSamples.Length - 2], b, 1e-5, bSamples);
        List<double> localMinima = new();
        for(int i = 0; i < aSamples.Length + 2; i++)
        {
            double squaredDistance = GetShortestDistanceBetweenPointAndOrbit(aSamples[i % aSamples.Length], b, 1e-5, bSamples);
            // Console.WriteLine($"i {i}: {squaredDistance}"); // DEBUG
            if (secondPreviousSquaredDistance > previousSquaredDistance && previousSquaredDistance < squaredDistance)
            {
                // Console.WriteLine("Local minima found!"); // DEBUG
                localMinima.Add((i - 2d) / aSamples.Length);
            }
            secondPreviousSquaredDistance = previousSquaredDistance;
            previousSquaredDistance = squaredDistance;
        }
        
        for(int i = 0; i < localMinima.Count; i++)
        {
            // Console.WriteLine($"Approaching local minima {i}"); // DEBUG
            localMinima[i] = ApproachLocalMinima(a, b, localMinima[i], localMinima[i] + 2d / aSamples.Length, targetRange, bSamples);
            // Console.WriteLine($"Approach complete: {localMinima[i]}"); // DEBUG
        }
        return localMinima.Min();
    }
    static double ApproachLocalMinima(Orbit a, Orbit b, double lowerBound, double upperBound, double targetRange = 1e-3, (double, double, double)[]? bSamples = null)
    {
        while(upperBound - lowerBound > targetRange)
        {
            (double, double, double)[] samples = a.SampleAngleBasedPoints(10, lowerBound, upperBound);

            double prevSquaredDistance = GetShortestDistanceBetweenPointAndOrbit(samples[samples.Length - 1], b, 1e-5, bSamples);
            double secondPrevSquaredDistance = GetShortestDistanceBetweenPointAndOrbit(samples[samples.Length - 1], b, 1e-5, bSamples);
            for (int i = 0; i <= samples.Length; i++)
            {
                double squaredDistance = GetShortestDistanceBetweenPointAndOrbit(samples[i % samples.Length], b, 1e-5, bSamples);
                if (secondPrevSquaredDistance > prevSquaredDistance && prevSquaredDistance < squaredDistance)
                {
                    (lowerBound, upperBound) = GetNextBounds(lowerBound, upperBound, (i - 2d) / samples.Length, (double)i / samples.Length);
                    break;
                }
                secondPrevSquaredDistance = prevSquaredDistance;
                prevSquaredDistance = squaredDistance;
            }
        }
        return GetShortestDistanceBetweenPointAndOrbit(a.GetPositionAtAngle((upperBound + lowerBound) / 2), b, 1e-5, bSamples);
    }
    static double GetShortestDistanceBetweenPointAndOrbit((double, double, double) point, Orbit orbit, double targetRange = 1e-5, (double, double, double)[]? samples = null)
    {
        if (samples == null)
            samples = orbit.SampleAngleBasedPoints(360);

        double minimumLowerBound = 0;
        double minimumUpperBound = 1;

        // Refine bounds
        do
        {
            double prevSquaredDistance = GetSquaredDistance(point, samples[samples.Length - 1]);
            double secondPrevSquaredDistance = GetSquaredDistance(point, samples[samples.Length - 2]);
            for (int i = 0; i <= samples.Length; i++)
            {
                double squaredDistance = GetSquaredDistance(point, samples[i % samples.Length]);
                if (secondPrevSquaredDistance >= prevSquaredDistance && prevSquaredDistance < squaredDistance)
                {
                    (minimumLowerBound, minimumUpperBound) = GetNextBounds(minimumLowerBound, minimumUpperBound, (i - 2d) / samples.Length, (double) i / samples.Length);
                    break;
                }
                if(squaredDistance == prevSquaredDistance && prevSquaredDistance == secondPrevSquaredDistance) // whoa, we hit a floating point accuracy limit! upper bound and lower bound are practically the same in this case
                    return squaredDistance;
                secondPrevSquaredDistance = prevSquaredDistance;
                prevSquaredDistance = squaredDistance;
            }

            samples = orbit.SampleAngleBasedPoints(10, minimumLowerBound, minimumUpperBound);
        } while (minimumUpperBound - minimumLowerBound > targetRange);
        return GetDistance(point, orbit.GetPositionAtAngle((minimumLowerBound + minimumUpperBound) / 2));
    }
    static (double, double) GetNextBounds(double curLowerBound, double curUpperBound, double lowerPercentage, double upperPercentage)
    {
        return (curLowerBound + lowerPercentage * (curUpperBound - curLowerBound), curLowerBound + upperPercentage * (curUpperBound - curLowerBound));
    }
    static double GetDistance((double, double, double) p1, (double, double, double) p2)
    {
        return Math.Sqrt(GetSquaredDistance(p1, p2));
    }
    static double GetSquaredDistance((double, double, double) p1, (double, double, double) p2)
    {
        (double, double, double) diff = (p1.Item1 - p2.Item1, p1.Item2 - p2.Item2, p1.Item3 - p2.Item3);
        return diff.Item1 * diff.Item1 + diff.Item2 * diff.Item2 + diff.Item3 * diff.Item3;
    }
    private struct OrbitArray
    {
        public Orbit[] orbits { get; set; }
    }
    private struct ApproachLocalMinimumArgs
    {
        public Orbit o;
        public Orbit p;
        public double oAngle;
        public double oRange;
    }
}