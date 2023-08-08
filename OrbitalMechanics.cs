// Shoutouts to Scott Anderson for these! I barely know anything about orbital mechanics LOL
// his repo on his implementation: https://github.com/ScottyRAnderson/Keplerian-Orbits
// his yt vid about this: https://www.youtube.com/watch?v=t89De819YMA (highly underrated btw, check him out)

// However his code is kinda incomplete and doesn't account for longitude of ascending node.
// I found an algorithm to account for it: https://downloads.rene-schwarz.com/download/M001-Keplerian_Orbit_Elements_to_Cartesian_State_Vectors.pdf

using System.Numerics;
using System.Text.Json.Serialization;

public struct Orbit
{
    public string Name { get; set; }
    public string ParentName { get; set; }
    public double Apoapsis
    {
        get { return Apoapsis; }
        set
        {
            Apoapsis = value;
            SemiMajorAxis = (Apoapsis + Periapsis) / 2;
            SemiMinorAxis = Math.Sqrt(Apoapsis * Periapsis);
            LinearEccentricity = SemiMajorAxis - Periapsis;
            Eccentricity = LinearEccentricity / SemiMajorAxis;
        }
    }
    public double Periapsis
    {
        get { return Periapsis; }
        set
        {
            Periapsis = value;
            SemiMajorAxis = (Apoapsis + Periapsis) / 2;
            SemiMinorAxis = Math.Sqrt(Apoapsis * Periapsis);
            LinearEccentricity = SemiMajorAxis - Periapsis;
            Eccentricity = LinearEccentricity / SemiMajorAxis;
        }
    }
    public double Inclination
    {
        get { return Inclination; }
        set
        {
            Inclination = value;
            transformationMatrix = GetTransformationMatrix();
        }
    }
    public double ArgumentOfPeriapsis
    {
        get { return ArgumentOfPeriapsis; }
        set
        {
            ArgumentOfPeriapsis = value;
            transformationMatrix = GetTransformationMatrix();
        }
    }
    public double LongitudeOfAscendingNode
    {
        get { return LongitudeOfAscendingNode; }
        set
        {
            LongitudeOfAscendingNode = value;
            transformationMatrix = GetTransformationMatrix();
        }
    }
    public double MeanAnomalyAtEpochRadians { get; set; }
    public double SemiMajorAxis { get; set; }
    public double SemiMinorAxis { get; set; }
    public double LinearEccentricity { get; set; }
    public double Eccentricity { get; set; }
    public double[,] transformationMatrix;
    public Orbit(double apoapsis = 1, double periapsis = 1, double inclination = 0, double argumentOfPeriapsis = 0, double longitudeOfAscendingNode = 0, double meanAnomalyAtEpochRadians = 0, string name = "Planet", string parentName = "Star")
    {
        Name = name;
        ParentName = parentName;
        if (periapsis > apoapsis)
            (apoapsis, periapsis) = (periapsis, apoapsis);
        Apoapsis = apoapsis;
        Periapsis = periapsis;
        Inclination = inclination;
        ArgumentOfPeriapsis = argumentOfPeriapsis;
        LongitudeOfAscendingNode = longitudeOfAscendingNode;
        MeanAnomalyAtEpochRadians = meanAnomalyAtEpochRadians;
        transformationMatrix = GetTransformationMatrix();
    }
    private double[,] GetTransformationMatrix()
    {
        double[,] matrix = new double[3,2];
        double sinInc   = Math.Sin(Inclination),              cosInc   = Math.Cos(Inclination);
        double sinArgPe = Math.Sin(ArgumentOfPeriapsis),      cosArgPe = Math.Cos(ArgumentOfPeriapsis);
        double sinLAN   = Math.Sin(LongitudeOfAscendingNode), cosLAN   = Math.Cos(LongitudeOfAscendingNode);

        // https://downloads.rene-schwarz.com/download/M001-Keplerian_Orbit_Elements_to_Cartesian_State_Vectors.pdf
        matrix[0, 0] = sinArgPe * cosLAN - sinArgPe * cosInc * sinLAN;
        matrix[0, 1] = -(sinArgPe * cosLAN + cosArgPe * cosInc * sinLAN);

        matrix[1, 0] = cosArgPe * sinLAN + sinArgPe * cosInc * cosLAN;
        matrix[1, 1] = cosArgPe * cosInc * cosLAN - sinArgPe * sinLAN;

        matrix[2, 0] = sinArgPe * sinInc;
        matrix[2, 1] = cosArgPe * sinInc;

        return matrix;
    }
    public double GetEccentricAnomaly(double meanAnomaly)
    {
        double targetAccuracy = 1e-9;
        int maxIterations = 1000;

        // Starting guess
        double eccentricAnomaly = Eccentricity > 0.8 ? Math.PI : Eccentricity;

        for (int i = 0; i < maxIterations; i++)
        {
            // NEWTON'S METHOD
            // x_n+1 = x_n + f(x_n)/f'(x_n)
            double nextValue = eccentricAnomaly - (OrbitalMechanics.KeplersEquation(meanAnomaly, eccentricAnomaly, Eccentricity) / OrbitalMechanics.KeplersEquation_dE(eccentricAnomaly, Eccentricity));
            double diff = Math.Abs(eccentricAnomaly - nextValue);
            eccentricAnomaly = nextValue;

            if (diff < targetAccuracy)
                break;
        }

        return eccentricAnomaly;
    }
    public (double, double, double) GetPositionAtTime(double t)
    {
        t = Modulo(t, 1);
        double meanAnomaly = t * Math.Tau + MeanAnomalyAtEpochRadians;
        double eccentricAnomaly = GetEccentricAnomaly(meanAnomaly);

        double x = SemiMajorAxis * (Math.Cos(eccentricAnomaly) - Eccentricity);
        double y = SemiMinorAxis * Math.Sin(eccentricAnomaly);

        return TiltPosition(x, y);
    }

    public (double, double, double) GetPositionAtAngle(double revolutions)
    {
        revolutions = Modulo(revolutions, 1);
        double trueAnomaly = revolutions * Math.Tau;

        double x = SemiMajorAxis * (Math.Cos(trueAnomaly) - Eccentricity);
        double y = SemiMinorAxis * Math.Sin(trueAnomaly);

        return TiltPosition(x, y);
    }
    public (double, double, double) TiltPosition(double x, double y)
    {
        double newX = x * transformationMatrix[0, 0] + y * transformationMatrix[0, 1];
        double newY = x * transformationMatrix[1, 0] + y * transformationMatrix[1, 1];
        double newZ = x * transformationMatrix[2, 0] + y * transformationMatrix[2, 1];
        return (newX, newY, newZ);
    }

    public (double, double, double)[] SampleAngleBasedPoints(int samples)
    {
        (double, double, double)[] sampleArray = new (double, double, double)[samples];
        double samplesAsDouble = samples;
        for(int i = 0; i < samples; i++)
            sampleArray[i] = GetPositionAtAngle(i / samplesAsDouble);
        return sampleArray;
    }
    public (double, double, double)[] SampleAngleBasedPoints(int samples, double start, double end)
    {
        (double, double, double)[] sampleArray = new (double, double, double)[samples];
        double samplesAsDouble = samples;
        for (int i = 0; i < samples; i++)
            sampleArray[i] = GetPositionAtAngle(start + (i / (samplesAsDouble * (end - start))));
        return sampleArray;
    }

    public override string ToString()
    {
        return $"{ParentName}>{Name}: {Apoapsis}x{Periapsis}";
    }
    public string ToStringDetail()
    {
        return $"{ParentName}>{Name}: {Apoapsis}x{Periapsis}, {Inclination}°i/{ArgumentOfPeriapsis}°ω/{LongitudeOfAscendingNode}°Ω";
    }

    // A modulo implementation that supports negatives.
    private double Modulo(double a, double b) { return (a % b + b) % b; }
}


// Static version
public static class OrbitalMechanics
{
    /// <param name="M">Mean Anomaly</param>
    /// <param name="E">Eccentric Anomaly</param>
    /// <param name="e">Eccentricity</param>
    /// <returns>0 if E is correct</returns>
    public static double KeplersEquation(double M, double E, double e)
    {
        // Original version: M = E - e sin(E)
        // Rearranged version: we want to guess for E; this returns 0 if E is correct
        return E - (e * Math.Sin(E)) - M;
    }

    // Differentiated version, for use in Newton's method
    public static double KeplersEquation_dE(double E, double e)
    {
        return 1 - (e * Math.Cos(E));
    }

    public static double GetEccentricAnomaly(double M, double e)
    {
        double targetAccuracy = 1e-9;
        int maxIterations = 1000;

        // Starting guess
        double E = e > 0.8 ? Math.PI : M;

        for(int i = 0; i < maxIterations; i++)
        {
            // NEWTON'S METHOD
            // x_n+1 = x_n + f(x_n)/f'(x_n)
            double nextValue = E - (KeplersEquation(M, E, e) / KeplersEquation_dE(E, e));
            double diff = Math.Abs(E - nextValue);
            E = nextValue;

            if (diff < targetAccuracy)
                break;
        }
        
        return E;
    }
}