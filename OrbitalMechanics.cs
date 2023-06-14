// Shoutouts to Scott Anderson for these! I barely know anything about orbital mechanics LOL
// his repo on his implementation: https://github.com/ScottyRAnderson/Keplerian-Orbits
// his yt vid about this: https://www.youtube.com/watch?v=t89De819YMA (highly underrated btw, check him out)

using System.Numerics;
using System.Text.Json.Serialization;

public struct Orbit
{
    public string Name { get; set; }
    public string ParentName { get; set; }
    public double Apoapsis { get; set; }
    public double Periapsis { get; set; }
    public double Inclination { get; set; }
    public double ArgumentOfPeriapsis { get; set; }
    public double SemiMajorAxis { get { return (Apoapsis + Periapsis) / 2; } }
    public double SemiMinorAxis { get { return Math.Sqrt(Apoapsis * Periapsis); } }
    public double Eccentricity { get { return LinearEccentricity / SemiMajorAxis; } }
    public double LinearEccentricity { get { return SemiMajorAxis - Periapsis; } }
    public Orbit(double apoapsis, double periapsis, double inclination, double argumentOfPeriapsis, string name = "Planet", string parentName = "Sun")
    {
        Name = name;
        ParentName = parentName;
        if (periapsis > apoapsis)
            (apoapsis, periapsis) = (periapsis, apoapsis);
        Apoapsis = apoapsis;
        Periapsis = periapsis;
        Inclination = inclination;
        ArgumentOfPeriapsis = argumentOfPeriapsis;
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
    public (double, double, double) GetTimeAdjustedPosition(double t)
    {
        t = (t % 1 + 1) % 1; // modulo that supports negatives
        double meanAnomaly = t * Math.Tau;
        double eccentricAnomaly = GetEccentricAnomaly(meanAnomaly);

        double x = SemiMajorAxis * (Math.Cos(eccentricAnomaly) - Eccentricity);
        double y = SemiMinorAxis * Math.Sin(eccentricAnomaly);

        return TiltPosition(x, y);
    }

    public (double, double, double) GetPositionAtAngle(double revolutions)
    {
        revolutions = (revolutions % 1 + 1) % 1; // modulo that supports negatives
        double trueAnomaly = revolutions * Math.Tau;

        double x = SemiMajorAxis * (Math.Cos(trueAnomaly) - Eccentricity);
        double y = SemiMinorAxis * Math.Sin(trueAnomaly);

        return TiltPosition(x, y);
    }

    public (double, double, double) TiltPosition(double x, double y)
    {
        double inclineRadians = Inclination * Math.PI / 180.0;
        double periapsisRadians = ArgumentOfPeriapsis * Math.PI / 180.0;

        double newX = x * Math.Cos(inclineRadians) - y * Math.Sin(inclineRadians);
        double newY = x * Math.Sin(inclineRadians) + y * Math.Cos(inclineRadians);

        double newZ = -newX * Math.Sin(periapsisRadians);
        newX *= Math.Cos(periapsisRadians);

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
    { return $"{ParentName}>{Name}: {Apoapsis}x{Periapsis}"; }
    public string ToStringDetail()
    { return $"{ParentName}>{Name}: {Apoapsis}x{Periapsis}, {Inclination}°/{ArgumentOfPeriapsis}°"; }
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

    public static double SolveKepler(double M, double e)
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

    /// <param name="ap">Apoapsis</param>
    /// <param name="pe">Periapsis</param>
    /// <param name="t">Time</param>
    /// <returns>a 2D vector</returns>
    public static (double, double, double) ConvertOrbitToXYZ(double ap, double pe, double t, double inclination, double argOfPe)
    {
        t %= 1;
        double semiMajorAxis = (ap + pe) / 2;
        double semiMinorAxis = Math.Sqrt(ap + pe);

        double meanAnomaly = t * Math.Tau;
        double linearEccentricity = semiMajorAxis - pe;
        double eccentricity = linearEccentricity / semiMajorAxis;

        double eccentricAnomaly = SolveKepler(meanAnomaly, eccentricity);

        double x = semiMajorAxis * (Math.Cos(eccentricAnomaly) - eccentricity);
        double y = semiMinorAxis * Math.Sin(eccentricAnomaly);

        return CalculateInclinedPosition(x, y, inclination, argOfPe);
    }

    public static (double, double, double) ConvertOrbitAngleToXYZ(double ap, double pe, double revolutions, double inclination, double argOfPe)
    {
        revolutions %= 1;
        double semiMajorAxis = (ap + pe) / 2;
        double semiMinorAxis = Math.Sqrt(ap + pe);
        double trueAnomaly = revolutions * Math.Tau;

        double linearEccentricity = semiMajorAxis - pe;
        double eccentricity = linearEccentricity / semiMajorAxis;

        double x = semiMajorAxis * (Math.Cos(trueAnomaly) - eccentricity);
        double y = semiMinorAxis * Math.Sin(trueAnomaly);

        return CalculateInclinedPosition(x, y, inclination, argOfPe);
    }

    public static (double, double, double) CalculateInclinedPosition(double x, double y, double inclination, double argOfPe)
    {

        double inclineRadians = inclination * Math.PI / 180.0;
        double periapsisRadians = argOfPe * Math.PI / 180.0;

        double newX = x * Math.Cos(inclineRadians) - y * Math.Sin(inclineRadians);
        double newY = x * Math.Sin(inclineRadians) + y * Math.Cos(inclineRadians);
        // double Z = 0.0; // Assuming z is initially 0

        // Rotate around the y-axis (argument of periapsis)
        double newZ = -newX * Math.Sin(periapsisRadians); // + Z * Math.Cos(periapsisRadians);
        newX *= Math.Cos(periapsisRadians); // + Z * Math.Sin(periapsisRadians);

        return (newX, newY, newZ);
    }
}
