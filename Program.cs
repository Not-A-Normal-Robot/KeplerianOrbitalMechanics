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

        double nearestApproach = GetNearestApproachDistance(orbits[0], orbits[1]);
        Console.WriteLine($"Nearest approach of {nearestApproach} calculated");
    }

    static Orbit[] GetOrbitsFromJSON()
    {
        string filePath = Path.Combine(AppDomain.CurrentDomain.BaseDirectory, "Orbits.json");
        string jsonContent = File.ReadAllText(filePath);

        return JsonSerializer.Deserialize<OrbitArray>(jsonContent).orbits;
    }
    static double GetNearestApproachDistance(Orbit a, Orbit b, double targetRange = 1e-3, (double, double, double)[]? aSamples = null, (double, double, double)[]? bSamples = null)
    {

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