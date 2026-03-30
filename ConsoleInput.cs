using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text.Json;

namespace UMF_lab2;

public record InteractiveSetup(IReadOnlyList<double> Nodes, int BasisOrder, 
    ICoefficients Coefficients, IBoundaryConditions BoundaryConditions);

public class SimpleCoefficients : ICoefficients
{
    public double LambdaConst { get; }
    public double SigmaA { get; }
    public double SigmaB { get; }
    public double FConst { get; }
    public double ThetaConst { get; }

    public SimpleCoefficients(double lambdaConst, double sigmaA, double sigmaB, double fConst, double thetaConst = 0.0)
    {
        LambdaConst = lambdaConst;
        SigmaA = sigmaA;
        SigmaB = sigmaB;
        FConst = fConst;
        ThetaConst = thetaConst;
    }

    public double Lambda(double u, double x, double t = 0) 
        => LambdaConst;

    public double Sigma(double u, double x, double dudx, double t = 0) 
        => SigmaA + SigmaB * Math.Abs(dudx);

    public double SigmaDuDt(double u, double x, double t = 0) 
        => SigmaB;

    public double F(double u, double x, double t = 0) 
        => FConst;

    public double Theta(double u, double x, double t = 0) 
        => ThetaConst;
}

public class SimpleBoundaryConditions : IBoundaryConditions
{
    public SimpleBoundaryConditions(BoundaryType leftType, BoundaryType rightType,
                                    double leftValue, double rightValue,
                                    Func<double, double, double>? thetaFunc = null)
    {
        LeftType = leftType;
        RightType = rightType;
        LeftValue = leftValue;
        RightValue = rightValue;
        Theta = thetaFunc ?? ((u, x) => 0.0);
    }

    public BoundaryType LeftType { get; }
    public BoundaryType RightType { get; }

    public double LeftValue { get; }
    public double RightValue { get; }

    public Func<double, double, double> Theta { get; }

    public bool IsBoundaryNode(int nodeId, int totalNodes) => nodeId == 0 || nodeId == totalNodes - 1;
}
public class TestDefinition
{
    public string? Name { get; set; }
    public double[]? Nodes { get; set; }
    public int BasisOrder { get; set; }
    public CoeffDef Coefficients { get; set; } = new CoeffDef();
    public BcDef BoundaryConditions { get; set; } = new BcDef();

    public class CoeffDef
    {
        public double Lambda { get; set; }
        public double SigmaA { get; set; }
        public double SigmaB { get; set; }
        public double F { get; set; }
        public double Theta { get; set; }
    }

    public class BcDef
    {
        public string? LeftType { get; set; }
        public string? RightType { get; set; }
        public double LeftValue { get; set; }
        public double RightValue { get; set; }
    }
}

public static class ConsoleInput
{
    // Выбор теста из папки 'tests'
    public static InteractiveSetup StartInteractiveSession()
    {
        var tests = LoadTestsFromDirectory("../../../tests");

        Console.WriteLine("Тесты:");
        for (int i = 0; i < tests.Count; i++)
        {
            Console.WriteLine($"[{i + 1}] {tests[i].FileName}");
        }

        Console.WriteLine("Выберите тест: ");
        var s = Console.ReadLine();
        int choice = int.Parse(s);

        var def = tests[choice - 1].Definition;

        var nodes = def.Nodes ?? [0.0, 1.0];
        ICoefficients coeffs = new SimpleCoefficients(def.Coefficients.Lambda, 
            def.Coefficients.SigmaA, def.Coefficients.SigmaB, def.Coefficients.F, 
            def.Coefficients.Theta);
        BoundaryType leftType = ParseBoundaryType(def.BoundaryConditions.LeftType);
        BoundaryType rightType = ParseBoundaryType(def.BoundaryConditions.RightType);
        IBoundaryConditions bc = new SimpleBoundaryConditions(leftType, rightType, 
            def.BoundaryConditions.LeftValue, def.BoundaryConditions.RightValue, 
            (u, x) => def.Coefficients.Theta);

        return new InteractiveSetup(nodes, def.BasisOrder, coeffs, bc);
    }

    private static BoundaryType ParseBoundaryType(string? s)
    {
        if (string.IsNullOrWhiteSpace(s)) return BoundaryType.Dirichlet;
        return s.Trim().ToLowerInvariant() switch
        {
            "d" or "dirichlet" => BoundaryType.Dirichlet,
            "n" or "neumann" => BoundaryType.Neumann,
            "r" or "robin" => BoundaryType.Robin,
            _ => BoundaryType.Dirichlet
        };
    }

    private static List<(string FileName, TestDefinition Definition)> LoadTestsFromDirectory(string path)
    {
        var list = new List<(string, TestDefinition)>();
        var files = Directory.GetFiles(path, "*.json");
        var options = new JsonSerializerOptions { PropertyNameCaseInsensitive = true };
        foreach (var f in files)
        {
            var json = File.ReadAllText(f);
            var def = JsonSerializer.Deserialize<TestDefinition>(json, options);
            list.Add((Path.GetFileName(f), def));
        }

        return list;
    }

}
