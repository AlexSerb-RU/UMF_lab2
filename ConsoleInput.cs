using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text.Json;

namespace UMF_lab2;

public record InteractiveSetup(
    IReadOnlyList<double> SpaceNodes,
    IReadOnlyList<double> TimePoints,
    ICoefficients Coefficients,
    IExactSolution ExactSolution,
    IBoundaryConditions BoundaryConditions,
    IInitialCondition InitialCondition,
    SolverDefinition Solver,
    string TestName);

public record CoefficientRegion(
    double XFrom,
    double XTo,
    double Lambda,
    double Sigma0,
    double Sigma1,
    double SigmaEps);

public class SolverDefinition
{
    public string PrimaryMethod { get; set; } = "both";
    public double Tolerance { get; set; } = 1e-8;
    public int PicardMaxIterations { get; set; } = 300;
    public int NewtonMaxIterations { get; set; } = 50;
    public double PicardRelaxation { get; set; } = 1.0;
    public double NewtonRelaxation { get; set; } = 1.0;
}

public class MeshDefinition
{
    public AxisMeshDefinition Space { get; set; } = new();
    public AxisMeshDefinition Time { get; set; } = new();
}

public class AxisMeshDefinition
{
    public double[]? Nodes { get; set; }
    public double[]? ReferencePoints { get; set; }
    public IntervalSplitDefinition[]? Splits { get; set; }
}

public class IntervalSplitDefinition
{
    public int Count { get; set; }
    public double Q { get; set; } = 1.0;
}

public class TestDefinition
{
    public string Name { get; set; } = "unnamed_test";
    public MeshDefinition Mesh { get; set; } = new();
    public string Dataset { get; set; } = "";
    public SolverDefinition Solver { get; set; } = new();
}

public class PiecewiseVariant11Coefficients : ICoefficients
{
    public class Region
    {
        public double XFrom { get; init; }
        public double XTo { get; init; }
        public double Lambda { get; init; } = 1.0;
        public double Sigma0 { get; init; } = 1.0;
        public double Sigma1 { get; init; }
        public double SigmaEps { get; init; } = 1e-6;
        public double F { get; init; }
    }

    private readonly List<Region> _regions;

    public PiecewiseVariant11Coefficients(IEnumerable<Region> regions)
    {
        _regions = [.. regions.OrderBy(r => r.XFrom)];
        if (_regions.Count == 0)
        {
            throw new ArgumentException("At least one coefficient region must be specified.");
        }
        foreach (var region in _regions)
        {
            if (region.XTo <= region.XFrom)
            {
                throw new ArgumentException("Coefficient region must satisfy x_to > x_from.");
            }
        }
    }

    private Region GetRegion(double x)
    {
        foreach (var region in _regions)
        {
            if (x >= region.XFrom - 1e-12 && x <= region.XTo + 1e-12)
            {
                return region;
            }
        }
        return _regions[^1];
    }

    public double Lambda(double x, double t = 0.0) => GetRegion(x).Lambda;

    public double Sigma(double dudx, double x, double t = 0.0)
    {
        var r = GetRegion(x);
        double eps = Math.Max(1e-12, r.SigmaEps);
        return r.Sigma0 + r.Sigma1 * dudx;
    }

    public double DSigmaDDuDx(double dudx, double x, double t = 0.0)
    {
        var r = GetRegion(x);
        double eps = Math.Max(1e-12, r.SigmaEps);
        return r.Sigma1 * dudx / Math.Sqrt(dudx * dudx + eps * eps);
    }

    // МЕНЯТЬ
    public double F(double x, double t = 0.0) => -2;
}

public class SimpleBoundaryConditions(
    BoundaryType leftType,
    BoundaryType rightType,
    Func<double, double> leftValue,
    Func<double, double> rightValue) : IBoundaryConditions
{
    public BoundaryType LeftType { get; } = leftType;
    public BoundaryType RightType { get; } = rightType;

    public double LeftValue(double t) => leftValue(t);
    public double RightValue(double t) => rightValue(t);
}

public class FunctionInitialCondition : IInitialCondition
{
    private readonly Func<double, double> _func;
    public FunctionInitialCondition(Func<double, double> func) => _func = func;
    public double Value(double x) => _func(x);
}

public class ManufacturedSolution : IManufacturedSolution, IExactSolution
{
    private readonly Func<double, double, double> _u;
    private readonly Func<double, double, double> _ux;
    private readonly Func<double, double, double> _uxx;
    private readonly Func<double, double, double> _ut;

    public ManufacturedSolution(
        Func<double, double, double> u,
        Func<double, double, double> ux,
        Func<double, double, double> uxx,
        Func<double, double, double> ut)
    {
        _u = u;
        _ux = ux;
        _uxx = uxx;
        _ut = ut;
    }

    public double U(double x, double t) => _u(x, t);
    public double Ux(double x, double t) => _ux(x, t);
    public double Uxx(double x, double t) => _uxx(x, t);
    public double Ut(double x, double t) => _ut(x, t);
    public double Value(double x, double t) => _u(x, t);
}

/// <summary>
/// Класс коэффициентов, который строит F по аналитическому решению, заданному в виде IManufacturedSolution. 
/// Позволяет тестировать численный метод на искусственно созданных решениях. 
/// В данном случае F(x,t) = Sigma(u_x, x, t) * u_t - Lambda(x, t) * u_xx, где u - заданное аналитическое решение. 
/// Коэффициенты Sigma и Lambda могут быть заданы как кусочно-постоянные функции от x (через CoefficientRegion), 
/// а также могут зависеть от u_x (для Sigma).
/// </summary>
public class ManufacturedCoefficients : ICoefficients
{
    private readonly IManufacturedSolution _solution;
    private readonly List<CoefficientRegion> _regions;

    public ManufacturedCoefficients(IManufacturedSolution solution, IEnumerable<CoefficientRegion> regions)
    {
        _solution = solution;
        _regions = regions.OrderBy(r => r.XFrom).ToList();
    }

    private CoefficientRegion Region(double x)
    {
        foreach (var r in _regions)
            if (x >= r.XFrom - 1e-12 && x <= r.XTo + 1e-12)
                return r;
        return _regions[^1];
    }

    public double Lambda(double x, double t = 0.0) => Region(x).Lambda;

    public double Sigma(double dudx, double x, double t = 0.0)
    {
        var r = Region(x);
        return r.Sigma0 + r.Sigma1 * Math.Sqrt(dudx * dudx + r.SigmaEps * r.SigmaEps);
    }

    public double DSigmaDDuDx(double dudx, double x, double t = 0.0)
    {
        var r = Region(x);
        return r.Sigma1 * dudx / Math.Sqrt(dudx * dudx + r.SigmaEps * r.SigmaEps);
    }

    public double F(double x, double t = 0.0)
    {
        var ux = _solution.Ux(x, t);
        var ut = _solution.Ut(x, t);
        var uxx = _solution.Uxx(x, t);

        return Sigma(ux, x, t) * ut - Lambda(x, t) * uxx;
    }
}

public static class DatasetFactory
{
    public static (IExactSolution exact, ICoefficients coeffs, IBoundaryConditions bc, IInitialCondition ic)
        Create(string datasetName, double xLeft, double xRight)
    {
        return datasetName.Trim().ToLowerInvariant() switch
        {
            "x2_const_lambda_sigma" => BuildX2Const(xLeft, xRight),
            "x_plus_t_const_lambda_sigma" => BuildXPlusTConst(xLeft, xRight),
            "sinxt_const_lambda_sigma" => BuildSinXtConst(xLeft, xRight),
            _ => throw new InvalidOperationException($"Неизвестный набор данных: {datasetName}")
        };
    }

    private static (IExactSolution, ICoefficients, IBoundaryConditions, IInitialCondition)
        BuildX2Const(double xLeft, double xRight)
    {
        var exact = new ManufacturedSolution(
            (x, t) => x * x,
            (x, t) => 2.0 * x,
            (x, t) => 2.0,
            (x, t) => 0.0);

        var coeffs = new ManufacturedCoefficients(exact, new[]
        {
            new CoefficientRegion(xLeft, xRight, Lambda: 1.0, Sigma0: 1.0, Sigma1: 0.0, SigmaEps: 1e-6)
        });

        var bc = new SimpleBoundaryConditions(
            BoundaryType.Dirichlet,
            BoundaryType.Dirichlet,
            t => exact.U(xLeft, t),
            t => exact.U(xRight, t));

        var ic = new FunctionInitialCondition(x => exact.U(x, 0.0));

        return (exact, coeffs, bc, ic);
    }

    private static (IExactSolution, ICoefficients, IBoundaryConditions, IInitialCondition)
        BuildXPlusTConst(double xLeft, double xRight)
    {
        var exact = new ManufacturedSolution(
            (x, t) => x + t,
            (x, t) => 1.0,
            (x, t) => 0.0,
            (x, t) => 1.0);

        var coeffs = new ManufacturedCoefficients(exact, new[]
        {
            new CoefficientRegion(xLeft, xRight, Lambda: 1.0, Sigma0: 1.0, Sigma1: 0.0, SigmaEps: 1e-6)
        });

        var bc = new SimpleBoundaryConditions(
            BoundaryType.Dirichlet,
            BoundaryType.Dirichlet,
            t => exact.U(xLeft, t),
            t => exact.U(xRight, t));

        var ic = new FunctionInitialCondition(x => exact.U(x, 0.0));

        return (exact, coeffs, bc, ic);
    }

    private static (IExactSolution, ICoefficients, IBoundaryConditions, IInitialCondition)
        BuildSinXtConst(double xLeft, double xRight)
    {
        var exact = new ManufacturedSolution(
            (x, t) => Math.Sin(x) * Math.Exp(t),
            (x, t) => Math.Cos(x) * Math.Exp(t),
            (x, t) => -Math.Sin(x) * Math.Exp(t),
            (x, t) => Math.Sin(x) * Math.Exp(t));

        var coeffs = new ManufacturedCoefficients(exact, new[]
        {
            new CoefficientRegion(xLeft, xRight, Lambda : 1.0, Sigma0 : 1.0, Sigma1 : 0.0, SigmaEps : 1e-6)
        });

        var bc = new SimpleBoundaryConditions(
            BoundaryType.Dirichlet,
            BoundaryType.Dirichlet,
            t => exact.U(xLeft, t),
            t => exact.U(xRight, t));

        var ic = new FunctionInitialCondition(x => exact.U(x, 0.0));

        return (exact, coeffs, bc, ic);
    }
}

public class FunctionExactSolution(Func<double, double, double> func) : IExactSolution
{
    public double Value(double x, double t) => func(x, t);
}

public static class ConsoleInput
{
    public static InteractiveSetup StartInteractiveSession()
    {
        const string testsDir = "../../../tests";
        if (!Directory.Exists(testsDir))
        {
            throw new DirectoryNotFoundException("Не найдена папка tests");
        }

        var files = Directory.GetFiles(testsDir, "*.json").OrderBy(Path.GetFileName).ToArray();
        if (files.Length == 0)
        {
            throw new InvalidOperationException("В папке tests нет JSON-файлов.");
        }

        Console.WriteLine("Доступные тесты:");
        for (int i = 0; i < files.Length; i++)
        {
            Console.WriteLine($"[{i + 1}] {Path.GetFileName(files[i])}");
        }

        Console.Write("Выберите тест: ");
        if (!int.TryParse(Console.ReadLine(), NumberStyles.Integer, CultureInfo.InvariantCulture, out int choice) ||
            choice < 1 || choice > files.Length)
        {
            throw new InvalidOperationException("Некорректный номер теста.");
        }

        string json = File.ReadAllText(files[choice - 1]);
        var test = JsonSerializer.Deserialize<TestDefinition>(json, new JsonSerializerOptions
        {
            PropertyNameCaseInsensitive = true
        }) ?? throw new InvalidOperationException("Не удалось прочитать JSON-тест.");

        var spaceNodes = BuildAxisNodes(test.Mesh.Space, "space");
        var timePoints = BuildAxisNodes(test.Mesh.Time, "time");
        ValidateIncreasing(timePoints, "TimePoints");
        ValidateIncreasing(spaceNodes, "SpaceNodes");

        var (exact, coeffs, bc, ic) = DatasetFactory.Create(test.Dataset, spaceNodes[0], spaceNodes[^1]);

        return new InteractiveSetup(
            spaceNodes,
            timePoints,
            coeffs,
            exact,
            bc,
            ic,
            test.Solver,
            test.Name);
    }

    private static IReadOnlyList<double> BuildAxisNodes(AxisMeshDefinition axis, string axisName)
    {
        if (axis.Nodes is { Length: > 1 })
        {
            ValidateIncreasing(axis.Nodes, axisName + ".nodes");
            return axis.Nodes;
        }

        if (axis.ReferencePoints is not { Length: > 1 })
        {
            throw new InvalidOperationException($"Для оси {axisName} должны быть заданы либо nodes, либо reference_points.");
        }
        if (axis.Splits is null || axis.Splits.Length != axis.ReferencePoints.Length - 1)
        {
            throw new InvalidOperationException($"Для оси {axisName} количество splits должно совпадать с количеством интервалов reference_points.");
        }

        ValidateIncreasing(axis.ReferencePoints, axisName + ".reference_points");
        var result = new List<double> { axis.ReferencePoints[0] };
        for (int i = 0; i < axis.ReferencePoints.Length - 1; i++)
        {
            double left = axis.ReferencePoints[i];
            double right = axis.ReferencePoints[i + 1];
            var split = axis.Splits[i];
            if (split.Count <= 0)
            {
                throw new InvalidOperationException($"Для оси {axisName} count должен быть положительным.");
            }
            var segment = BuildSegment(left, right, split.Count, split.Q);
            for (int k = 1; k < segment.Count; k++)
            {
                result.Add(segment[k]);
            }
        }
        return result;
    }

    private static List<double> BuildSegment(double left, double right, int count, double q)
    {
        double length = right - left;
        if (length <= 0.0) throw new InvalidOperationException("Интервал разбиения должен иметь положительную длину.");
        var nodes = new List<double>(count + 1) { left };

        if (Math.Abs(q - 1.0) < 1e-14)
        {
            double h = length / count;
            for (int i = 1; i <= count; i++)
            {
                nodes.Add(left + i * h);
            }
            nodes[^1] = right;
            return nodes;
        }

        double sum = (Math.Pow(q, count) - 1.0) / (q - 1.0);
        double h0 = length / sum;
        double acc = left;
        for (int i = 1; i <= count; i++)
        {
            acc += h0 * Math.Pow(q, i - 1);
            nodes.Add(acc);
        }
        nodes[^1] = right;
        return nodes;
    }

    private static void ValidateIncreasing(IReadOnlyList<double> values, string name)
    {
        if (values.Count < 2) throw new InvalidOperationException($"{name}: требуется как минимум две точки.");
        for (int i = 1; i < values.Count; i++)
        {
            if (values[i] <= values[i - 1]) throw new InvalidOperationException($"{name}: точки должны строго возрастать.");
        }
    }

    private static BoundaryType ParseBoundaryType(string s)
        => s.Trim().ToLowerInvariant() switch
        {
            "dirichlet" or "d" => BoundaryType.Dirichlet,
            "neumann" or "n" => BoundaryType.Neumann,
            _ => throw new InvalidOperationException($"Неизвестный тип краевого условия: {s}")
        };
}
