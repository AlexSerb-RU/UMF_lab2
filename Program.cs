using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.Diagnostics;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text.Json;

namespace UMF_lab2;

public record LinearBasisFunction(int LocalId, double LeftX, double RightX) : IBasisFunction
{
    public double Value(double x)
    {
        if (x < LeftX || x > RightX) return 0.0;
        double h = RightX - LeftX;
        return LocalId == 0 ? (RightX - x) / h : (x - LeftX) / h;
    }

    public double Derivative(double x)
    {
        if (x < LeftX || x > RightX) return 0.0;
        double h = RightX - LeftX;
        return LocalId == 0 ? -1.0 / h : 1.0 / h;
    }
}

public class FiniteElement1D : IFiniteElement
{
    public FiniteElement1D(int id, int leftNode, int rightNode, double leftX, double rightX)
    {
        if (rightX <= leftX) throw new ArgumentException("Invalid element interval.");
        Id = id;
        GlobalNodes = new ReadOnlyCollection<int>([leftNode, rightNode]);
        Length = rightX - leftX;
        BasisFunctions = new ReadOnlyCollection<IBasisFunction>(
        [
            new LinearBasisFunction(0, leftX, rightX),
            new LinearBasisFunction(1, leftX, rightX)
        ]);
        LeftX = leftX;
        RightX = rightX;
    }

    public int Id { get; }
    public IReadOnlyList<int> GlobalNodes { get; }
    public IReadOnlyList<IBasisFunction> BasisFunctions { get; }
    public double Length { get; }
    public double LeftX { get; }
    public double RightX { get; }
}

public class SpaceMesh : ISpaceMesh
{
    private readonly List<double> _nodes;
    private readonly List<IFiniteElement> _elements;

    public SpaceMesh(IReadOnlyList<double> nodes)
    {
        if (nodes.Count < 2) throw new ArgumentException("Need at least two nodes.");
        _nodes = nodes.ToList();
        for (int i = 1; i < _nodes.Count; i++)
        {
            if (_nodes[i] <= _nodes[i - 1]) throw new ArgumentException("Nodes must be strictly increasing.");
        }

        _elements = new List<IFiniteElement>(_nodes.Count - 1);
        for (int i = 0; i < _nodes.Count - 1; i++)
        {
            _elements.Add(new FiniteElement1D(i, i, i + 1, _nodes[i], _nodes[i + 1]));
        }
    }

    public IReadOnlyList<double> Nodes => _nodes;
    public IReadOnlyList<IFiniteElement> Elements => _elements;
    public int NodeCount => _nodes.Count;
    public int ElementCount => _elements.Count;
    public double LeftBoundary => _nodes[0];
    public double RightBoundary => _nodes[^1];
    public double GetNodeCoordinate(int nodeId) => _nodes[nodeId];
    public IFiniteElement GetElement(int elementId) => _elements[elementId];
}

public class TimeMesh : ITimeMesh
{
    private readonly List<double> _timePoints;

    public TimeMesh(IReadOnlyList<double> timePoints)
    {
        if (timePoints.Count < 2) throw new ArgumentException("Need at least two time points.");
        _timePoints = timePoints.ToList();
        for (int i = 1; i < _timePoints.Count; i++)
        {
            if (_timePoints[i] <= _timePoints[i - 1])
            {
                throw new ArgumentException("Time points must be strictly increasing.");
            }
        }
    }

    public IReadOnlyList<double> TimePoints => _timePoints;
    public int NodeCount => _timePoints.Count;
    public int TimeStepCount => _timePoints.Count - 1;
    public double TimeStart => _timePoints[0];
    public double TimeEnd => _timePoints[^1];
    public double GetTimeStep(int stepIndex) => _timePoints[stepIndex + 1] - _timePoints[stepIndex];
    public double GetTimePoint(int stepIndex) => _timePoints[stepIndex];
}

public class BandedMatrix
{
    private readonly double[] _diag;
    private readonly double[] _lower;
    private readonly double[] _upper;

    public BandedMatrix(int size)
    {
        if (size <= 0) throw new ArgumentOutOfRangeException(nameof(size));
        Size = size;
        _diag = new double[size];
        _lower = new double[Math.Max(0, size - 1)];
        _upper = new double[Math.Max(0, size - 1)];
    }

    public int Size { get; }

    public double this[int i, int j]
    {
        get
        {
            ValidateIndices(i, j);
            if (i == j) return _diag[i];
            if (i == j + 1) return _lower[j];
            if (j == i + 1) return _upper[i];
            return 0.0;
        }
        set
        {
            ValidateIndices(i, j);
            if (i == j) _diag[i] = value;
            else if (i == j + 1) _lower[j] = value;
            else if (j == i + 1) _upper[i] = value;
            else if (Math.Abs(value) > 1e-14) throw new InvalidOperationException("Only tridiagonal storage is supported.");
        }
    }

    public void Add(int i, int j, double value) => this[i, j] = this[i, j] + value;

    public void Clear()
    {
        Array.Clear(_diag, 0, _diag.Length);
        Array.Clear(_lower, 0, _lower.Length);
        Array.Clear(_upper, 0, _upper.Length);
    }

    private void ValidateIndices(int i, int j)
    {
        if (i < 0 || i >= Size || j < 0 || j >= Size) throw new ArgumentOutOfRangeException();
    }
}

public class TridiagonalLUSolver : ILinearSolver
{
    public double[] Solve(BandedMatrix matrix, double[] rhs)
    {
        int n = matrix.Size;
        if (rhs.Length != n) throw new ArgumentException("RHS size mismatch.");

        if (n == 1)
        {
            return [rhs[0] / matrix[0, 0]];
        }

        var a = new double[n - 1];
        var b = new double[n];
        var c = new double[n - 1];
        var d = (double[])rhs.Clone();

        for (int i = 0; i < n; i++)
        {
            b[i] = matrix[i, i];
            if (i < n - 1) c[i] = matrix[i, i + 1];
            if (i > 0) a[i - 1] = matrix[i, i - 1];
        }

        for (int i = 1; i < n; i++)
        {
            if (Math.Abs(b[i - 1]) < 1e-16) throw new InvalidOperationException("Zero pivot in tridiagonal LU.");
            double m = a[i - 1] / b[i - 1];
            b[i] -= m * c[i - 1];
            d[i] -= m * d[i - 1];
        }

        var x = new double[n];
        x[n - 1] = d[n - 1] / b[n - 1];
        for (int i = n - 2; i >= 0; i--)
        {
            x[i] = (d[i] - c[i] * x[i + 1]) / b[i];
        }

        return x;
    }
}

public class NonlinearSolverResult
{
    public NonlinearSolverResult(double[] solution, double[] residualHistory, int iterationsCount, bool converged, double finalResidual, double computationTime)
    {
        Solution = solution;
        ResidualHistory = residualHistory;
        IterationsCount = iterationsCount;
        Converged = converged;
        FinalResidual = finalResidual;
        ComputationTime = computationTime;
    }

    public double[] Solution { get; }
    public double[] ResidualHistory { get; }
    public int IterationsCount { get; }
    public bool Converged { get; }
    public double FinalResidual { get; }
    public double ComputationTime { get; }
}

public class SimpleIterationSolver : INonlinearSolver
{
    private readonly ILinearSolver _linearSolver;
    private readonly double _tolerance;
    private readonly int _maxIterations;

    public SimpleIterationSolver(ILinearSolver linearSolver, double tolerance, int maxIterations, double relaxation = 1.0)
    {
        _linearSolver = linearSolver;
        _tolerance = tolerance;
        _maxIterations = maxIterations;
        RelaxationParameter = relaxation;
    }

    public double RelaxationParameter { get; set; }

    public NonlinearSolverResult Solve(INonlinearProblem problem, double[] initialGuess)
    {
        var sw = Stopwatch.StartNew();
        var q = (double[])initialGuess.Clone();
        var history = new List<double>();

        for (int iter = 1; iter <= _maxIterations; iter++)
        {
            problem.ComputeSystem(q, out var a, out var b);
            var x = _linearSolver.Solve(a, b);
            var qNew = new double[q.Length];
            for (int i = 0; i < q.Length; i++)
            {
                qNew[i] = q[i] + RelaxationParameter * (x[i] - q[i]);
            }

            double rr = problem.ComputeRelativeResidual(qNew);
            history.Add(rr);
            q = qNew;
            if (rr < _tolerance)
            {
                sw.Stop();
                return new NonlinearSolverResult(q, [.. history], iter, true, rr, sw.Elapsed.TotalSeconds);
            }
        }

        sw.Stop();
        double finalResidual = history.Count == 0 ? double.PositiveInfinity : history[^1];
        return new NonlinearSolverResult(q, [.. history], _maxIterations, false, finalResidual, sw.Elapsed.TotalSeconds);
    }
}

public class NewtonSolver : INonlinearSolver
{
    private readonly ILinearSolver _linearSolver;
    private readonly double _tolerance;
    private readonly int _maxIterations;

    public NewtonSolver(ILinearSolver linearSolver, double tolerance, int maxIterations, double relaxation = 1.0)
    {
        _linearSolver = linearSolver;
        _tolerance = tolerance;
        _maxIterations = maxIterations;
        RelaxationParameter = relaxation;
    }

    public double RelaxationParameter { get; set; }

    public NonlinearSolverResult Solve(INonlinearProblem problem, double[] initialGuess)
    {
        if (problem is not ILinearizableNonlinearProblem linearizable)
        {
            throw new InvalidOperationException("Newton solver requires ILinearizableNonlinearProblem.");
        }

        var sw = Stopwatch.StartNew();
        var q = (double[])initialGuess.Clone();
        var history = new List<double>();

        for (int iter = 1; iter <= _maxIterations; iter++)
        {
            linearizable.ComputeLinearizeMatrixAndResidual(q, out var lin_matrix, out var residual);
            var minusResidual = residual.Select(v => -v).ToArray();
            var delta = _linearSolver.Solve(lin_matrix, minusResidual);
            for (int i = 0; i < q.Length; i++)
            {
                q[i] += RelaxationParameter * delta[i];
            }

            double rr = problem.ComputeRelativeResidual(q);
            history.Add(rr);
            if (rr < _tolerance)
            {
                sw.Stop();
                return new NonlinearSolverResult(q, history.ToArray(), iter, true, rr, sw.Elapsed.TotalSeconds);
            }
        }

        sw.Stop();
        double finalResidual = history.Count == 0 ? double.PositiveInfinity : history[^1];
        return new NonlinearSolverResult(q, history.ToArray(), _maxIterations, false, finalResidual, sw.Elapsed.TotalSeconds);
    }
}

public class ImplicitEulerTimeIntegrator : ITimeIntegrator
{
    private readonly INonlinearSolver _solver;

    public ImplicitEulerTimeIntegrator(INonlinearSolver solver)
    {
        _solver = solver;
        LastResidualHistory = [];
    }

    public double[] LastResidualHistory { get; private set; }
    public int LastIterationsCount { get; private set; }

    public double[] TimeStep(ILinearizableNonlinearProblem problem, double[] initialGuess)
    {
        var result = _solver.Solve(problem, initialGuess);
        LastResidualHistory = result.ResidualHistory;
        LastIterationsCount = result.IterationsCount;
        return result.Solution;
    }
}

public class ParabolicProblem : ILinearizableNonlinearProblem
{
    private static readonly double[,] LocalMassTemplate = { { 2.0, 1.0 }, { 1.0, 2.0 } };
    private static readonly double[,] LocalStiffTemplate = { { 1.0, -1.0 }, { -1.0, 1.0 } };

    private readonly ISpaceMesh _mesh;
    private readonly ICoefficients _coeffs;
    private readonly IBoundaryConditions _bc;
    private readonly IReadOnlyList<double> _uPrev;
    private readonly double _dt;
    private readonly double _time;

    public ParabolicProblem(ISpaceMesh mesh, ICoefficients coeffs, IBoundaryConditions bc, IReadOnlyList<double> uPrev, double dt, double time)
    {
        if (dt <= 0.0) throw new ArgumentOutOfRangeException(nameof(dt));
        _mesh = mesh;
        _coeffs = coeffs;
        _bc = bc;
        _uPrev = uPrev;
        _dt = dt;
        _time = time;
    }

    public int Size => _mesh.NodeCount;

    public void ComputeSystem(IReadOnlyList<double> q, out BandedMatrix a, out double[] b)
    {
        AssemblePicardSystem(q, out a, out b);
    }

    public double ComputeRelativeResidual(IReadOnlyList<double> q)
    {
        ComputeSystem(q, out var a, out var b);
        var residual = Multiply(a, q);
        double normR = 0.0;
        double normB = 0.0;
        for (int i = 0; i < residual.Length; i++)
        {
            residual[i] -= b[i];
            normR += residual[i] * residual[i];
            normB += b[i] * b[i];
        }

        return Math.Sqrt(normR) / Math.Max(Math.Sqrt(normB), 1e-30);
    }

    public void ComputeLinearizeMatrixAndResidual(IReadOnlyList<double> q, out BandedMatrix lin_matrix, out double[] residual)
    {
        lin_matrix = new BandedMatrix(Size);
        residual = new double[Size];

        for (int e = 0; e < _mesh.ElementCount; e++)
        {
            var element = (FiniteElement1D)_mesh.GetElement(e);
            int i0 = element.GlobalNodes[0];
            int i1 = element.GlobalNodes[1];
            double h = element.Length;
            double xMid = 0.5 * (element.LeftX + element.RightX);
            double lambda = _coeffs.Lambda(xMid, _time);

            double q0 = q[i0];
            double q1 = q[i1];
            double p0 = _uPrev[i0];
            double p1 = _uPrev[i1];
            double dudx = (q1 - q0) / h;
            double sigma = _coeffs.Sigma(dudx, xMid, _time);
            double dsigma = _coeffs.DSigmaDDuDx(dudx, xMid, _time);
            double fMid = _coeffs.F(xMid, _time);

            var k = new double[2, 2];
            var m = new double[2, 2];
            var mBase = new double[2, 2];
            for (int i = 0; i < 2; i++)
            for (int j = 0; j < 2; j++)
            {
                k[i, j] = lambda / h * LocalStiffTemplate[i, j];
                mBase[i, j] = h / 6.0 * LocalMassTemplate[i, j] / _dt;
                m[i, j] = sigma * mBase[i, j];
            }

            var w = new[] { q0 - p0, q1 - p1 };
            var v = new double[2];
            for (int i = 0; i < 2; i++)
            {
                for (int j = 0; j < 2; j++)
                {
                    v[i] += mBase[i, j] * w[j];
                }
            }
            var d = new[] { -1.0 / h, 1.0 / h };

            var localJ = new double[2, 2];
            var localR = new double[2];
            var qLocal = new[] { q0, q1 };
            var prevLocal = new[] { p0, p1 };
            var localF = new[] { fMid * h / 2.0, fMid * h / 2.0 };

            for (int i = 0; i < 2; i++)
            {
                for (int j = 0; j < 2; j++)
                {
                    localJ[i, j] = k[i, j] + m[i, j] + dsigma * v[i] * d[j];
                }

                double stiffPart = 0.0;
                double massPart = 0.0;
                for (int j = 0; j < 2; j++)
                {
                    stiffPart += k[i, j] * qLocal[j];
                    massPart += m[i, j] * (qLocal[j] - prevLocal[j]);
                }
                localR[i] = stiffPart + massPart - localF[i];
            }

            AssembleLocalSystem(lin_matrix, residual, i0, i1, localJ, localR);
        }

        ApplyBoundaryConditionsToLinearizeMatrixAndResidual(q, lin_matrix, residual);
    }

    private void AssemblePicardSystem(IReadOnlyList<double> q, out BandedMatrix a, out double[] b)
    {
        a = new BandedMatrix(Size);
        b = new double[Size];

        for (int e = 0; e < _mesh.ElementCount; e++)
        {
            var element = (FiniteElement1D)_mesh.GetElement(e);
            int i0 = element.GlobalNodes[0];
            int i1 = element.GlobalNodes[1];
            double h = element.Length;
            double xMid = 0.5 * (element.LeftX + element.RightX);
            double dudx = (q[i1] - q[i0]) / h;
            double lambda = _coeffs.Lambda(xMid, _time);
            double sigma = _coeffs.Sigma(dudx, xMid, _time);
            double fMid = _coeffs.F(xMid, _time);

            var localA = new double[2, 2];
            var localB = new double[2];
            for (int i = 0; i < 2; i++)
            {
                for (int j = 0; j < 2; j++)
                {
                    localA[i, j] = lambda / h * LocalStiffTemplate[i, j]
                                 + sigma * h / 6.0 * LocalMassTemplate[i, j] / _dt;
                    localB[i] += sigma * h / 6.0 * LocalMassTemplate[i, j] * _uPrev[element.GlobalNodes[j]] / _dt;
                }
                localB[i] += fMid * h / 2.0;
            }

            AssembleLocal(a, b, i0, i1, localA, localB);
        }

        ApplyBoundaryConditionsToSystem(a, b);
    }

    private void ApplyBoundaryConditionsToSystem(BandedMatrix a, double[] b)
    {
        ApplyBoundaryConditionSide(true, a, b);
        ApplyBoundaryConditionSide(false, a, b);
    }

    private void ApplyBoundaryConditionSide(bool isLeft, BandedMatrix a, double[] b)
    {
        int boundary = isLeft ? 0 : Size - 1;
        int adjacent = isLeft ? 1 : Size - 2;
        BoundaryType type = isLeft ? _bc.LeftType : _bc.RightType;
        double value = isLeft ? _bc.LeftValue(_time) : _bc.RightValue(_time);
        double beta = isLeft ? _bc.LeftBeta(_time) : _bc.RightBeta(_time);

        switch (type)
        {
            case BoundaryType.Dirichlet:
                if (Size > 1)
                {
                    b[adjacent] -= a[adjacent, boundary] * value;
                    a[adjacent, boundary] = 0.0;
                    a[boundary, adjacent] = 0.0;
                }
                a[boundary, boundary] = 1.0;
                b[boundary] = value;
                break;

            case BoundaryType.Neumann:
                b[boundary] += isLeft ? -value : value;
                break;

            case BoundaryType.Robin:
                a[boundary, boundary] += beta;
                b[boundary] += beta * value;
                break;
        }
    }

    private void ApplyBoundaryConditionsToLinearizeMatrixAndResidual(IReadOnlyList<double> q, BandedMatrix lin_matrix, double[] residual)
    {
        ApplyBoundaryForNewtonSide(true, q, lin_matrix, residual);
        ApplyBoundaryForNewtonSide(false, q, lin_matrix, residual);
    }

    private void ApplyBoundaryForNewtonSide(bool isLeft, IReadOnlyList<double> q, BandedMatrix lin_matrix, double[] residual)
    {
        int boundary = isLeft ? 0 : Size - 1;
        int adjacent = isLeft ? 1 : Size - 2;
        BoundaryType type = isLeft ? _bc.LeftType : _bc.RightType;
        double value = isLeft ? _bc.LeftValue(_time) : _bc.RightValue(_time);
        double beta = isLeft ? _bc.LeftBeta(_time) : _bc.RightBeta(_time);

        switch (type)
        {
            case BoundaryType.Dirichlet:
                if (Size > 1)
                {
                    lin_matrix[adjacent, boundary] = 0.0;
                    lin_matrix[boundary, adjacent] = 0.0;
                }
                lin_matrix[boundary, boundary] = 1.0;
                residual[boundary] = q[boundary] - value;
                break;
            case BoundaryType.Neumann:
                break;
            case BoundaryType.Robin:
                lin_matrix[boundary, boundary] += beta;
                residual[boundary] += beta * (q[boundary] - value);
                break;
        }
    }

    private static void AssembleLocal(BandedMatrix a, double[] b, int i0, int i1, double[,] localA, double[] localB)
    {
        a.Add(i0, i0, localA[0, 0]);
        a.Add(i0, i1, localA[0, 1]);
        a.Add(i1, i0, localA[1, 0]);
        a.Add(i1, i1, localA[1, 1]);
        b[i0] += localB[0];
        b[i1] += localB[1];
    }

    private static void AssembleLocalSystem(BandedMatrix a, double[] residual, int i0, int i1, double[,] localJ, double[] localR)
    {
        a.Add(i0, i0, localJ[0, 0]);
        a.Add(i0, i1, localJ[0, 1]);
        a.Add(i1, i0, localJ[1, 0]);
        a.Add(i1, i1, localJ[1, 1]);
        residual[i0] += localR[0];
        residual[i1] += localR[1];
    }

    private static double[] Multiply(BandedMatrix a, IReadOnlyList<double> q)
    {
        int n = a.Size;
        var result = new double[n];
        for (int i = 0; i < n; i++)
        {
            result[i] += a[i, i] * q[i];
            if (i > 0) result[i] += a[i, i - 1] * q[i - 1];
            if (i < n - 1) result[i] += a[i, i + 1] * q[i + 1];
        }
        return result;
    }
}

public static class ProgramRunner
{
    public static void Main()
    {
        var setup = ConsoleInput.StartInteractiveSession();
        var spaceMesh = new SpaceMesh(setup.SpaceNodes);
        var timeMesh = new TimeMesh(setup.TimePoints);

        var linearSolver = new TridiagonalLUSolver();
        var picard = new SimpleIterationSolver(
            linearSolver,
            setup.Solver.Tolerance,
            setup.Solver.PicardMaxIterations,
            setup.Solver.PicardRelaxation);
        var newton = new NewtonSolver(
            linearSolver,
            setup.Solver.Tolerance,
            setup.Solver.NewtonMaxIterations,
            setup.Solver.NewtonRelaxation);

        var picardIntegrator = new ImplicitEulerTimeIntegrator(picard);
        var newtonIntegrator = new ImplicitEulerTimeIntegrator(newton);

        var u0 = new double[spaceMesh.NodeCount];
        for (int i = 0; i < u0.Length; i++)
        {
            u0[i] = setup.InitialCondition.Value(spaceMesh.GetNodeCoordinate(i));
        }

        object? picardHistory = null;
        object? newtonHistory = null;
        string mode = setup.Solver.PrimaryMethod.Trim().ToLowerInvariant();

        if (mode is "both" or "picard")
        {
            picardHistory = SolveAllLayers("Picard", spaceMesh, timeMesh, setup.Coefficients, setup.BoundaryConditions, picardIntegrator, u0);
        }
        if (mode is "both" or "newton")
        {
            newtonHistory = SolveAllLayers("Newton", spaceMesh, timeMesh, setup.Coefficients, setup.BoundaryConditions, newtonIntegrator, u0);
        }

        var output = new
        {
            TestName = setup.TestName,
            SpaceNodes = spaceMesh.Nodes,
            TimePoints = timeMesh.TimePoints,
            Picard = picardHistory,
            Newton = newtonHistory
        };

        string outputFile = $"../../../results/{setup.TestName}_results.json";
        File.WriteAllText(outputFile, JsonSerializer.Serialize(output, new JsonSerializerOptions { WriteIndented = true }));
        Console.WriteLine($"Результаты сохранены в {outputFile}");
    }

    private static object SolveAllLayers(
        string name,
        ISpaceMesh spaceMesh,
        ITimeMesh timeMesh,
        ICoefficients coeffs,
        IBoundaryConditions bc,
        ITimeIntegrator integrator,
        double[] u0)
    {
        var solutions = new List<double[]> { (double[])u0.Clone() };
        var iterations = new List<int> { 0 };
        var residuals = new List<double[]> { Array.Empty<double>() };
        var current = (double[])u0.Clone();

        Console.WriteLine($"Метод {name}");
        for (int step = 0; step < timeMesh.TimeStepCount; step++)
        {
            double tPrev = timeMesh.GetTimePoint(step);
            double tNext = timeMesh.GetTimePoint(step + 1);
            double dt = timeMesh.GetTimeStep(step);

            var problem = new ParabolicProblem(spaceMesh, coeffs, bc, current, dt, tNext);
            current = integrator.TimeStep(problem, (double[])current.Clone());
            solutions.Add((double[])current.Clone());
            iterations.Add(integrator.LastIterationsCount);
            residuals.Add((double[])integrator.LastResidualHistory.Clone());

            double lastResidual = integrator.LastResidualHistory.Length == 0 ? 0.0 : integrator.LastResidualHistory[^1];
            Console.WriteLine(
                $"  слой {step + 1}: t_prev={tPrev.ToString("F6", CultureInfo.InvariantCulture)}, t={tNext.ToString("F6", CultureInfo.InvariantCulture)}, dt={dt.ToString("F6", CultureInfo.InvariantCulture)}, итераций={integrator.LastIterationsCount}, невязка={lastResidual:E3}");
        }

        return new
        {
            Method = name,
            Solutions = solutions,
            Iterations = iterations,
            ResidualHistories = residuals
        };
    }
}
