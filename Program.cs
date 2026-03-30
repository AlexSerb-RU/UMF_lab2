using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Text.Json;

namespace UMF_lab2;

public record LinearBasisFunction(int LocalId, double LeftX, double RightX) : IBasisFunction
{
    public double Value(double x)
    {
        if (x < LeftX || x > RightX) return 0.0;
        
        if (LocalId == 0)
            return (RightX - x) / (RightX - LeftX);
        else
            return (x - LeftX) / (RightX - LeftX);
    }

    public double Derivative(double x)
    {
        if (x < LeftX || x > RightX) return 0.0;

        if (LocalId == 0)
            return -1.0 / (RightX - LeftX);
        else
            return 1.0 / (RightX - LeftX);
    }
}

public class QuadraticBasisFunction : IBasisFunction
{
    public QuadraticBasisFunction(int localId, double x1, double x2, double x3)
    {
        // Инициализация
    }

    public double Value(double x)
    {
        // Реализация
        throw new NotImplementedException();
    }

    public double Derivative(double x)
    {
        // Реализация
        throw new NotImplementedException();
    }
}

public class FiniteElement1D : IFiniteElement
{
    public int Id { get; }
    public IReadOnlyList<int> GlobalNodes { get; }
    public IReadOnlyList<IBasisFunction> BasisFunctions { get; }
    public double LeftX { get; }
    public double RightX { get; }
    public double Length { get; }

    public FiniteElement1D(int id, IReadOnlyList<int> globalNodes, double leftX, double rightX)
    {
        if (globalNodes.Count != 2) throw new ArgumentException("Linear element requires exactly 2 global nodes.");
        if (rightX <= leftX) throw new ArgumentException("rightX must be greater than leftX.");

        RightX = rightX;
        LeftX = leftX;
        Id = id;
        GlobalNodes = new ReadOnlyCollection<int>([.. globalNodes]);
        Length = rightX - leftX;

        var basis = new List<IBasisFunction>
        {
            new LinearBasisFunction(0, leftX, rightX),
            new LinearBasisFunction(1, leftX, rightX)
        };
        BasisFunctions = new ReadOnlyCollection<IBasisFunction>(basis);
    }

    // Локальная матрица жесткости для оператора -d^2/dx^2 (без коэффициентов)
    // K_local = (1 / L) * [ [1, -1], [-1, 1] ]
    public double[,] StiffnessMatrix()
    {
        var k = 1.0 / Length;
        return new double[,] { { k, -k }, { -k, k } };
    }

    // Локальная матрица масс для линейного элемента
    // M_local = (L / 6) * [ [2, 1], [1, 2] ]
    public double[,] MassMatrix()
    {
        var c = Length / 6.0;
        return new double[,] { { 2 * c, 1 * c }, { 1 * c, 2 * c } };
    }
}

public class Mesh : IMesh
{
    private readonly List<double> nodes;
    private readonly List<IFiniteElement> elements;

    public Mesh(IReadOnlyList<double> nodes, IElementFactory? elementFactory, int basisOrder)
    {
        if (nodes == null) throw new ArgumentNullException(nameof(nodes));
        this.nodes = new List<double>(nodes);
        // For now elements are not constructed (assembly uses node-based scheme)
        this.elements = new List<IFiniteElement>();
        Nodes = this.nodes.AsReadOnly();
        Elements = this.elements.AsReadOnly();
        NodeCount = this.nodes.Count;
        ElementCount = Math.Max(0, NodeCount - 1);
    }

    public IReadOnlyList<double> Nodes { get; private set; }
    public IReadOnlyList<IFiniteElement> Elements { get; private set; }
    public int NodeCount { get; private set; }
    public int ElementCount { get; private set; }

    public double GetNodeCoordinate(int nodeId)
    {
        if (nodeId < 0 || nodeId >= nodes.Count) throw new ArgumentOutOfRangeException(nameof(nodeId));
        return nodes[nodeId];
    }

    public IFiniteElement GetElement(int elementId)
    {
        if (elementId < 0 || elementId >= elements.Count) throw new ArgumentOutOfRangeException(nameof(elementId));
        return elements[elementId];
    }
}

public class GalerkinLocalAssembly(IQuadrature quadrature) : ILocalAssembly
{
    private readonly IQuadrature quadrature = quadrature;
    private double[,] localMatrix;
    private double[] localVector;
    private int size;

    public double[,] LocalMatrix => localMatrix;
    public double[] LocalVector => localVector;
    public int Size => size;

    public void Compute(IFiniteElement element, ICoefficients coeffs, IReadOnlyList<double> globalSolution)
    {
        if (element.BasisFunctions.Count != element.GlobalNodes.Count)
            throw new InvalidOperationException("Mismatch between basis functions and global nodes count");

        size = element.BasisFunctions.Count;

        // Инициализация локальной матрицы и вектора
        localMatrix = new double[size, size];
        localVector = new double[size];

        // Информация об элементе
        double leftX, rightX;
        double elementLength = element.Length;
        IReadOnlyList<int> globalNodes = element.GlobalNodes;
        IReadOnlyList<IBasisFunction> basis = element.BasisFunctions;

        // Приведение типа к FiniteElement1D для получения координат
        if (element is FiniteElement1D elem1D)
        {
            leftX = elem1D.LeftX;
            rightX = elem1D.RightX;
        }
        else
        {
            throw new NotSupportedException($"Element type {element.GetType().Name} is not supported");
        }

        // Численное интегрирование по точкам Гаусса
        GaussIntegrate(coeffs, globalSolution, leftX, rightX, elementLength, globalNodes, basis);
    }

    private void GaussIntegrate(ICoefficients coeffs, IReadOnlyList<double> globalSolution, double leftX, double rightX, double elementLength, IReadOnlyList<int> globalNodes, IReadOnlyList<IBasisFunction> basis)
    {
        var quadPoints = quadrature.Points;
        var quadWeights = quadrature.Weights;

        for (int qp = 0; qp < quadPoints.Count; qp++)
        {
            // Преобразование с отрезка [-1, 1] на [leftX, rightX]
            double xiRef = quadPoints[qp];  // Опорный элемент [-1, 1]
            double x = 0.5 * (rightX + leftX) + 0.5 * (rightX - leftX) * xiRef;
            double weight = quadWeights[qp] * 0.5 * elementLength;  // Якобиан

            // Интерполяция решения и его производной в точке квадратуры
            double uAtQp = 0.0;
            double duDxAtQp = 0.0;
            for (int i = 0; i < size; i++)
            {
                uAtQp += globalSolution[globalNodes[i]] * basis[i].Value(x);
                duDxAtQp += globalSolution[globalNodes[i]] * basis[i].Derivative(x);
            }

            // Вычисление коэффициентов в точке квадратуры
            double sigma = coeffs.Sigma(uAtQp, x, duDxAtQp);
            double lambda = coeffs.Lambda(uAtQp, x);
            double f = coeffs.F(uAtQp, x);

            // Вычисление значений базисных функций и их производных
            var phiValues = new double[size];
            var phiDerivValues = new double[size];
            for (int i = 0; i < size; i++)
            {
                phiValues[i] = basis[i].Value(x);
                phiDerivValues[i] = basis[i].Derivative(x);
            }

            // Сборка локальной матрицы жесткости и масс
            for (int i = 0; i < size; i++)
            {
                for (int j = 0; j < size; j++)
                {
                    // Вклад жёсткости
                    double stiffness = sigma * phiDerivValues[i] * phiDerivValues[j];
                    // Вклад массы
                    double mass = lambda * phiValues[i] * phiValues[j];

                    localMatrix[i, j] += (stiffness + mass) * weight;
                }

                // Сборка локального вектора правой части
                localVector[i] += f * phiValues[i] * weight;
            }
        }
    }
}

public class GlobalAssembly : IGlobalAssembly
{
    private BandedMatrix globalMatrix;
    private double[] globalVector;
    private readonly int nodeCount;
    private readonly int bandwidth;

    public BandedMatrix GlobalMatrix => globalMatrix;
    public double[] GlobalVector => globalVector;

    public GlobalAssembly(int nodeCount, int bandwidth = 1)
    {
        this.nodeCount = nodeCount;
        this.bandwidth = bandwidth;  // Учитывает соседние узлы для 1D задачи
        globalMatrix = new BandedMatrix(nodeCount, bandwidth);
        globalVector = new double[nodeCount];
    }

    public void Assemble(IMesh mesh, ILocalAssembly localAssembly, ICoefficients coeffs, IReadOnlyList<double> solution, IBoundaryConditions bc)
    {
        // Очистка глобальных структур
        globalMatrix.Clear();
        Array.Clear(globalVector, 0, globalVector.Length);

        int elementCount = mesh.ElementCount;

        // Проход по всем элементам
        for (int elemId = 0; elemId < elementCount; elemId++)
        {
            IFiniteElement element = mesh.GetElement(elemId);

            // Вычисление локальной матрицы и вектора для текущего элемента
            if (localAssembly != null)
            {
                localAssembly.Compute(element, coeffs, solution);
                var localMatrix = localAssembly.LocalMatrix;
                var localVector = localAssembly.LocalVector;
                int localSize = localAssembly.Size;

                // Перевод локальных индексов в глобальные
                IReadOnlyList<int> globalNodeIndices = element.GlobalNodes;

                // Вставка локальной матрицы в глобальную матрицу
                for (int i = 0; i < localSize; i++)
                {
                    int globalI = globalNodeIndices[i];
                    if (globalI < 0 || globalI >= nodeCount)
                        throw new InvalidOperationException($"Global node index {globalI} out of bounds");

                    for (int j = 0; j < localSize; j++)
                    {
                        int globalJ = globalNodeIndices[j];
                        if (globalJ < 0 || globalJ >= nodeCount)
                            throw new InvalidOperationException($"Global node index {globalJ} out of bounds");

                        // Накопление значений в глобальной матрице
                        globalMatrix[globalI, globalJ] += localMatrix[i, j];
                    }
                }

                // Вставка локального вектора в глобальный вектор
                for (int i = 0; i < localSize; i++)
                {
                    int globalI = globalNodeIndices[i];
                    globalVector[globalI] += localVector[i];
                }
            }
        }

        // Применение граничных условий
        ApplyBoundaryConditions(bc, solution);
    }

    public void ApplyBoundaryConditions(IBoundaryConditions bc, IReadOnlyList<double> solution)
    {
        // Левое краевое условие (узел 0)
        if (bc.LeftType == BoundaryType.Dirichlet)
        {
            // 1-е краевое условие: u(0) = u_left
            // Обнулить строку и установить диагональный элемент = 1
            for (int j = 0; j < nodeCount; j++)
                globalMatrix[0, j] = 0.0;
            globalMatrix[0, 0] = 1.0;
            globalVector[0] = bc.LeftValue;
        }
        else if (bc.LeftType == BoundaryType.Neumann)
        {
            // 2-е краевое условие: sigma * du/dx|_(x=0) = q_left
            globalVector[0] -= bc.LeftValue;
        }
        else if (bc.LeftType == BoundaryType.Robin)
        {
            // Условие Робина: theta * u + sigma * du/dx = g_left
            // Требует специального учёта (обычно в локальной сборке)
        }

        // Правое краевое условие (узел N-1)
        int lastNode = nodeCount - 1;
        if (bc.RightType == BoundaryType.Dirichlet)
        {
            // 1-е краевое условие: u(L) = u_right
            for (int j = 0; j < nodeCount; j++)
                globalMatrix[lastNode, j] = 0.0;
            globalMatrix[lastNode, lastNode] = 1.0;
            globalVector[lastNode] = bc.RightValue;
        }
        else if (bc.RightType == BoundaryType.Neumann)
        {
            // 2-е краевое условие: sigma * du/dx|_(x=L) = q_right
            globalVector[lastNode] += bc.RightValue;
        }
        else if (bc.RightType == BoundaryType.Robin)
        {
            // 3-е краевое
        }
    }

    public void Clear()
    {
        globalMatrix.Clear();
        Array.Clear(globalVector, 0, globalVector.Length);
    }

    public override string ToString()
    {
        return $"GlobalAssembly: NodeCount={nodeCount}, Bandwidth={bandwidth}, MatrixSize={globalMatrix.Size}";
    }
}

public class BandedMatrix
{
    private readonly double[,] data;
    public BandedMatrix(int size, int bandwidth)
    {
        Size = size;
        Bandwidth = Math.Max(0, bandwidth);
        data = new double[size, size];  // Хранится полная матрица
    }

    public double this[int i, int j]
    {
        get
        {
            return data[i, j];
        }
        set
        {
            data[i, j] = value;
        }
    }

    public int Size { get; private set; }
    public int Bandwidth { get; private set; }

    public void Clear()
    {
        Array.Clear(data, 0, data.Length);
    }
}

public class BandedLUSolver : ILinearSolver
{
    public double[] Solve(BandedMatrix matrix, double[] rhs)
    {
        if (matrix == null) throw new ArgumentNullException(nameof(matrix));
        int n = matrix.Size;
        if (rhs == null) throw new ArgumentNullException(nameof(rhs));
        if (rhs.Length != n) throw new ArgumentException("RHS length mismatch", nameof(rhs));

        // Convert to dense array
        var a = new double[n, n];
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                a[i, j] = matrix[i, j];

        var b = new double[n];
        Array.Copy(rhs, b, n);

        // Gaussian elimination with partial pivoting
        for (int k = 0; k < n; k++)
        {
            // find pivot
            int piv = k;
            double max = Math.Abs(a[k, k]);
            for (int i = k + 1; i < n; i++)
            {
                double val = Math.Abs(a[i, k]);
                if (val > max) { max = val; piv = i; }
            }
            if (Math.Abs(a[piv, k]) < 1e-15) throw new InvalidOperationException("Matrix is singular or nearly singular");
            if (piv != k)
            {
                // swap rows k and piv
                for (int j = k; j < n; j++)
                {
                    double tmp = a[k, j]; a[k, j] = a[piv, j]; a[piv, j] = tmp;
                }
                double tb = b[k]; b[k] = b[piv]; b[piv] = tb;
            }

            double akk = a[k, k];
            for (int i = k + 1; i < n; i++)
            {
                double factor = a[i, k] / akk;
                a[i, k] = 0.0;
                for (int j = k + 1; j < n; j++) a[i, j] -= factor * a[k, j];
                b[i] -= factor * b[k];
            }
        }

        // Back substitution
        var x = new double[n];
        for (int i = n - 1; i >= 0; i--)
        {
            double sum = b[i];
            for (int j = i + 1; j < n; j++) sum -= a[i, j] * x[j];
            x[i] = sum / a[i, i];
        }

        return x;
    }
}

public class SimpleIterationSolver : INonlinearSolver
{
    private readonly ILinearSolver linearSolver;
    private readonly double tolerance;
    private readonly int maxIterations;
    public double RelaxationParameter { get; set; }

    public SimpleIterationSolver(ILinearSolver linearSolver, double tolerance, int maxIterations, double relaxation = 1.0)
    {
        this.linearSolver = linearSolver;
        this.tolerance = tolerance > 0 ? tolerance : 1e-8;
        this.maxIterations = Math.Max(1, maxIterations);
        this.RelaxationParameter = relaxation;
    }

    public NonlinearSolverResult Solve(INonlinearProblem problem, double[] initialGuess)
    {
        var sw = Stopwatch.StartNew();
        int n = initialGuess.Length;
        var q = new double[n];
        Array.Copy(initialGuess, q, n);
        var residualHistory = new List<double>();
        bool converged = false;
        double finalResidual = double.PositiveInfinity;

        for (int iter = 0; iter < maxIterations; iter++)
        {
            // Build system A(q) and b(q)
            problem.ComputeSystem(q, out BandedMatrix A, out double[] b);

            // Solve A * x = b
            double[] x;
            try
            {
                x = linearSolver.Solve(A, b);
            }
            catch (Exception ex)
            {
                throw new InvalidOperationException("Linear solver failed", ex);
            }

            // Relaxation: q_new = q + omega*(x - q)
            double omega = RelaxationParameter;
            var qNew = new double[n];
            for (int i = 0; i < n; i++) qNew[i] = q[i] + omega * (x[i] - q[i]);

            // Compute residual R(qNew) = A(qNew)*qNew - b(qNew)
            double[] R = problem.ComputeResidual(qNew);
            double resNorm = 0.0;
            if (R != null && R.Length > 0) resNorm = Math.Sqrt(R.Select(v => v * v).Sum());
            residualHistory.Add(resNorm);
            finalResidual = resNorm;

            if (resNorm <= tolerance)
            {
                converged = true;
                q = qNew;
                sw.Stop();
                return new NonlinearSolverResult(q, residualHistory.ToArray(), iter + 1, true, finalResidual, sw.Elapsed.TotalSeconds);
            }

            // update
            q = qNew;
        }

        sw.Stop();
        return new NonlinearSolverResult(q, residualHistory.ToArray(), maxIterations, converged, finalResidual, sw.Elapsed.TotalSeconds);
    }
}

public class NewtonSolver : INonlinearSolver
{
    public NewtonSolver(ILinearSolver linearSolver, ILinearizationStrategy linearization,
                        double tolerance, int maxIterations, double relaxation = 1.0)
    {
        // Инициализация
    }

    public NonlinearSolverResult Solve(INonlinearProblem problem, double[] initialGuess)
    {
        throw new NotImplementedException();
    }

    public double RelaxationParameter { get; set; }
}

public class NewtonLinearization : ILinearizationStrategy
{
    public NewtonLinearization(IQuadrature quadrature)
    {
        // Инициализация
    }

    public void Linearize(ILocalAssembly localAssembly, IFiniteElement element,
                            ICoefficients coeffs, IReadOnlyList<double> solution,
                            out double[,] linearizedMatrix, out double[] linearizedVector)
    {
        // Реализация: ∂(A(q)q)/∂q и ∂b(q)/∂q
        throw new NotImplementedException();
    }
}

public class GaussQuadrature : IQuadrature
{
    private readonly List<double> points;
    private readonly List<double> weights;

    public IReadOnlyList<double> Points => points.AsReadOnly();
    public IReadOnlyList<double> Weights => weights.AsReadOnly();
    public int Order { get; private set; }

    public GaussQuadrature(int order)
    {
        Order = order;
        points = [];
        weights = [];

        InitializeQuadratureRule(order);
    }

    private void InitializeQuadratureRule(int order)
    {
        switch (order)
        {
            case 1:
                // 1-точечная формула
                points.Add(0.0);
                weights.Add(2.0);
                break;

            case 2:
                // 2-точечная формула
                points.Add(-1.0 / Math.Sqrt(3));
                points.Add(1.0 / Math.Sqrt(3));
                weights.Add(1.0);
                weights.Add(1.0);
                break;

            case 3:
                // 3-точечная формула
                points.Add(-Math.Sqrt(3.0 / 5.0));
                points.Add(0.0);
                points.Add(Math.Sqrt(3.0 / 5.0));
                weights.Add(5.0 / 9.0);
                weights.Add(8.0 / 9.0);
                weights.Add(5.0 / 9.0);
                break;

            case 4:
                // 4-точечная формула
                double sqrt525 = Math.Sqrt(525.0 - 70.0 * Math.Sqrt(30.0)) / 35.0;
                double sqrt525_2 = Math.Sqrt(525.0 + 70.0 * Math.Sqrt(30.0)) / 35.0;
                
                points.Add(-sqrt525_2);
                points.Add(-sqrt525);
                points.Add(sqrt525);
                points.Add(sqrt525_2);
                
                double w1 = (630.0 - 42.0 * Math.Sqrt(30.0)) / 1260.0;
                double w2 = (630.0 + 42.0 * Math.Sqrt(30.0)) / 1260.0;
                
                weights.Add(w1);
                weights.Add(w2);
                weights.Add(w2);
                weights.Add(w1);
                break;

            case 5:
                // 5-точечная формула
                double sqrt245 = Math.Sqrt(245.0 - 14.0 * Math.Sqrt(70.0)) / 21.0;
                double sqrt245_2 = Math.Sqrt(245.0 + 14.0 * Math.Sqrt(70.0)) / 21.0;
                
                points.Add(-sqrt245_2);
                points.Add(-sqrt245);
                points.Add(0.0);
                points.Add(sqrt245);
                points.Add(sqrt245_2);
                
                double w1_5 = (322.0 - 13.0 * Math.Sqrt(70.0)) / 900.0;
                double w2_5 = (322.0 + 13.0 * Math.Sqrt(70.0)) / 900.0;
                
                weights.Add(w1_5);
                weights.Add(w2_5);
                weights.Add(128.0 / 225.0);
                weights.Add(w2_5);
                weights.Add(w1_5);
                break;

            case 6:
                // 6-точечная формула
                double sqrt231 = Math.Sqrt(231.0);
                double sqrt385_1 = Math.Sqrt(385.0 - 14.0 * sqrt231) / 33.0;
                double sqrt385_2 = Math.Sqrt(385.0 + 14.0 * sqrt231) / 33.0;
                
                points.Add(-sqrt385_2);
                points.Add(-sqrt385_1);
                points.Add(-Math.Sqrt(1.0 / 7.0));
                points.Add(Math.Sqrt(1.0 / 7.0));
                points.Add(sqrt385_1);
                points.Add(sqrt385_2);
                
                double w1_6 = (1280.0 - 56.0 * sqrt231) / 8960.0;
                double w2_6 = (1280.0 + 56.0 * sqrt231) / 8960.0;
                
                weights.Add(w1_6);
                weights.Add(w2_6);
                weights.Add(189.0 / 1120.0);
                weights.Add(189.0 / 1120.0);
                weights.Add(w2_6);
                weights.Add(w1_6);
                break;

            default:
                throw new InvalidOperationException($"Order {order} not implemented");
        }
    }

    public int Accuracy => 2 * Order - 1;

    public bool Validate()
    {
        double sumWeights = weights.Sum();
        return Math.Abs(sumWeights - 2.0) < 1e-10;
    }

    // Интегрирование функции f(x) на отрезке [-1, 1] по квадратуре
    public double Integrate(Func<double, double> f)
    {
        double result = 0.0;
        for (int i = 0; i < points.Count; i++)
        {
            result += weights[i] * f(points[i]);
        }
        return result;
    }

    public override string ToString()
    {
        return $"Gauss-Legendre Quadrature: Order={Order}, Points={points.Count}, Accuracy=O(x^{Accuracy})";
    }
}

public class LinearElementFactory : IElementFactory
{
    public IFiniteElement CreateElement(int id, int node1, int node2)
    {
        // Реализация
        throw new NotImplementedException();
    }

    public IBasisFunction CreateLinearBasisFunction(int localId, double x1, double x2)
    {
        // Реализация
        throw new NotImplementedException();
    }

    public IBasisFunction CreateQuadraticBasisFunction(int localId, double x1, double x2, double x3)
    {
        // Реализация (возвращает null для линейных элементов)
        throw new NotImplementedException();
    }
}

public class QuadraticElementFactory : IElementFactory
{
    public IFiniteElement CreateElement(int id, int node1, int node2)
    {
        // Реализация: создание квадратичного элемента с промежуточным узлом
        throw new NotImplementedException();
    }

    public IBasisFunction CreateLinearBasisFunction(int localId, double x1, double x2)
    {
        // Реализация (возвращает null для квадратичных элементов)
        throw new NotImplementedException();
    }

    public IBasisFunction CreateQuadraticBasisFunction(int localId, double x1, double x2, double x3)
    {
        // Реализация
        throw new NotImplementedException();
    }
}

public class NonlinearSolverResult : ISolutionResult
{
    public NonlinearSolverResult(double[] solution, double[] residualHistory,
                                    int iterations, bool converged, double finalResidual,
                                    double computationTime)
    {
        Solution = solution ?? Array.Empty<double>();
        ResidualHistory = residualHistory ?? Array.Empty<double>();
        IterationsCount = iterations;
        Converged = converged;
        FinalResidual = finalResidual;
        ComputationTime = computationTime;
    }

    public double[] Solution { get; private set; }
    public double[] ResidualHistory { get; private set; }
    public int IterationsCount { get; private set; }
    public bool Converged { get; private set; }
    public double FinalResidual { get; private set; }
    public double ComputationTime { get; private set; }
}

public class NonlinearFEMSolver
{
    public NonlinearFEMSolver() { }

    public void ConfigureMesh(IReadOnlyList<double> nodes, int basisOrder) { }
    public void ConfigureCoefficients(ICoefficients coefficients) { }
    public void ConfigureBoundaryConditions(IBoundaryConditions boundaryConditions) { }
    public void ConfigureSolver(INonlinearSolver solver) { }
    public ISolutionResult Solve(double[] initialGuess)
    {
        throw new NotImplementedException();
    }
    public void ExportResults(string filename, ISolutionResult result) { }

    private IMesh mesh;
    private ICoefficients coefficients;
    private IBoundaryConditions boundaryConditions;
    private INonlinearSolver nonlinearSolver;
    private ILinearSolver linearSolver;
    private IGlobalAssembly globalAssembly;
    private IQuadrature quadrature;
    private IElementFactory elementFactory;
}

public static class ProgramRunner
{
    public static void Main(string[] args)
    {
        var setup = ConsoleInput.StartInteractiveSession();

        // Создание пространственной сетки
        ISpaceMesh spaceMesh = new UniformSpaceMesh(0.0, 1.0, setup.Nodes.Count);
        
        // Создание временной сетки
        ITimeMesh timeMesh = new UniformTimeMesh(0.0, 1.0, dt: 0.01);

        Console.WriteLine($"Пространственная сетка: {spaceMesh.NodeCount} узлов");
        Console.WriteLine($"Временная сетка: {timeMesh.TimePoints.Count} точек");

        // Начальное условие
        IInitialCondition initialCondition = new FunctionInitialCondition(x => Math.Sin(Math.PI * x));

        // Параболическая задача
        var problem = new ParabolicProblem(spaceMesh, setup.Coefficients, setup.BoundaryConditions);

        // Решатели
        var linearSolver = new BandedLUSolver();
        var nonlinearSolver = new SimpleIterationSolver(linearSolver, tolerance: 1e-6, maxIterations: 100, relaxation: 1.0);
        var timeIntegrator = new ImplicitEulerTimeIntegrator(nonlinearSolver);

        // Получение начального вектора
        var u0 = new double[problem.Size];
        for (int i = 0; i < problem.Size; i++)
        {
            u0[i] = initialCondition.Value(spaceMesh.GetNodeCoordinate(i));
        }

        // Решение параболической задачи
        Console.WriteLine("Решение параболической задачи...");
        var sw = Stopwatch.StartNew();
        
        var solutions = new List<double[]>();
        var timePoints = new List<double>();
        var uCurrent = u0;
        
        solutions.Add((double[])uCurrent.Clone());
        timePoints.Add(timeMesh.TimeStart);
        
        int stepCount = 0;
        for (int step = 0; step < timeMesh.TimeStepCount; step++)
        {
            double t = timeMesh.GetTimePoint(step);
            double dt = timeMesh.GetTimeStep(step);
            
            uCurrent = timeIntegrator.TimeStep(problem, uCurrent, dt, t);
            solutions.Add((double[])uCurrent.Clone());
            timePoints.Add(t + dt);
            stepCount++;
            
            if (stepCount % 10 == 0)
                Console.WriteLine($"Шаг {stepCount}, t = {t + dt:F4}");
        }

        sw.Stop();

        // Сохранение результатов
        var jsonOptions = new JsonSerializerOptions { WriteIndented = true };
        var output = new
        {
            TimePoints = timePoints.ToArray(),
            Solutions = solutions.Select(s => s.Take(Math.Min(10, s.Length)).ToArray()).ToArray(),
            TotalTimeSteps = stepCount,
            ComputationTime = sw.Elapsed.TotalSeconds
        };
        
        File.WriteAllText("parabolic_results.json", JsonSerializer.Serialize(output, jsonOptions));
        
        Console.WriteLine($"Задача решена за {sw.Elapsed.TotalSeconds:F2} сек");
        Console.WriteLine($"Выполнено временных шагов: {stepCount}");
    }
}

// Равномерная пространственная сетка
public class UniformSpaceMesh : ISpaceMesh
{
    private readonly List<double> nodes;
    private readonly List<IFiniteElement> elements;
    private readonly double leftBoundary;
    private readonly double rightBoundary;

    public IReadOnlyList<double> Nodes => nodes.AsReadOnly();
    public IReadOnlyList<IFiniteElement> Elements => elements.AsReadOnly();
    public int NodeCount => nodes.Count;
    public int ElementCount => elements.Count;
    public double LeftBoundary => leftBoundary;
    public double RightBoundary => rightBoundary;

    // Создать равномерную пространственную сетку:
    // x₀ = leftBoundary, x₁ = leftBoundary + h, ..., x_{N-1} = rightBoundary
    public UniformSpaceMesh(double leftBoundary, double rightBoundary, int nodeCount)
    {
        if (nodeCount < 2) 
            throw new ArgumentException("nodeCount must be >= 2", nameof(nodeCount));
        if (rightBoundary <= leftBoundary) 
            throw new ArgumentException("rightBoundary must be > leftBoundary", nameof(rightBoundary));

        this.leftBoundary = leftBoundary;
        this.rightBoundary = rightBoundary;
        this.nodes = new List<double>();
        this.elements = new List<IFiniteElement>();

        // Создание узлов
        double step = (rightBoundary - leftBoundary) / (nodeCount - 1);
        for (int i = 0; i < nodeCount; i++)
        {
            nodes.Add(leftBoundary + i * step);
        }

        // Создание линейных элементов
        for (int i = 0; i < nodeCount - 1; i++)
        {
            var element = new FiniteElement1D(i, [i, i + 1], nodes[i], nodes[i + 1]);
            elements.Add(element);
        }
    }

    public double GetNodeCoordinate(int nodeId)
    {
        if (nodeId < 0 || nodeId >= nodes.Count)
            throw new ArgumentOutOfRangeException(nameof(nodeId));
        return nodes[nodeId];
    }

    public IFiniteElement GetElement(int elementId)
    {
        if (elementId < 0 || elementId >= elements.Count)
            throw new ArgumentOutOfRangeException(nameof(elementId));
        return elements[elementId];
    }

    public override string ToString()
    {
        return $"UniformSpaceMesh: [{LeftBoundary}, {RightBoundary}], Nodes={NodeCount}, Elements={ElementCount}";
    }
}

// Равномерная временная сетка
public class UniformTimeMesh : ITimeMesh
{
    private readonly List<double> timePoints;
    private readonly double timeStart;
    private readonly double timeEnd;
    private readonly double dt;

    public IReadOnlyList<double> TimePoints => timePoints.AsReadOnly();
    public int NodeCount => timePoints.Count;
    public int TimeStepCount => timePoints.Count - 1;
    public double TimeStart => timeStart;
    public double TimeEnd => timeEnd;

    // Создать равномерную временную сетку:
    // t₀ = timeStart, t₁ = timeStart + dt, ..., t_{M} = timeEnd
    public UniformTimeMesh(double timeStart, double timeEnd, double dt)
    {
        if (timeEnd <= timeStart) 
            throw new ArgumentException("timeEnd must be > timeStart", nameof(timeEnd));
        if (dt <= 0) 
            throw new ArgumentException("dt must be positive", nameof(dt));

        this.timeStart = timeStart;
        this.timeEnd = timeEnd;
        this.dt = dt;
        this.timePoints = new List<double>();

        // Количество шагов
        int nSteps = (int)Math.Ceiling((timeEnd - timeStart) / dt);

        // Создание точек
        for (int n = 0; n <= nSteps; n++)
        {
            double t = timeStart + n * dt;
            if (t > timeEnd + 1e-10) break;
            timePoints.Add(t);
        }

        // Убедиться, что последняя точка близка к timeEnd
        if (Math.Abs(timePoints[timePoints.Count - 1] - timeEnd) > 1e-10)
        {
            timePoints[timePoints.Count - 1] = timeEnd;
        }
    }

    public double GetTimeStep(int stepIndex)
    {
        if (stepIndex < 0 || stepIndex >= TimeStepCount)
            throw new ArgumentOutOfRangeException(nameof(stepIndex));
        return dt;
    }

    public double GetTimePoint(int stepIndex)
    {
        if (stepIndex < 0 || stepIndex >= timePoints.Count)
            throw new ArgumentOutOfRangeException(nameof(stepIndex));
        return timePoints[stepIndex];
    }

    public override string ToString()
    {
        return $"UniformTimeMesh: [{TimeStart}, {TimeEnd}], dt={dt}, TimePoints={NodeCount}, TimeSteps={TimeStepCount}";
    }
}

/// <summary>
/// Полная параболическая задача с дискретизацией в пространстве и времени
/// Используется вместе с ITimeIntegrator для решения нестационарных задач
/// </summary>
public class ParabolicProblem : IParabolicProblem
{
    private readonly ISpaceMesh spaceMesh;
    private readonly ICoefficients coeffs;
    private readonly IBoundaryConditions bc;

    public int Size => spaceMesh.NodeCount;

    public ParabolicProblem(ISpaceMesh spaceMesh, ICoefficients coeffs, IBoundaryConditions bc)
    {
        this.spaceMesh = spaceMesh ?? throw new ArgumentNullException(nameof(spaceMesh));
        this.coeffs = coeffs ?? throw new ArgumentNullException(nameof(coeffs));
        this.bc = bc ?? throw new ArgumentNullException(nameof(bc));
    }

    public virtual double GetMassCoefficient(double u, double x, double t)
    {
        return 1.0;  // Коэффициент при ∂u/∂t
    }

    /// <summary>
    /// Вычисляет систему для одного временного шага неявного метода Эйлера:
    /// [M/dt + A(u^(n+1))]·u^(n+1) = M·u^n/dt + b(u^(n+1))
    /// где M — матрица масс, A — матрица жесткости, u^n — решение на предыдущем слое
    /// </summary>
    public void ComputeParabolicSystem(IReadOnlyList<double> uPrev, double dt, double time,
        out BandedMatrix A, out double[] b)
    {
        if (uPrev == null) throw new ArgumentNullException(nameof(uPrev));
        if (dt <= 0) throw new ArgumentException("dt must be positive", nameof(dt));

        int n = spaceMesh.NodeCount;
        A = new BandedMatrix(n, 1);
        b = new double[n];

        var x = new double[n];
        for (int i = 0; i < n; i++) x[i] = spaceMesh.GetNodeCoordinate(i);

        double dtInv = 1.0 / dt;

        // Внутренние узлы
        for (int i = 1; i < n - 1; i++)
        {
            double hL = x[i] - x[i - 1];
            double hR = x[i + 1] - x[i];
            if (hL <= 0 || hR <= 0)
                throw new InvalidOperationException("Mesh nodes must be strictly increasing");

            double dudxL = (uPrev[i] - uPrev[i - 1]) / hL;
            double dudxR = (uPrev[i + 1] - uPrev[i]) / hR;
            double xLmid = 0.5 * (x[i] + x[i - 1]);
            double xRmid = 0.5 * (x[i] + x[i + 1]);
            double uLmid = 0.5 * (uPrev[i] + uPrev[i - 1]);
            double uRmid = 0.5 * (uPrev[i] + uPrev[i + 1]);

            double sigmaL = coeffs.Sigma(uLmid, xLmid, dudxL, time);
            double sigmaR = coeffs.Sigma(uRmid, xRmid, dudxR, time);

            double aL = sigmaL / hL;
            double aR = sigmaR / hR;

            double massCoeff = GetMassCoefficient(uPrev[i], x[i], time);

            A[i, i - 1] = -aL;
            A[i, i] = aL + aR + coeffs.Lambda(uPrev[i], x[i], time) + massCoeff * dtInv;
            A[i, i + 1] = -aR;
            
            b[i] = massCoeff * uPrev[i] * dtInv + coeffs.F(uPrev[i], x[i], time);
        }

        ApplyBoundaryConditionsParabolic(uPrev, dt, time, A, b, dtInv);
    }

    private void ApplyBoundaryConditionsParabolic(IReadOnlyList<double> uPrev, double dt, double time,
        BandedMatrix A, double[] b, double dtInv)
    {
        int n = spaceMesh.NodeCount;

        // Левое граничное условие
        if (bc.LeftType == BoundaryType.Dirichlet)
        {
            for (int j = 0; j < n; j++) A[0, j] = 0.0;
            A[0, 0] = 1.0;
            b[0] = bc.LeftValue;
        }
        else if (bc.LeftType == BoundaryType.Neumann)
        {
            double x0 = spaceMesh.GetNodeCoordinate(0);
            double x1 = spaceMesh.GetNodeCoordinate(1);
            double hR = x1 - x0;
            double massCoeff = GetMassCoefficient(uPrev[0], x0, time);

            A[0, 0] = coeffs.Lambda(uPrev[0], x0, time) + massCoeff * dtInv;
            A[0, 1] = -coeffs.Sigma(0.5 * (uPrev[0] + uPrev[1]), 0.5 * (x0 + x1), 
                                    (uPrev[1] - uPrev[0]) / hR, time) / hR;
            b[0] = massCoeff * uPrev[0] * dtInv + coeffs.F(uPrev[0], x0, time) - bc.LeftValue;
        }

        // Правое граничное условие
        int lastNode = n - 1;
        if (bc.RightType == BoundaryType.Dirichlet)
        {
            for (int j = 0; j < n; j++) A[lastNode, j] = 0.0;
            A[lastNode, lastNode] = 1.0;
            b[lastNode] = bc.RightValue;
        }
        else if (bc.RightType == BoundaryType.Neumann)
        {
            double xN = spaceMesh.GetNodeCoordinate(lastNode);
            double xN1 = spaceMesh.GetNodeCoordinate(lastNode - 1);
            double hL = xN - xN1;
            double massCoeff = GetMassCoefficient(uPrev[lastNode], xN, time);

            A[lastNode, lastNode] = coeffs.Lambda(uPrev[lastNode], xN, time) + massCoeff * dtInv;
            A[lastNode, lastNode - 1] = -coeffs.Sigma(0.5 * (uPrev[lastNode] + uPrev[lastNode - 1]),
                                                       0.5 * (xN + xN1),
                                                       (uPrev[lastNode] - uPrev[lastNode - 1]) / hL, time) / hL;
            b[lastNode] = massCoeff * uPrev[lastNode] * dtInv + coeffs.F(uPrev[lastNode], xN, time) + bc.RightValue;
        }
    }

    /// <summary>
    /// Для интеграции с методом простой итерации (решение стационарной части)
    /// </summary>
    public void ComputeSystem(IReadOnlyList<double> q, out BandedMatrix A, out double[] b)
    {
        throw new NotImplementedException(
            "Use ComputeParabolicSystem for parabolic problems. " +
            "For elliptic problems, use NonlinearSpatialProblem instead.");
    }

    public double[] ComputeResidual(IReadOnlyList<double> q)
    {
        throw new NotImplementedException("Use ComputeParabolicSystem for parabolic problems.");
    }

    public double ComputeResidualNorm(IReadOnlyList<double> q)
    {
        throw new NotImplementedException("Use ComputeParabolicSystem for parabolic problems.");
    }
}
