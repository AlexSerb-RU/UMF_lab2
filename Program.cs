using System;
using System.Collections.Generic;

public interface IBasisFunction
{
    double Value(double x);
    double Derivative(double x);
}

public interface IFiniteElement
{
    int Id { get; }
    IReadOnlyList<int> GlobalNodes { get; }
    IReadOnlyList<IBasisFunction> BasisFunctions { get; }
    double Length { get; }
}

public interface IMesh
{
    IReadOnlyList<double> Nodes { get; }
    IReadOnlyList<IFiniteElement> Elements { get; }
    int NodeCount { get; }
    int ElementCount { get; }
    double GetNodeCoordinate(int nodeId);
    IFiniteElement GetElement(int elementId);
}

public interface ICoefficients
{
    double Lambda(double u, double x);

    double Gamma(double u, double x);

    double F(double u, double x);

    // для краевого условия 2-го рода
    double Theta(double u, double x);
}

public interface IBoundaryConditions
{
    // Тип условий
    BoundaryType LeftType { get; }

    BoundaryType RightType { get; }

    // значения на границе (для 1-го рода)
    double LeftValue { get; }

    double RightValue { get; }

    // Функции для граничных условий
    Func<double, double> LeftBoundaryFunction { get; }

    Func<double, double> RightBoundaryFunction { get; }

    bool IsBoundaryNode(int nodeId, int totalNodes);
}

public enum BoundaryType
{
    Dirichlet,   // 1-го рода
    Neumann,     // 2-го рода
    Robin        // 3-го рода
}

public interface ILocalAssembly
{
    double[,] LocalMatrix { get; }
    double[] LocalVector { get; }
    int Size { get; }

    void Compute(IFiniteElement element, ICoefficients coeffs,
                    IReadOnlyList<double> globalSolution);
}

public interface IGlobalAssembly
{
    BandedMatrix GlobalMatrix { get; }
    double[] GlobalVector { get; }

    void Assemble(IMesh mesh, ILocalAssembly localAssembly,
        ICoefficients coeffs, IReadOnlyList<double> solution,
        IBoundaryConditions bc);

    void ApplyBoundaryConditions(IBoundaryConditions bc, IReadOnlyList<double> solution);
}

public interface ILinearSolver
{
    double[] Solve(BandedMatrix matrix, double[] rhs);
}

public interface INonlinearSolver
{
    NonlinearSolverResult Solve(INonlinearProblem problem, double[] initialGuess);
}

public interface INonlinearProblem
{
    void ComputeSystem(IReadOnlyList<double> q, out BandedMatrix A, out double[] b);

    double[] ComputeResidual(IReadOnlyList<double> q);

    double ComputeResidualNorm(IReadOnlyList<double> q);

    int Size { get; }
}

public interface ILinearizationStrategy
{
    void Linearize(ILocalAssembly localAssembly, IFiniteElement element,
        ICoefficients coeffs, IReadOnlyList<double> solution,
        out double[,] linearizedMatrix, out double[] linearizedVector);
}

public interface IQuadrature
{
    IReadOnlyList<double> Points { get; }
    IReadOnlyList<double> Weights { get; }
    int Order { get; }
}

public interface IElementFactory
{
    IFiniteElement CreateElement(int id, int node1, int node2);
    IBasisFunction CreateLinearBasisFunction(int localId, double x1, double x2);
    IBasisFunction CreateQuadraticBasisFunction(int localId, double x1, double x2, double x3);
}

public interface ISolutionResult
{
    double[] Solution { get; }
    double[] ResidualHistory { get; }
    int IterationsCount { get; }
    bool Converged { get; }
    double FinalResidual { get; }
    double ComputationTime { get; }
}

public class LinearBasisFunction : IBasisFunction
{
    public LinearBasisFunction(int localId, double leftX, double rightX)
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

public class FiniteElement : IFiniteElement
{
    public FiniteElement(int id, IReadOnlyList<int> globalNodes,
                            IReadOnlyList<IBasisFunction> basisFunctions)
    {
        // Инициализация
    }

    public int Id { get; private set; }
    public IReadOnlyList<int> GlobalNodes { get; private set; }
    public IReadOnlyList<IBasisFunction> BasisFunctions { get; private set; }
    public double Length { get; private set; }
}

public class Mesh : IMesh
{
    public Mesh(IReadOnlyList<double> nodes, IElementFactory elementFactory, int basisOrder)
    {
        // Инициализация
    }

    public IReadOnlyList<double> Nodes { get; private set; }
    public IReadOnlyList<IFiniteElement> Elements { get; private set; }
    public int NodeCount { get; private set; }
    public int ElementCount { get; private set; }

    public double GetNodeCoordinate(int nodeId)
    {
        // Реализация
        throw new NotImplementedException();
    }

    public IFiniteElement GetElement(int elementId)
    {
        // Реализация
        throw new NotImplementedException();
    }
}

public class GalerkinLocalAssembly : ILocalAssembly
{
    public GalerkinLocalAssembly(IQuadrature quadrature)
    {
        // Инициализация
    }

    public double[,] LocalMatrix { get; private set; }
    public double[] LocalVector { get; private set; }
    public int Size { get; private set; }

    public void Compute(IFiniteElement element, ICoefficients coeffs,
                        IReadOnlyList<double> globalSolution)
    {
        // вычисление интегралов
        throw new NotImplementedException();
    }
}

public class GlobalAssembly : IGlobalAssembly
{
    public GlobalAssembly(int nodeCount, int bandwidth)
    {
        // Инициализация
    }

    public BandedMatrix GlobalMatrix { get; private set; }
    public double[] GlobalVector { get; private set; }

    public void Assemble(IMesh mesh, ILocalAssembly localAssembly,
                            ICoefficients coeffs, IReadOnlyList<double> solution,
                            IBoundaryConditions bc)
    {
        // разборка локальных матриц и векторов
        throw new NotImplementedException();
    }

    public void ApplyBoundaryConditions(IBoundaryConditions bc, IReadOnlyList<double> solution)
    {
        // Реализация: модификация матрицы и вектора с учётом КУ
        throw new NotImplementedException();
    }
}

public class BandedMatrix
{
    public BandedMatrix(int size, int bandwidth)
    {
        // Инициализация
    }

    public double this[int i, int j]
    {
        get
        {
            // Реализация
            throw new NotImplementedException();
        }
        set
        {
            // Реализация
            throw new NotImplementedException();
        }
    }

    public int Size { get; private set; }
    public int Bandwidth { get; private set; }

    public void Clear()
    {
        // Реализация
        throw new NotImplementedException();
    }
}

public class BandedLUSolver : ILinearSolver
{
    public double[] Solve(BandedMatrix matrix, double[] rhs)
    {
        // Реализация: LU-разложение для ленточной матрицы
        throw new NotImplementedException();
    }
}

public class SimpleIterationSolver : INonlinearSolver
{
    public SimpleIterationSolver(ILinearSolver linearSolver, double tolerance,
                                    int maxIterations, double relaxation = 1.0)
    {
        // Инициализация
    }

    public NonlinearSolverResult Solve(INonlinearProblem problem, double[] initialGuess)
    {
        // Реализация: q(k+1) = (A(q(k)))^(-1) * b(q(k))
        // с возможностью релаксации
        throw new NotImplementedException();
    }

    public double RelaxationParameter { get; set; }
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

public class NonlinearEllipticProblem : INonlinearProblem
{
    public NonlinearEllipticProblem(IMesh mesh, ICoefficients coeffs,
        IBoundaryConditions bc, IGlobalAssembly assembler,
        int bandwidth)
    {
        // Инициализация
    }

    public int Size { get; private set; }

    public void ComputeSystem(IReadOnlyList<double> q, out BandedMatrix A, out double[] b)
    {
        // сборка A(q) и b(q)
        throw new NotImplementedException();
    }

    public double[] ComputeResidual(IReadOnlyList<double> q)
    {
        // R(q) = A(q)·q - b(q)
        throw new NotImplementedException();
    }

    public double ComputeResidualNorm(IReadOnlyList<double> q)
    {
        // ||R(q)|| / ||b(q)||
        throw new NotImplementedException();
    }

    protected IMesh Mesh { get; set; }
    protected ICoefficients Coefficients { get; set; }
    protected IBoundaryConditions BoundaryConditions { get; set; }
    protected IGlobalAssembly Assembler { get; set; }
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
    public GaussQuadrature(int order)
    {
        // Инициализация узлов и весов Гаусса
    }

    public IReadOnlyList<double> Points { get; private set; }
    public IReadOnlyList<double> Weights { get; private set; }
    public int Order { get; private set; }
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
        // Инициализация
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